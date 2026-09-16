#!/usr/bin/env python3
"""zBET: Fast Blade Element Theory (BET) solver for rapid rotor aerodynamic evaluation,
conceptual trade studies, and advance ratio (mu) sweeps.

Key Model Characteristics:
  - Methodology: Fast Blade Element Theory (BET) with closed-form radial-moment integrals,
    coupled with global actuator-disk momentum theory for mean induced inflow (lambda_i),
    and classical spatial inflow gradient models (Glauert, Coleman, Drees).
    (Note: This is a fast analytical/semi-empirical BET solver designed for rapid sizing,
    NOT an iterative discretized multi-annulus BEMT strip-theory solver).

Coordinate System & Sign Conventions:
  - Hub Coordinate Frame (standard helicopter body/shaft convention):
      * x-axis: Points FORWARD (along nominal vehicle flight direction / aircraft nose).
      * y-axis: Points to the RIGHT (starboard side).
      * z-axis: Points DOWNWARD (through the bottom of the rotor disk).
  - Rotor Rotation:
      * Viewed from ABOVE (looking down along +z): The rotor blades rotate COUNTER-CLOCKWISE (CCW).
      * Advancing blade: Starboard / right side (azimuth psi = 90 deg), where tangential velocity is u_T = x + mu*sin(psi).
      * Retreating blade: Port / left side (azimuth psi = 270 deg), where tangential velocity is u_T = x - mu.
  - Shaft Torque & Drive Direction:
      * By the right-hand rule about the downward +z axis, CCW rotor rotation corresponds to a vector along -z.
      * Aerodynamic blade drag resists rotation in the CLOCKWISE direction (+z, viewed from above).
      * Drive shaft torque (Q) overcomes blade drag, delivering driving power in the -z direction.
      * By convention, rotor shaft torque CQ is positive (CQ > 0) when power is delivered to the rotor (P = Q*Omega > 0).
      * The reaction torque exerted by the rotor on the fuselage acts in the CLOCKWISE (+z) direction.
  - Axial Flow & mu_z (Wind Direction):
      * Inflow velocity lambda = mu_z + lambda_i points DOWNWARD (+z direction).
      * Mean downwash lambda_i >= 0 is always directed downward (+z).
      * mu_z > 0: Relative oncoming wind is flowing DOWNWARD (wind coming from ABOVE the rotor disk,
        e.g., helicopter in vertical climb or top-down relative airflow).
      * mu_z < 0: Relative oncoming wind is flowing UPWARD (wind coming from BELOW the rotor disk,
        e.g., forward flight with forward rotor tilt alpha > 0, where mu_z = -mu*tan(alpha) < 0,
        or helicopter in vertical descent / autorotation).

Theoretical and empirical references:
  1. Glauert, H. (1926) - Classic uniform inflow (Kx = 0, Ky = 0).
  2. Coleman, R. P., Feingold, A. M., Stempin, C. W. (1945) - NACA ARR L5E10 / WR L-101:
     Longitudinal inflow gradient as a function of wake skew angle chi:
     Kx = tan(chi/2) = mu / (sqrt(mu^2 + lambda^2) + |lambda|), Ky = 0.
  3. Johnson, Wayne - NASA/TP-2009-215402 (NDARC Theory) and Rotorcraft Aeromechanics (2013):
     Coleman-Feingold inflow with 15*pi/32 first harmonic and lateral gradient Ky = -fy * 2*mu:
     Kx = fx * (15*pi/32) * tan(chi/2), Ky = -fy * 2*mu.
  4. Drees, J. M. (1949) - A Theory of Airflow Through Rotors:
     Longitudinal and lateral inflow gradients:
     Kx = (4/3) * (1 - 1.8*mu^2) * tan(chi/2), Ky = -2*mu.
  5. Leishman, J. G. - Principles of Helicopter Aerodynamics (2nd ed., 2006):
     Linear induced velocity modeling, Figure of Merit (FoM), and Blade Element Theory moments.
  6. Johnson, Wayne - Helicopter Theory (1980) and Rotorcraft Aeromechanics (2013):
     Rotor effective lift-to-drag ratio: (L/D)_eff = mu * CT / (mu * CH + CQ) = mu * CT / CPair.
"""

from dataclasses import dataclass
from functools import lru_cache
import math
from pathlib import Path
import sys

# Ensure UTF-8 output on Windows for printing symbols (alpha, mu, deg, etc.)
if sys.stdout and hasattr(sys.stdout, "reconfigure"):
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# =============================================================================
# USER CONFIGURATION — Edit this block to set up a new rotor case
# =============================================================================

# -----------------------------------------------------------------------------
# 1. ATMOSPHERE & AMBIENT CONDITIONS
# -----------------------------------------------------------------------------
RHO = 1.225               # Air density [kg/m^3] (ISA standard sea level)
SPEED_OF_SOUND = 340.3    # Speed of sound in air [m/s]

# -----------------------------------------------------------------------------
# 2. ROTOR GEOMETRY & BASIC AERODYNAMICS (Standard Utility Helicopter)
# -----------------------------------------------------------------------------
RPM = 390.0              # Rotor rotational speed [RPM]
R = 5.5                  # Rotor blade outer radius R [m]
R0_BAR = 0.15            # Blade root cutout ratio x0 = r0 / R [-] (first aerodynamic station)
A_LIFT = 5.73            # Incompressible 2D section lift curve slope 'a0' [1/rad]
CD0 = 0.009              # Profile parasite drag coefficient [-] (typical NACA 0012)

# -----------------------------------------------------------------------------
# 3. SOLIDITY & CHORD CONFIGURATION (Choose Option 1, Option 2, or Option 3)
#   "sigma_ref"  -> Option 1: Provide reference solidity extrapolated to hub (r=0)
#   "sigma_geom" -> Option 2: Provide true geometric solidity (actual blade area / disk area)
#   "chords"     -> Option 3: Provide blade count N, root chord at x0, and tip chord at x=1.0
# -----------------------------------------------------------------------------
SOLIDITY_MODE = "sigma_ref"   # "sigma_ref", "sigma_geom", or "chords" (taper)

# Option 1 parameter (SOLIDITY_MODE = "sigma_ref"):
SIGMA_REF = 0.075        # Reference solidity extrapolated to hub [-]

# Option 2 parameter (SOLIDITY_MODE = "sigma_geom"):
SIGMA_GEOM = 0.06375     # True geometric solidity [-] (actual physical blade area / disk area)

# Option 3 parameters (SOLIDITY_MODE = "chords"):
# Note: CHORD_ROOT and CHORD_TIP are defined at root cutout x0 and tip x=1.0 [m]
N_BLADES = 4             # Number of blades [-]
CHORD_ROOT = 0.324       # Chord at root cutout station x0 = R0_BAR [m]
CHORD_TIP = 0.324        # Chord at blade tip station x = 1.0 [m]

# -----------------------------------------------------------------------------
# 4. BLADE PITCH, TWIST & HOVER TRIM CONFIGURATION
#   PITCH_MODE:
#     "constant"     -> Uniform collective pitch theta0 (Option A)
#     "linear_twist" -> Linear twist defined by root and tip angles (Option B)
#
#   HOVER_TRIM_MODE (strictly evaluated at hover mu=0, then frozen for the sweep):
#     "collective"   -> Fixed RPM; adjusts collective theta (preserving total twist constant)
#     "rpm"          -> Fixed pitch/twist; adjusts RPM to generate target Thrust in Newtons
#     "none"         -> Direct input without hover trim (uses prescribed RPM and pitch directly)
# -----------------------------------------------------------------------------
PITCH_MODE = "constant"         # "constant" or "linear_twist"
HOVER_TRIM_MODE = "collective"  # "collective", "rpm", or "none"
COLLECTIVE_MODE = HOVER_TRIM_MODE  # Compatibility alias

CT_HOVER_TARGET = 0.0065        # Target hover thrust coefficient [-] (typical helicopter loading)
THRUST_HOVER_N = None           # Target hover thrust in Newtons [N] (e.g., 35000.0)
THETA0 = None                   # rad — direct pitch value when HOVER_TRIM_MODE="none"

# Twist parameters for PITCH_MODE = "linear_twist":
# Geometric pitch angles in degrees at root cutout station x0 and blade tip x=1.0.
# In "collective" trim mode, the total twist (THETA_TIP - THETA_ROOT) is kept strictly constant
# and a uniform collective offset delta_theta0 is added across the entire blade.
THETA_ROOT_DEG = 12.0    # Blade pitch at root cutout station x0 [deg]
THETA_TIP_DEG = 4.0      # Blade pitch at tip station x=1.0 [deg] (-8 deg washout)

# -----------------------------------------------------------------------------
# 5. ADVANCED AERODYNAMIC MODELS (OPTIONAL — DISABLED BY DEFAULT)
# -----------------------------------------------------------------------------
# Blade tip loss factor (B):
#   "none":     disabled (B = 1.0, classical rotor disk)
#   "fixed":    uses prescribed value from TIP_LOSS_B (e.g. 0.97)
#   "sissingh": B = 1 - sqrt(2*CT)/N (classic self-adjusting formula)
TIP_LOSS_MODE = "none"          # "none", "fixed", or "sissingh"
TIP_LOSS_B = 0.97               # Prescribed tip loss factor when TIP_LOSS_MODE="fixed"

# Prandtl-Glauert compressibility correction on lift curve slope 'a':
# False: a = A_LIFT fixed incompressible (classic baseline)
# True:  a = A_LIFT / sqrt(1 - M^2) based on effective Mach number at 0.75R in forward flight
USE_PRANDTL_GLAUERT = False     # False (default) or True

# -----------------------------------------------------------------------------
# 6. INFLOW MODELS TO EVALUATE IN ADVANCE RATIO SWEEP
# Available models:
#   "uniform":          Classic uniform inflow (Glauert: Kx = 0, Ky = 0)
#   "coleman_simple":   Simple Coleman (Coleman et al. 1945: Kx = tan(chi/2), Ky = 0)
#   "coleman_feingold": NDARC Coleman-Feingold (Kx = fx*(15*pi/32)*tan(chi/2), Ky = -fy*2*mu)
#   "drees":            Drees (1949: Kx = (4/3)*(1 - 1.8*mu^2)*tan(chi/2), Ky = -2*mu)
# -----------------------------------------------------------------------------
INFLOW_MODELS = ["uniform", "coleman_simple", "coleman_feingold", "drees"]

# Calibration scaling factors for NDARC Coleman-Feingold (official standard: fx=1.0, fy=1.0)
FX_COLEMAN = 1.0         # Longitudinal factor fx
FY_COLEMAN = 1.0         # Lateral factor fy (set to 0.0 to disable Ky)

# Semi-empirical induced shaft power factor (default 1.15)
K_IND = 1.15

# -----------------------------------------------------------------------------
# 7. OPERATING CONDITIONS (ADVANCE RATIO & AXIAL FLOW SWEEP)
# -----------------------------------------------------------------------------
MU_MIN = 0.0
MU_MAX = 0.40
MU_STEP = 0.05

# Axial flow condition: "alpha" [deg], "mu_z" [-], or "w" [m/s].
# Sign convention for axial flow through the downward +z axis:
#   - mu_z > 0 (or w > 0): Wind coming from ABOVE the rotor disk (flowing downward along +z, e.g. climb).
#   - mu_z < 0 (or w < 0): Wind coming from BELOW the rotor disk (flowing upward along -z, e.g. descent).
#   - alpha > 0: Propulsive forward disk tilt where relative wind comes from BELOW (mu_z = -mu*tan(alpha) < 0).
AXIAL_FLOW = "alpha"
AXIAL_VALUES = [0.0, -4.0, 4.0]

# Models for induced torque and profile drag:
# "complete": includes reverse flow integration, radial velocity u_R, and root cutout
# "simple_bet": classic closed-form small-angle BET equations
CQ_MODEL = "complete"         # "complete" or "simple_bet"
PROFILE_MODEL = "complete"    # "complete" or "simple_bet"

OUTPUT_DIR = Path("outputs")
# =============================================================================


OUTPUTS = [
    "CT", "CQ", "CQi", "CQ0", "CH", "CHi", "CH0", "CY",
    "CMy", "CMx", "CPair", "lambda", "lambda_i", "L_D_eff",
]
MODELS = ["uniform", "coleman_simple", "coleman_feingold", "drees"]
MODEL_LABELS = {
    "uniform": "Uniform",
    "coleman_simple": "Coleman Simple",
    "coleman_feingold": "Coleman-Feingold (NDARC)",
    "coleman": "Coleman-Feingold (NDARC)",
    "drees": "Drees",
}


@dataclass(frozen=True)
class BladeSolidity:
    """Represents rotor solidity, computing three classical aerodynamic and geometric definitions:
      - sigma_ref: Reference solidity extrapolated to hub r=0 (stretched_area / (pi*R^2)).
      - sigma_geom: True geometric solidity (actual physical blade area from x0 to 1 / (pi*R^2)).
      - sigma_thrust: Thrust-weighted equivalent solidity (3 * integral_{x0}^1 x^2 sigma(x) dx).

    Theoretical note on BET integration:
      In Blade Element Theory, differential lift dL is scaled by local chord c(r), or local
      solidity sigma(x) = N * c(x) / (pi * R). For rectangular blades, sigma(x) = sigma_ref
      throughout the active blade span. The root cutout x0 is taken into account directly by
      the integration bounds J_m = integral_{x0}^1 x^m dx.
      Therefore, sigma_ref parameterizes the analytical equations. Using sigma_geom in the
      integrals would count the root cutout penalty twice.
    """
    mode: str                  # "sigma_ref", "sigma_geom", or "chords"
    sigma_ref: float           # Reference solidity extrapolated to hub [-]
    sigma_geom: float          # True geometric physical solidity [-]
    sigma_thrust: float        # Thrust-weighted equivalent solidity [-]
    n_blades: int = 4          # Number of blades
    chord_root: float = 0.0    # Chord at root cutout station x0 [m]
    chord_tip: float = 0.0     # Chord at tip station x=1.0 [m]
    chord_center: float = 0.0  # Extrapolated chord at hub x=0 [m]
    radius: float = 1.0        # Rotor radius R [m]
    root_cutout: float = 0.0   # Root cutout ratio x0 = r0 / R

    def sigma(self, x):
        """Returns local solidity sigma(x) = N * c(x) / (pi * R).
        Accepts either float scalar or NumPy array.
        """
        x0 = self.root_cutout
        if self.chord_root == self.chord_tip:
            if isinstance(x, np.ndarray):
                return np.full_like(x, self.sigma_ref, dtype=float)
            return float(self.sigma_ref)
        c = self.chord_root + (self.chord_tip - self.chord_root) * (x - x0) / (1.0 - x0)
        return self.n_blades * c / (math.pi * self.radius)

    @property
    def linear_coeffs(self):
        """Returns (s0, s1) such that sigma(x) = s0 + s1 * x."""
        if self.chord_root == self.chord_tip:
            return self.sigma_ref, 0.0
        x0 = self.root_cutout
        sig_root = self.n_blades * self.chord_root / (math.pi * self.radius)
        sig_tip = self.n_blades * self.chord_tip / (math.pi * self.radius)
        s1 = (sig_tip - sig_root) / (1.0 - x0)
        s0 = sig_root - s1 * x0
        return s0, s1


@dataclass(frozen=True)
class BladePitch:
    """Represents blade pitch distribution, either constant collective theta0 or linearly twisted."""
    mode: str                  # "constant" or "linear_twist"
    theta0: float              # Reference or mean pitch angle [rad]
    theta_root: float          # Pitch at root cutout station x0 [rad]
    theta_tip: float           # Pitch at tip station x=1.0 [rad]
    root_cutout: float = 0.0   # x0 = r0 / R

    def theta(self, x):
        """Returns local blade pitch angle theta(x) in radians."""
        if self.mode == "constant":
            if isinstance(x, np.ndarray):
                return np.full_like(x, self.theta0, dtype=float)
            return float(self.theta0)
        x0 = self.root_cutout
        return self.theta_root + (self.theta_tip - self.theta_root) * (x - x0) / (1.0 - x0)

    @property
    def linear_coeffs(self):
        """Returns (t0, t1) such that theta(x) = t0 + t1 * x in radians."""
        if self.mode == "constant":
            return self.theta0, 0.0
        x0 = self.root_cutout
        t1 = (self.theta_tip - self.theta_root) / (1.0 - x0)
        t0 = self.theta_root - t1 * x0
        return t0, t1


class Geometry:
    """Global rotor geometry and aerodynamic properties.

    Supports both the legacy positional signature:
      Geometry(sigma, rpm, radius, lift_curve_slope, root_cutout, cd0)
    and the structured signature with BladeSolidity:
      Geometry(rpm, radius, lift_curve_slope, root_cutout, cd0, solidity)
    """
    rpm: float
    radius: float
    lift_curve_slope: float
    root_cutout: float
    cd0: float
    solidity: BladeSolidity
    rho: float
    speed_of_sound: float
    tip_loss_mode: str
    tip_loss_b: float
    use_prandtl_glauert: bool

    def __init__(self, *args, **kwargs):
        if len(args) == 6 and isinstance(args[0], (int, float)) and isinstance(args[5], (int, float)):
            # Legacy signature: (sigma, rpm, radius, lift_curve_slope, root_cutout, cd0)
            sig, rpm, rad, a_lift, r0, cd0 = args
            sol = resolve_solidity("sigma_ref", sigma_ref=sig, radius=rad, root_cutout=r0)
            self.rpm = float(rpm)
            self.radius = float(rad)
            self.lift_curve_slope = float(a_lift)
            self.root_cutout = float(r0)
            self.cd0 = float(cd0)
            self.solidity = sol
        elif len(args) == 6 and isinstance(args[5], BladeSolidity):
            rpm, rad, a_lift, r0, cd0, sol = args
            self.rpm = float(rpm)
            self.radius = float(rad)
            self.lift_curve_slope = float(a_lift)
            self.root_cutout = float(r0)
            self.cd0 = float(cd0)
            self.solidity = sol
        else:
            self.rpm = float(kwargs.get("rpm", args[0] if len(args) > 0 else RPM))
            self.radius = float(kwargs.get("radius", args[1] if len(args) > 1 else R))
            self.lift_curve_slope = float(kwargs.get("lift_curve_slope", args[2] if len(args) > 2 else A_LIFT))
            self.root_cutout = float(kwargs.get("root_cutout", args[3] if len(args) > 3 else R0_BAR))
            self.cd0 = float(kwargs.get("cd0", args[4] if len(args) > 4 else CD0))
            self.solidity = kwargs.get("solidity", args[5] if len(args) > 5 else DEFAULT_SOLIDITY)

        self.rho = float(kwargs.get("rho", RHO))
        self.speed_of_sound = float(kwargs.get("speed_of_sound", SPEED_OF_SOUND))
        self.tip_loss_mode = str(kwargs.get("tip_loss_mode", TIP_LOSS_MODE))
        self.tip_loss_b = float(kwargs.get("tip_loss_b", TIP_LOSS_B))
        self.use_prandtl_glauert = bool(kwargs.get("use_prandtl_glauert", USE_PRANDTL_GLAUERT))

        if not 0.0 <= self.root_cutout < 1.0:
            raise ValueError("root_cutout must satisfy 0 <= root_cutout < 1")
        if min(self.rpm, self.radius, self.lift_curve_slope) <= 0.0:
            raise ValueError("rpm, radius, and lift_curve_slope must be positive")
        if self.cd0 < 0.0:
            raise ValueError("cd0 must be non-negative")

    @property
    def sigma(self):
        """Reference solidity of the rotor."""
        return self.solidity.sigma_ref

    @property
    def omega(self):
        """Rotor angular velocity [rad/s]."""
        return self.rpm * 2.0 * math.pi / 60.0

    @property
    def vtip(self):
        """Blade tip rotational speed [m/s]."""
        return self.omega * self.radius

    @property
    def disk_area(self):
        """Rotor disk area A = pi * R^2 [m^2]."""
        return math.pi * (self.radius ** 2)

    @property
    def tip_mach(self):
        """Hover blade tip Mach number (Vtip / a_sound)."""
        return self.vtip / self.speed_of_sound if self.speed_of_sound > 0.0 else 0.0

    def with_rpm(self, new_rpm):
        """Creates a copy of geometry with an updated RPM (used in RPM hover trim)."""
        return Geometry(
            new_rpm,
            self.radius,
            self.lift_curve_slope,
            self.root_cutout,
            self.cd0,
            self.solidity,
            rho=self.rho,
            speed_of_sound=self.speed_of_sound,
            tip_loss_mode=self.tip_loss_mode,
            tip_loss_b=self.tip_loss_b,
            use_prandtl_glauert=self.use_prandtl_glauert,
        )

    def b_factor(self, ct=None):
        """Returns the blade tip loss factor B in (x0, 1.0]."""
        if self.tip_loss_mode == "none":
            return 1.0
        if self.tip_loss_mode == "fixed":
            return max(self.root_cutout + 0.01, min(1.0, float(self.tip_loss_b)))
        if self.tip_loss_mode == "sissingh":
            if ct is not None and ct > 0.0:
                nb = self.solidity.n_blades
                b_val = 1.0 - math.sqrt(2.0 * float(ct)) / nb
                return max(self.root_cutout + 0.01, min(1.0, b_val))
            return 1.0
        return 1.0

    def lift_slope(self, mu=0.0):
        """Returns 2D section lift curve slope 'a', corrected for Prandtl-Glauert compressibility if enabled."""
        if not self.use_prandtl_glauert:
            return self.lift_curve_slope
        # Representative Mach number at 75% radius in forward flight:
        m_eff = self.tip_mach * math.sqrt(0.75 * 0.75 + 0.5 * mu * mu)
        m_eff = min(m_eff, 0.85)  # Subsonic safety cap
        return self.lift_curve_slope / math.sqrt(max(0.01, 1.0 - m_eff * m_eff))

    def thrust_from_ct(self, ct):
        """Calculates dimensional thrust in Newtons: T = CT * rho * A * Vtip^2."""
        return float(ct) * self.rho * self.disk_area * (self.vtip ** 2)

    def ct_from_thrust(self, thrust_n):
        """Calculates non-dimensional CT from thrust in Newtons: CT = T / (rho * A * Vtip^2)."""
        denom = self.rho * self.disk_area * (self.vtip ** 2)
        if denom <= 0.0:
            raise ValueError("Thrust non-dimensionalization denominator must be positive")
        return float(thrust_n) / denom


def resolve_solidity(
    mode=SOLIDITY_MODE,
    sigma_ref=SIGMA_REF,
    sigma_geom=SIGMA_GEOM,
    n_blades=N_BLADES,
    chord_root=CHORD_ROOT,
    chord_tip=CHORD_TIP,
    radius=R,
    root_cutout=R0_BAR,
):
    """Constructs the BladeSolidity object, analytically calculating the three solidezes:
      1. sigma_ref: Reference solidity extrapolated to hub (stretched area from 0 to 1 / (pi*R^2)).
      2. sigma_geom: True geometric physical solidity (actual area from x0 to 1 / (pi*R^2)).
      3. sigma_thrust: Thrust-weighted solidity (3 * integral_{x0}^1 x^2 sigma(x) dx).
    """
    x0 = float(root_cutout)
    rad = float(radius)
    nb = int(n_blades)

    # Case 1: Input by reference solidity (constant chord)
    # Guarantees: output.sigma_ref == input sigma_ref exactly
    if mode in ("sigma_ref", "reference"):
        s_ref = float(sigma_ref)
        if s_ref <= 0.0:
            raise ValueError("SIGMA_REF must be positive")
        c = s_ref * math.pi * rad / nb
        s_geom = (1.0 - x0) * s_ref
        s_thrust = (1.0 - x0 ** 3) * s_ref
        return BladeSolidity(
            mode="sigma_ref",
            sigma_ref=s_ref,          # Identical to input SIGMA_REF
            sigma_geom=s_geom,        # s_geom = (1 - x0) * s_ref
            sigma_thrust=s_thrust,    # s_thrust = (1 - x0^3) * s_ref
            n_blades=nb,
            chord_root=c,
            chord_tip=c,
            chord_center=c,
            radius=rad,
            root_cutout=x0,
        )

    # Case 2: Input by true geometric physical solidity (constant chord)
    # Guarantees: output.sigma_geom == input sigma_geom exactly
    if mode == "sigma_geom":
        s_geom = float(sigma_geom)
        if s_geom <= 0.0:
            raise ValueError("SIGMA_GEOM must be positive")
        s_ref = s_geom / (1.0 - x0)
        c = s_ref * math.pi * rad / nb
        s_thrust = (1.0 - x0 ** 3) * s_ref
        return BladeSolidity(
            mode="sigma_geom",
            sigma_ref=s_ref,          # s_ref = s_geom / (1 - x0)
            sigma_geom=s_geom,        # Identical to input SIGMA_GEOM
            sigma_thrust=s_thrust,
            n_blades=nb,
            chord_root=c,
            chord_tip=c,
            chord_center=c,
            radius=rad,
            root_cutout=x0,
        )

    # Case 3: Input by root chord (at x0) and tip chord (at x=1.0)
    if mode in ("chords", "taper"):
        c_root = float(chord_root)
        c_tip = float(chord_tip)
        if min(nb, c_root, c_tip, rad) <= 0:
            raise ValueError("N_BLADES, CHORD_ROOT, CHORD_TIP, and R must be positive")

        # Linear chord distribution c(x) and extrapolation to center (x=0)
        c_center = c_root - (c_tip - c_root) * x0 / (1.0 - x0)
        sig_root = nb * c_root / (math.pi * rad)
        sig_tip = nb * c_tip / (math.pi * rad)
        s1 = (sig_tip - sig_root) / (1.0 - x0)
        s0 = sig_root - s1 * x0

        # 1. Reference solidity extrapolated to hub: integral from 0 to 1 of sigma(x) dx
        s_ref = s0 + 0.5 * s1
        # 2. True geometric solidity: integral from x0 to 1 of sigma(x) dx
        s_geom = (1.0 - x0) * (s0 + 0.5 * s1 * (1.0 + x0))
        # 3. Thrust-weighted solidity: 3 * integral_{x0}^1 x^2 sigma(x) dx
        s_thrust = 3.0 * (s0 * (1.0 - x0 ** 3) / 3.0 + s1 * (1.0 - x0 ** 4) / 4.0)

        return BladeSolidity(
            mode="chords",
            sigma_ref=s_ref,
            sigma_geom=s_geom,
            sigma_thrust=s_thrust,
            n_blades=nb,
            chord_root=c_root,
            chord_tip=c_tip,
            chord_center=c_center,
            radius=rad,
            root_cutout=x0,
        )

    raise ValueError(
        f"Invalid SOLIDITY_MODE: '{mode}'. Use 'sigma_ref', 'sigma_geom', or 'chords'."
    )


DEFAULT_SOLIDITY = resolve_solidity()
DEFAULT_GEOMETRY = Geometry(RPM, R, A_LIFT, R0_BAR, CD0, DEFAULT_SOLIDITY)


def radial_integrals(geometry, b=None):
    """Returns J_n = integral(x**n dx) from x0 to B (where B is tip loss factor)."""
    r0 = geometry.root_cutout
    b_val = geometry.b_factor() if b is None else float(b)
    return {n: (b_val ** (n + 1) - r0 ** (n + 1)) / (n + 1) for n in range(8)}


def radial_moments(geometry, pitch, b=None):
    """Calculates weighted radial moments of solidity and pitch analytically from x0 to B:
      J_n = integral_{x0}^B x^n dx
      I_m = integral_{x0}^B sigma(x) * x^m dx
      T_m = integral_{x0}^B sigma(x) * theta(x) * x^m dx
    """
    r0 = geometry.root_cutout
    b_val = geometry.b_factor() if b is None else float(b)
    j = {n: (b_val ** (n + 1) - r0 ** (n + 1)) / (n + 1) for n in range(8)}

    # Linear solidity coefficients sigma(x) = s0 + s1 * x
    s0, s1 = geometry.solidity.linear_coeffs
    i_mom = {m: s0 * j[m] + s1 * j[m + 1] for m in range(6)}

    # Linear pitch coefficients theta(x) = t0 + t1 * x
    t0, t1 = pitch.linear_coeffs
    p0 = s0 * t0
    p1 = s0 * t1 + s1 * t0
    p2 = s1 * t1
    t_mom = {m: p0 * j[m] + p1 * j[m + 1] + p2 * j[m + 2] for m in range(6)}

    return j, i_mom, t_mom


def axial_condition(mu, value, mode, geometry):
    """Converts the axial flow input into non-dimensional mu_z and provides a descriptive label."""
    if mode == "alpha":
        return -mu * math.tan(math.radians(value)), f"α={value:g}°"
    if mode == "mu_z":
        return float(value), f"μ_z={value:g}"
    if mode == "w":
        return float(value) / geometry.vtip, f"w={value:g} m/s"
    raise ValueError("AXIAL_FLOW must be 'alpha', 'mu_z', or 'w'")


def trim_hover(
    geometry,
    hover_trim_mode=HOVER_TRIM_MODE,
    pitch_mode=PITCH_MODE,
    ct_hover_target=CT_HOVER_TARGET,
    thrust_hover_n=THRUST_HOVER_N,
    theta0=THETA0,
    theta_root_deg=THETA_ROOT_DEG,
    theta_tip_deg=THETA_TIP_DEG,
):
    """Performs hover trim strictly at hover (mu=0, alpha=0) and returns (geometry, pitch).

    Hover trim modes:
      1. 'collective': Constant RPM. Adjusts collective to match target CT or Thrust [N].
         For linear twist, preserves total twist delta_twist = (theta_tip - theta_root)
         STRICTLY CONSTANT and shifts the blade pitch curve by delta_theta0.
      2. 'rpm': Fixed pitch and twist. Calculates rotor CT_hover and adjusts RPM to match
         the specified target Thrust in Newtons (THRUST_HOVER_N).
      3. 'none': No hover trim; uses prescribed RPM and pitch angles directly.
    """
    x0 = geometry.root_cutout

    # Mode 3: Direct mode (no trim)
    if hover_trim_mode == "none":
        if pitch_mode == "constant":
            t0 = float(theta0) if theta0 is not None else 0.0
            pitch = BladePitch("constant", theta0=t0, theta_root=t0, theta_tip=t0, root_cutout=x0)
        elif pitch_mode == "linear_twist":
            th_root_rad = math.radians(theta_root_deg)
            th_tip_rad = math.radians(theta_tip_deg)
            mean_theta = 0.5 * (th_root_rad + th_tip_rad)
            pitch = BladePitch("linear_twist", theta0=mean_theta, theta_root=th_root_rad, theta_tip=th_tip_rad, root_cutout=x0)
        else:
            raise ValueError(f"Unknown PITCH_MODE: '{pitch_mode}'")
        return geometry, pitch

    # Mode 2: RPM trim (fixed pitch/twist, adjusts RPM to generate target thrust in Newtons)
    if hover_trim_mode == "rpm":
        if thrust_hover_n is None or float(thrust_hover_n) <= 0.0:
            raise ValueError("For HOVER_TRIM_MODE='rpm', set THRUST_HOVER_N to a positive value in Newtons.")
        t_target = float(thrust_hover_n)

        # Base fixed pitch:
        if pitch_mode == "constant":
            t0 = float(theta0) if theta0 is not None else math.radians(theta_root_deg)
            pitch = BladePitch("constant", theta0=t0, theta_root=t0, theta_tip=t0, root_cutout=x0)
        elif pitch_mode == "linear_twist":
            th_root_rad = math.radians(theta_root_deg)
            th_tip_rad = math.radians(theta_tip_deg)
            mean_theta = 0.5 * (th_root_rad + th_tip_rad)
            pitch = BladePitch("linear_twist", theta0=mean_theta, theta_root=th_root_rad, theta_tip=th_tip_rad, root_cutout=x0)
        else:
            raise ValueError(f"Unknown PITCH_MODE: '{pitch_mode}'")

        # Solve hover CT produced by this pitch:
        b_val = geometry.b_factor()
        moments = radial_moments(geometry, pitch, b=b_val)
        lift_slope = geometry.lift_slope(mu=0.0)

        def hover_residual(lambda_i):
            momentum_ct = 2.0 * (b_val * b_val) * (lambda_i ** 2)
            bet_ct = 0.5 * lift_slope * (moments[2][2] - lambda_i * moments[1][1])
            return bet_ct - momentum_ct

        # Root search for lambda_hover:
        lo, hi = 0.0, 0.5
        for _ in range(100):
            mid = 0.5 * (lo + hi)
            if hover_residual(mid) > 0.0:
                lo = mid
            else:
                hi = mid
        lam_hover = 0.5 * (lo + hi)
        ct_hover = 2.0 * (b_val * b_val) * (lam_hover ** 2)

        if ct_hover <= 1e-6:
            raise ValueError(
                f"Base pitch results in non-positive hover CT ({ct_hover:g}). "
                "Adjust pitch angles to produce positive thrust before RPM trimming."
            )

        # Required RPM: T = rho * A * Vtip^2 * CT_hover
        vtip_req = math.sqrt(t_target / (geometry.rho * geometry.disk_area * ct_hover))
        omega_req = vtip_req / geometry.radius
        rpm_req = omega_req * 60.0 / (2.0 * math.pi)
        geometry_trimmed = geometry.with_rpm(rpm_req)
        return geometry_trimmed, pitch

    # Mode 1: Collective trim (fixed RPM, adjusts pitch to match target thrust or CT)
    if hover_trim_mode == "collective":
        if thrust_hover_n is not None:
            target_ct = geometry.ct_from_thrust(thrust_hover_n)
        elif ct_hover_target is not None:
            target_ct = float(ct_hover_target)
        else:
            raise ValueError("For HOVER_TRIM_MODE='collective', set CT_HOVER_TARGET or THRUST_HOVER_N.")

        if target_ct <= 0.0:
            raise ValueError("Hover target thrust (CT or Newtons) must be positive.")

        b_val = geometry.b_factor(ct=target_ct)
        a_hover = geometry.lift_slope(mu=0.0)
        lambda_hover = math.sqrt(target_ct / (2.0 * (b_val ** 2)))
        j = radial_integrals(geometry, b=b_val)
        s0, s1 = geometry.solidity.linear_coeffs
        i_mom = {m: s0 * j[m] + s1 * j[m + 1] for m in range(5)}

        if pitch_mode == "constant":
            # CT_target = a/2 * [ theta0 * I2 - lambda_hover * I1 ]
            t0 = (2.0 * target_ct / a_hover + lambda_hover * i_mom[1]) / i_mom[2]
            pitch = BladePitch("constant", theta0=t0, theta_root=t0, theta_tip=t0, root_cutout=x0)
            return geometry, pitch

        if pitch_mode == "linear_twist":
            th_root_rad = math.radians(theta_root_deg)
            th_tip_rad = math.radians(theta_tip_deg)
            delta_twist = th_tip_rad - th_root_rad  # Total twist kept invariant!

            # theta(x) = theta_base(x) + delta_theta0
            # where theta_base(x) = th_root_rad + t1_twist * (x - x0)
            t1_twist = delta_twist / (1.0 - x0)
            t2_twist = t1_twist * (i_mom[3] - x0 * i_mom[2])

            # Collective is adjusted by adding delta_theta0 equally across entire blade:
            # T2 = (th_root_rad + delta_theta0) * I2 + t2_twist
            th_root_calc = (
                (2.0 * target_ct / a_hover + lambda_hover * i_mom[1] - t2_twist)
                / i_mom[2]
            )
            delta_theta0 = th_root_calc - th_root_rad
            th_tip_calc = th_tip_rad + delta_theta0  # Total twist (th_tip_calc - th_root_calc) is preserved!
            mean_theta = 0.5 * (th_root_calc + th_tip_calc)
            pitch = BladePitch(
                "linear_twist",
                theta0=mean_theta,
                theta_root=th_root_calc,
                theta_tip=th_tip_calc,
                root_cutout=x0,
            )
            return geometry, pitch

        raise ValueError(f"Unknown PITCH_MODE: '{pitch_mode}'")

    raise ValueError(
        f"Unknown HOVER_TRIM_MODE: '{hover_trim_mode}'. Options: 'collective', 'rpm', or 'none'."
    )


def resolve_pitch(
    geometry,
    pitch_mode=PITCH_MODE,
    ct_hover_target=CT_HOVER_TARGET,
    theta0=THETA0,
    theta_root_deg=THETA_ROOT_DEG,
    theta_tip_deg=THETA_TIP_DEG,
    hover_trim_mode=None,
    thrust_hover_n=THRUST_HOVER_N,
):
    """Compatibility function maintaining full backwards compatibility with existing calls."""
    if pitch_mode == "constant" and ct_hover_target is not None and theta0 is not None and hover_trim_mode is None:
        raise ValueError("Specify exactly one of CT_HOVER_TARGET or THETA0 for PITCH_MODE='constant'")

    trim_mode = hover_trim_mode or ("collective" if (ct_hover_target is not None or thrust_hover_n is not None) else "none")
    _, pitch = trim_hover(
        geometry,
        hover_trim_mode=trim_mode,
        pitch_mode=pitch_mode,
        ct_hover_target=ct_hover_target,
        thrust_hover_n=thrust_hover_n,
        theta0=theta0,
        theta_root_deg=theta_root_deg,
        theta_tip_deg=theta_tip_deg,
    )
    return pitch


def collective_pitch(geometry, ct_hover_target=None, theta0=None):
    """Legacy helper function returning collective pitch angle theta0."""
    pitch = resolve_pitch(geometry, "constant", ct_hover_target, theta0)
    return pitch.theta0


def inflow_gradients(mu, lam, model, fx=FX_COLEMAN, fy=FY_COLEMAN):
    """Calculates non-dimensional spatial inflow gradient coefficients (Kx, Ky).

    Classical first-order linear induced velocity distribution:
      lambda_d(x, psi) = lambda + x * (lambda_1c * cos(psi) + lambda_1s * sin(psi))
    where:
      lambda_1c = Kx * lambda_i  (longitudinal gradient, positive aft)
      lambda_1s = Ky * lambda_i  (lateral gradient, negative retreating side)
      tan(chi/2) = mu / (sqrt(mu^2 + lambda^2) + |lambda|)

    Models:
      - "uniform":          Kx = 0, Ky = 0
      - "coleman_simple":   Kx = tan(chi/2), Ky = 0 (Coleman et al. 1945)
      - "coleman_feingold": NDARC Coleman-Feingold (Kx = fx*(15*pi/32)*tan(chi/2), Ky = -fy*2*mu)
      - "drees":            Kx = (4/3)*(1 - 1.8*mu^2)*tan(chi/2), Ky = -2*mu (Drees 1949)
    """
    if model == "uniform":
        return 0.0, 0.0

    denom = math.sqrt(mu * mu + lam * lam) + abs(lam)
    tan_chi_half = (mu / denom) if denom > 1e-15 else 0.0

    if model == "coleman_simple":
        return tan_chi_half, 0.0
    if model in ("coleman_feingold", "coleman"):
        kx = fx * (15.0 * math.pi / 32.0) * tan_chi_half
        ky = -fy * 2.0 * mu
        return kx, ky
    if model == "drees":
        kx = (4.0 / 3.0) * (1.0 - 1.8 * mu * mu) * tan_chi_half
        ky = -2.0 * mu
        return kx, ky

    raise ValueError(f"Unknown inflow model: '{model}'. Options: {MODELS}")


def coleman_kx(mu, lam, fx=FX_COLEMAN):
    """NDARC Coleman-Feingold longitudinal harmonic factor (maintained for compatibility)."""
    kx, _ = inflow_gradients(mu, lam, "coleman_feingold", fx=fx, fy=0.0)
    return kx


def ct_bet(mu, lam, lambda_1s, moments, geometry, a=None):
    """Calculates analytical BET thrust coefficient CT using weighted radial moments."""
    _, i_mom, t_mom = moments
    lift_slope = geometry.lift_slope(mu) if a is None else float(a)
    return 0.5 * lift_slope * (
        t_mom[2] + 0.5 * mu * mu * t_mom[0]
        - (lam + 0.5 * mu * lambda_1s) * i_mom[1]
    )


def solve_inflow(mu, mu_z, pitch, geometry, model, fx=FX_COLEMAN, fy=FY_COLEMAN):
    """Solves for lambda_i >= 0 via momentum theory and Blade Element Theory equilibrium."""
    b_val = geometry.b_factor()
    moments = radial_moments(geometry, pitch, b=b_val)
    lift_slope = geometry.lift_slope(mu)

    def residual(lambda_i):
        lam = mu_z + lambda_i
        _, ky = inflow_gradients(mu, lam, model, fx=fx, fy=fy)
        lambda_1s = ky * lambda_i
        momentum_ct = 2.0 * (b_val * b_val) * lambda_i * math.sqrt(mu * mu + lam * lam)
        bet_ct = ct_bet(mu, lam, lambda_1s, moments, geometry, a=lift_slope)
        return bet_ct - momentum_ct

    lo = 0.0
    flo = residual(lo)
    if flo < 0.0:
        raise ValueError(
            f"No physical root with positive thrust for mu={mu:g}, mu_z={mu_z:g}, model={model}"
        )
    if abs(flo) < 1e-14:
        return 0.0

    hi = 0.1
    fhi = residual(hi)
    while fhi > 0.0 and hi < 100.0:
        hi *= 2.0
        fhi = residual(hi)
    if fhi > 0.0:
        raise ValueError(f"Could not bracket inflow root for mu={mu:g}, mu_z={mu_z:g}")

    for _ in range(200):
        mid = 0.5 * (lo + hi)
        fmid = residual(mid)
        if abs(fmid) < 1e-13 or hi - lo < 1e-13:
            return mid
        if fmid > 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


@lru_cache(maxsize=None)
def _gauss_nodes(order):
    return np.polynomial.legendre.leggauss(order)


def profile_coefficients(mu, mu_z, geometry, profile_model=PROFILE_MODEL):
    """Calculates profile drag force CH0 and profile torque CQ0.

    Under 'complete', numerically integrates local total velocity (accounting for
    radial velocity u_R, reverse flow, mu_z, and local solidity sigma(r)).
    """
    b_val = geometry.b_factor()
    j = radial_integrals(geometry, b=b_val)
    s0, s1 = geometry.solidity.linear_coeffs
    i_mom = {m: s0 * j[m] + s1 * j[m + 1] for m in range(5)}
    cd0 = geometry.cd0

    if profile_model == "simple_bet":
        ch0 = cd0 * mu * i_mom[1] / 2.0
        cq0 = cd0 / 2.0 * (i_mom[3] + 0.5 * mu * mu * i_mom[1])
        return ch0, cq0

    if profile_model != "complete":
        raise ValueError("profile_model must be 'complete' or 'simple_bet'")

    # 2D Gauss-Legendre quadrature (radial and azimuthal)
    xr, wr = _gauss_nodes(48)
    xp, wp = _gauss_nodes(96)
    r0 = geometry.root_cutout
    r = 0.5 * (b_val - r0) * (xr + 1.0) + r0
    radial_weights = 0.5 * (b_val - r0) * wr
    psi = math.pi * (xp + 1.0)
    azimuth_average_weights = 0.5 * wp

    rr, pp = np.meshgrid(r, psi, indexing="ij")
    weights = radial_weights[:, None] * azimuth_average_weights[None, :]
    u_t = rr + mu * np.sin(pp)
    u_r = mu * np.cos(pp)
    total_speed = np.sqrt(u_t * u_t + u_r * u_r + mu_z * mu_z)

    # Local solidity sigma(r) at radial station
    sigma_r = geometry.solidity.sigma(rr)
    factor = sigma_r * cd0 / 2.0
    ch0 = np.sum(weights * factor * total_speed * (rr * np.sin(pp) + mu))
    cq0 = np.sum(weights * factor * total_speed * u_t * rr)
    return float(ch0), float(cq0)


def coefficients(
    mu,
    mu_z,
    pitch_input,
    geometry,
    model,
    profile_model=PROFILE_MODEL,
    cq_model=CQ_MODEL,
    k_ind=K_IND,
    fx=FX_COLEMAN,
    fy=FY_COLEMAN,
):
    """Calculates non-dimensional aerodynamic coefficients and rotor performance metrics."""
    if isinstance(pitch_input, (int, float)):
        pitch = BladePitch(
            "constant",
            theta0=float(pitch_input),
            theta_root=float(pitch_input),
            theta_tip=float(pitch_input),
            root_cutout=geometry.root_cutout,
        )
    else:
        pitch = pitch_input

    b_val = geometry.b_factor()
    j, i_mom, t_mom = radial_moments(geometry, pitch, b=b_val)
    a = geometry.lift_slope(mu)

    lambda_i = solve_inflow(mu, mu_z, pitch, geometry, model, fx=fx, fy=fy)
    lam = mu_z + lambda_i
    kx, ky = inflow_gradients(mu, lam, model, fx=fx, fy=fy)
    lambda_1c = kx * lambda_i
    lambda_1s = ky * lambda_i

    ct = ct_bet(mu, lam, lambda_1s, (j, i_mom, t_mom), geometry, a=a)
    ch0, cq0 = profile_coefficients(mu, mu_z, geometry, profile_model)

    # Induced longitudinal H-force CHi:
    chi = 0.25 * a * (lam * mu * t_mom[0] + lambda_1s * t_mom[2])

    # Simplified BET induced torque:
    cqi_bet = 0.5 * a * (
        (lam + 0.5 * mu * lambda_1s) * t_mom[2]
        - lam * lam * i_mom[1]
        - 0.5 * (lambda_1c * lambda_1c + lambda_1s * lambda_1s) * i_mom[3]
    )

    if cq_model == "complete":
        if k_ind <= 0.0:
            raise ValueError("K_IND must be positive")
        # Rotor shaft torque energy balance
        cqi = k_ind * lambda_i * ct + mu_z * ct - mu * chi
    elif cq_model == "simple_bet":
        cqi = cqi_bet
    else:
        raise ValueError("cq_model must be 'complete' or 'simple_bet'")

    # Side force CY (lateral projection of normal force):
    cy = -0.25 * a * lambda_1c * (t_mom[2] - 2.0 * lam * i_mom[1])

    # Rolling moment CMx (body axis: +Mx right wing down / advancing side down):
    cmx = -0.5 * a * mu * (t_mom[2] - 0.5 * lam * i_mom[1]) + 0.25 * a * lambda_1s * i_mom[3]

    # Pitching moment CMy (body axis: +My nose up):
    cmy = 0.25 * a * lambda_1c * i_mom[3]

    cq_total = cqi + cq0
    ch_total = chi + ch0

    # 1. Figure of Merit (FoM) in hover: FoM = CT^(3/2) / (sqrt(2) * CQ)
    fom = (ct ** 1.5) / (math.sqrt(2.0) * cq_total) if (ct > 0 and cq_total > 0) else 0.0

    # 2. Total aerodynamic power against the air: CPair = CQ + mu * CH
    cp_air = cq_total + mu * ch_total

    # 3. Effective rotor L/D ratio in forward flight: (L/D)_eff = mu * CT / CPair
    l_d_eff = (mu * ct / cp_air) if (mu > 1e-6 and cp_air > 1e-12) else 0.0

    # 4. Dimensional outputs (Thrust [N] and Powers [kW])
    thrust_n = geometry.thrust_from_ct(ct)
    power_shaft_kw = cq_total * geometry.rho * geometry.disk_area * (geometry.vtip ** 3) / 1000.0
    power_air_kw = cp_air * geometry.rho * geometry.disk_area * (geometry.vtip ** 3) / 1000.0

    return {
        "CT": ct,
        "CQ": cq_total,
        "CQi": cqi,
        "CQ0": cq0,
        "CH": ch_total,
        "CHi": chi,
        "CH0": ch0,
        "CY": cy,
        "CMy": cmy,
        "CMx": cmx,
        "CPair": cp_air,
        "lambda": lam,
        "lambda_i": lambda_i,
        "L_D_eff": l_d_eff,
        "FoM": fom,
        "T_N": thrust_n,
        "P_kW": power_shaft_kw,
        "P_air_kw": power_air_kw,
    }


def run_sweep(
    geometry,
    mus,
    axial_flow,
    axial_values,
    pitch,
    inflow_models=INFLOW_MODELS,
    profile_model=PROFILE_MODEL,
    cq_model=CQ_MODEL,
    k_ind=K_IND,
    fx=FX_COLEMAN,
    fy=FY_COLEMAN,
):
    """Executes advance ratio (mu) sweep across axial conditions for configured inflow models."""
    if isinstance(pitch, (int, float)):
        pitch = BladePitch(
            "constant",
            theta0=float(pitch),
            theta_root=float(pitch),
            theta_tip=float(pitch),
            root_cutout=geometry.root_cutout,
        )

    rows = []
    for axial_input in axial_values:
        for mu in mus:
            mu = float(mu)
            mu_z, condition = axial_condition(mu, axial_input, axial_flow, geometry)
            for model in inflow_models:
                row = {
                    "model": model,
                    "axial_flow": axial_flow,
                    "axial_input": float(axial_input),
                    "condition": condition,
                    "mu": mu,
                    "mu_z": mu_z,
                    "V_m_s": mu * geometry.vtip,
                    "V_km_h": mu * geometry.vtip * 3.6,
                    "theta0_rad": pitch.theta0,
                    "theta0_deg": math.degrees(pitch.theta0),
                    "theta_root_deg": math.degrees(pitch.theta_root),
                    "theta_tip_deg": math.degrees(pitch.theta_tip),
                    "sigma_ref": geometry.solidity.sigma_ref,
                    "sigma_geom": geometry.solidity.sigma_geom,
                    "sigma_thrust": geometry.solidity.sigma_thrust,
                    "profile_model": profile_model,
                    "cq_model": cq_model,
                    "K_ind": k_ind,
                    "RPM": geometry.rpm,
                    "B_tip_loss": geometry.b_factor(),
                    "lift_slope_a": geometry.lift_slope(mu),
                    "Mach_tip": geometry.tip_mach,
                }
                res = coefficients(
                    mu,
                    mu_z,
                    pitch,
                    geometry,
                    model,
                    profile_model=profile_model,
                    cq_model=cq_model,
                    k_ind=k_ind,
                    fx=fx,
                    fy=fy,
                )
                kx_actual, ky_actual = inflow_gradients(mu, res["lambda"], model, fx=fx, fy=fy)
                row["Kx"] = kx_actual
                row["Ky"] = ky_actual
                row.update(res)
                rows.append(row)
    return pd.DataFrame(rows)


def save_results(df, output_dir):
    """Saves results table to CSV with UTF-8 BOM encoding.

    Generates:
      1. outputs/zBET.csv: consolidated dataset containing all inflow models.
      2. outputs/zBET_<model>.csv: individual dataset per inflow model.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. Consolidated CSV file
    consolidated_path = output_dir / "zBET.csv"
    df.to_csv(consolidated_path, index=False, encoding="utf-8-sig")

    # 2. Individual CSV files per inflow model
    for model, group in df.groupby("model", sort=False):
        model_csv = output_dir / f"zBET_{model}.csv"
        group.to_csv(model_csv, index=False, encoding="utf-8-sig")
        if model == "coleman_feingold":
            group.to_csv(output_dir / "zBET_coleman.csv", index=False, encoding="utf-8-sig")

    return consolidated_path


def plot_results(df, output_dir, models_to_plot=None):
    """Generates plots for each aerodynamic coefficient vs advance ratio (mu) for each inflow model."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    display_names = {
        "CPair": "C_Pair (Total Air Power: CQ + μ·CH)",
        "lambda": "λ (Total Mean Inflow)",
        "lambda_i": "λ_i (Induced Mean Inflow)",
        "L_D_eff": "(L/D)_eff (Effective Rotor L/D Ratio)",
    }
    if models_to_plot is None:
        models_to_plot = df["model"].unique()

    for model in models_to_plot:
        model_df = df[df["model"] == model]
        if model_df.empty:
            continue
        model_label = MODEL_LABELS.get(model, model)
        for output in OUTPUTS:
            fig, ax = plt.subplots(figsize=(7.2, 4.6))
            for condition, group in model_df.groupby("condition", sort=False):
                group = group.sort_values("mu")
                ax.plot(group["mu"], group[output], marker="o", ms=3.5, label=condition)
            ax.set_xlabel("Advance ratio, μ")
            ylabel = display_names.get(output, output)
            ax.set_ylabel(ylabel)
            ax.set_title(f"{ylabel} vs. μ — {model_label} Inflow")
            ax.grid(True, alpha=0.35)
            ax.legend()
            fig.tight_layout()
            fig.savefig(output_dir / f"{output.lower()}_vs_mu_{model}.png", dpi=180)
            if model == "coleman_feingold":
                fig.savefig(output_dir / f"{output.lower()}_vs_mu_coleman.png", dpi=180)
            plt.close(fig)


def main():
    solidity = resolve_solidity()
    base_geometry = Geometry(
        RPM,
        R,
        A_LIFT,
        R0_BAR,
        CD0,
        solidity,
        rho=RHO,
        speed_of_sound=SPEED_OF_SOUND,
        tip_loss_mode=TIP_LOSS_MODE,
        tip_loss_b=TIP_LOSS_B,
        use_prandtl_glauert=USE_PRANDTL_GLAUERT,
    )
    geometry, pitch = trim_hover(base_geometry)

    mus = np.round(np.arange(MU_MIN, MU_MAX + 0.5 * MU_STEP, MU_STEP), 10)
    df = run_sweep(
        geometry,
        mus,
        AXIAL_FLOW,
        AXIAL_VALUES,
        pitch,
        inflow_models=INFLOW_MODELS,
        profile_model=PROFILE_MODEL,
        cq_model=CQ_MODEL,
        k_ind=K_IND,
        fx=FX_COLEMAN,
        fy=FY_COLEMAN,
    )
    csv_path = save_results(df, OUTPUT_DIR)
    plot_results(df, OUTPUT_DIR)

    # Hover results summary (mu = 0, alpha = 0):
    hover_row = df[(df["mu"] == 0.0) & (df["axial_input"] == 0.0) & (df["model"] == "uniform")]
    fom_hover = hover_row["FoM"].values[0] if not hover_row.empty else 0.0
    ct_hover = hover_row["CT"].values[0] if not hover_row.empty else 0.0
    cq_hover = hover_row["CQ"].values[0] if not hover_row.empty else 0.0
    thrust_hover = hover_row["T_N"].values[0] if not hover_row.empty else 0.0
    power_hover_kw = hover_row["P_kW"].values[0] if not hover_row.empty else 0.0

    print("=" * 75)
    print("zBET — Execution completed successfully")
    print("=" * 75)
    print(f"Solidity input mode: {solidity.mode}")
    print(f"  sigma_ref    (extrapolated to hub r=0): {solidity.sigma_ref:.4f}")
    print(f"  sigma_geom   (physical area x0 to 1):   {solidity.sigma_geom:.4f}")
    print(f"  sigma_thrust (thrust-weighted):         {solidity.sigma_thrust:.4f}")
    print(f"  Root chord (x0={geometry.root_cutout:.2f}): {solidity.chord_root:.3f} m | Tip chord (x=1.0): {solidity.chord_tip:.3f} m")
    print(f"Hover trim mode: {HOVER_TRIM_MODE}")
    print(f"  Rotor speed: {geometry.rpm:.1f} RPM (Vtip = {geometry.vtip:.1f} m/s | Mach_tip = {geometry.tip_mach:.3f})")
    print(f"Blade pitch mode: {pitch.mode}")
    print(f"  theta_root:   {pitch.theta_root:.6f} rad = {math.degrees(pitch.theta_root):.3f}°")
    print(f"  theta_tip:    {pitch.theta_tip:.6f} rad = {math.degrees(pitch.theta_tip):.3f}°")
    delta_twist_deg = math.degrees(pitch.theta_tip - pitch.theta_root)
    print(f"  total twist:  {delta_twist_deg:.3f}° (total twist invariant)")
    print(f"  theta0/mean:  {pitch.theta0:.6f} rad = {math.degrees(pitch.theta0):.3f}°")
    print(f"Advanced aerodynamic models:")
    print(f"  Tip Loss: {geometry.tip_loss_mode} (B = {geometry.b_factor():.3f})")
    print(f"  Prandtl-Glauert: {geometry.use_prandtl_glauert} (a0 = {geometry.lift_curve_slope:.2f})")
    print(f"Hover performance (mu=0, alpha=0):")
    print(f"  CT = {ct_hover:.5f} | Thrust T = {thrust_hover:.1f} N")
    print(f"  CQ = {cq_hover:.6f} | Shaft Power P = {power_hover_kw:.2f} kW")
    print(f"  Figure of Merit (FoM): {fom_hover:.4f}")
    print(f"Inflow models evaluated: {list(INFLOW_MODELS)}")
    print(f"Sweep rows generated: {len(df)}")
    print(f"Consolidated CSV: {csv_path}")
    print(f"Per-model CSVs: {OUTPUT_DIR}/zBET_<model>.csv")
    print(f"Plots saved to: {OUTPUT_DIR}")
    print("=" * 75)


if __name__ == "__main__":
    main()
