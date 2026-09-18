#!/usr/bin/env python3
"""zBET: Fast Blade Element Theory (BET) rotor solver for rapid conceptual sizing,
parametric trade studies, and forward flight advance ratio (mu) sweeps.

For theoretical derivations, coordinate frame diagrams, sign conventions,
and analytical moment formulations, see 'zBET-documentation.md'.
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
# USER CONFIGURATION
# =============================================================================

# --- 1. Atmosphere & Ambient Fluid -------------------------------------------
RHO = 1.225                 # Air density [kg/m^3]
SPEED_OF_SOUND = 340.3      # Speed of sound [m/s] (used for Mach number & Prandtl-Glauert)

# --- 2. Rotor Geometry & Aerodynamics ----------------------------------------
RPM = 390.0                 # Rotational speed [RPM]
R = 5.5                     # Rotor radius [m]
R0_BAR = 0.15               # Non-dimensional root cutout r0/R [0.0 to 1.0)
A_LIFT = 5.73               # 2D section lift curve slope a0 [1/rad]
CD0 = 0.009                 # 2D section profile drag coefficient [-]

# --- 3. Solidity & Chord -----------------------------------------------------
# Options: "sigma_ref" (reference solidity) | "sigma_geom" (geometric solidity) | "chords" (root & tip chords)
SOLIDITY_MODE = "sigma_ref"
SIGMA_REF = 0.075           # Reference solidity [-] (used if SOLIDITY_MODE = "sigma_ref")
SIGMA_GEOM = 0.06375        # Geometric physical solidity [-] (used if SOLIDITY_MODE = "sigma_geom")
N_BLADES = 4                # Number of blades [-] (used if SOLIDITY_MODE = "chords" or tip loss "sissingh")
CHORD_ROOT = 0.324          # Root chord at station x0 [m] (used if SOLIDITY_MODE = "chords")
CHORD_TIP = 0.324           # Tip chord at station x=1.0 [m] (used if SOLIDITY_MODE = "chords")

# --- 4. Pitch & Hover Trim ---------------------------------------------------
# Options: "constant" (uniform pitch) | "linear_twist" (linear twist from x0 to tip)
PITCH_MODE = "constant"

# Options: "collective" (adjusts theta0 to match hover target) | "rpm" (adjusts RPM to match target thrust) | "none" (prescribed theta0/RPM)
HOVER_TRIM_MODE = "collective"
COLLECTIVE_MODE = HOVER_TRIM_MODE  # Compatibility alias

CT_HOVER_TARGET = 0.0065    # Target hover thrust coefficient [-] (for "collective" mode; set None if using THRUST_HOVER_N)
THRUST_HOVER_N = None       # Target hover thrust [N] (required for "rpm" mode, overrides CT_HOVER_TARGET in "collective")
THETA0 = None               # Prescribed collective pitch [rad] (used if HOVER_TRIM_MODE = "none" & PITCH_MODE = "constant")
THETA_ROOT_DEG = 12.0       # Root pitch [deg] at x0 (used if PITCH_MODE = "linear_twist"; twist preserved during collective trim)
THETA_TIP_DEG = 4.0         # Tip pitch [deg] at x=1.0 (used if PITCH_MODE = "linear_twist")

# --- 5. Aerodynamic Corrections ----------------------------------------------
# Options: "none" (B = 1.0) | "fixed" (uses TIP_LOSS_B) | "sissingh" (B = 1 - sqrt(2*CT)/N_BLADES)
TIP_LOSS_MODE = "none"
TIP_LOSS_B = 0.97           # Tip loss factor B [-] (used if TIP_LOSS_MODE = "fixed")
USE_PRANDTL_GLAUERT = False # Subsonic compressibility correction on lift curve slope (True | False)

# --- 6. Inflow Models & Induced Power ----------------------------------------
# Available models: "uniform", "coleman_simple", "coleman_feingold", "drees"
INFLOW_MODELS = ["uniform", "coleman_simple", "coleman_feingold", "drees"]
FX_COLEMAN = 1.0            # Longitudinal inflow gradient scale factor [-]
FY_COLEMAN = 1.0            # Lateral inflow gradient scale factor [-] (set 0.0 to disable lateral inflow)
K_IND = 1.15                # Induced power factor [-] (CQi = K_IND * lambda_i * CT + ...)

# --- 7. Operating Sweep & Output ---------------------------------------------
MU_MIN = 0.0                # Minimum advance ratio [-]
MU_MAX = 0.40               # Maximum advance ratio [-]
MU_STEP = 0.05              # Advance ratio step size [-]

# Axial flow mode: "alpha" (disk angle of attack [deg]) | "mu_z" (axial velocity [-]) | "w" (vertical climb speed [m/s])
# Convention: alpha > 0 is wind from below (tilted nose-up, increases CT); mu_z > 0 is wind from above (climb)
AXIAL_FLOW = "alpha"
AXIAL_VALUES = [0.0, -4.0, 4.0]  # Values corresponding to AXIAL_FLOW mode

# Aerodynamic model selectors.
# Induced shaft torque:
#   "analytical_bet" -> direct closed-form BET torque
#   "energy_balance" -> energy-balance closure with K_IND
# Profile drag:
#   "analytical_tangential" -> closed-form tangential-only BET profile drag
#   "analytical_vectorial"  -> low-order closed-form vector profile drag
#   "numerical_vectorial"   -> radial/azimuthal vector profile-drag quadrature
# CT from lift, CHi, CY, CMx, and CMy use the analytical weighted-moment formulation.
# Vectorial profile models also compute CT0, the normal profile-drag contribution.
INDUCED_TORQUE_MODEL = "energy_balance"
PROFILE_DRAG_MODEL = "numerical_vectorial"

OUTPUT_DIR = Path("outputs") # Directory for CSV and plot outputs
# =============================================================================


OUTPUTS = [
    "CT", "CT0", "CQ", "CQi", "CQ0", "CH", "CHi", "CH0", "CY",
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
    """Constructs the BladeSolidity object, analytically calculating three solidity metrics:
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


def profile_drag_coefficients(mu, mu_z, geometry, profile_drag_model=PROFILE_DRAG_MODEL):
    """Calculates profile-drag contributions CT0, CH0, and CQ0.

    "analytical_tangential" uses the classical tangential-only small-angle BET
    expressions. It neglects radial and axial profile-drag projections, so CT0=0.

    "analytical_vectorial" retains the leading vector corrections from the local
    velocity components. It is a low-order expansion: CH0 is first order in mu,
    CQ0 is second order in mu and mu_z, and CT0 is first order in mu_z.

    "numerical_vectorial" performs radial/azimuthal Gauss-Legendre quadrature
    of the profile-drag vector using u_T, u_R, imposed axial velocity mu_z,
    and local solidity sigma(r). It evaluates the Johnson profile-force
    projections directly and therefore captures reverse-flow sign changes.

    Profile drag is integrated over the physical blade span x0..1.0. The
    effective tip-loss radius B applies to lift loading, not to blade skin drag.
    """
    # Profile drag exists on the physical blade all the way to x=1, even when
    # an effective tip-loss radius B<1 is used for lift-induced loading.
    j = radial_integrals(geometry, b=1.0)
    s0, s1 = geometry.solidity.linear_coeffs
    i_mom = {m: s0 * j[m] + s1 * j[m + 1] for m in range(5)}
    cd0 = geometry.cd0

    if profile_drag_model == "analytical_tangential":
        ct0 = 0.0
        ch0 = cd0 * mu * i_mom[1] / 2.0
        cq0 = cd0 / 2.0 * (i_mom[3] + 0.5 * mu * mu * i_mom[1])
        return ct0, ch0, cq0

    if profile_drag_model == "analytical_vectorial":
        # Low-order expansion of W = sqrt(u_T^2 + u_R^2 + mu_z^2):
        #   <W (x sin psi + mu)> = (3/2) mu x + O(mu^3)
        #   <W u_T x> = x^3 + (3/4 mu^2 + 1/2 mu_z^2) x + O(4)
        #   <W> (-mu_z) = -mu_z x + O(3)
        ct0 = -0.5 * cd0 * mu_z * i_mom[1]
        ch0 = 0.75 * cd0 * mu * i_mom[1]
        cq0 = 0.5 * cd0 * (
            i_mom[3] + (0.75 * mu * mu + 0.5 * mu_z * mu_z) * i_mom[1]
        )
        return ct0, ch0, cq0

    if profile_drag_model != "numerical_vectorial":
        raise ValueError(
            "profile_drag_model must be 'analytical_tangential', "
            "'analytical_vectorial', or 'numerical_vectorial'"
        )

    # 2D Gauss-Legendre quadrature (radial and azimuthal).
    xr, wr = _gauss_nodes(48)
    xp, wp = _gauss_nodes(96)
    r0 = geometry.root_cutout
    r = 0.5 * (1.0 - r0) * (xr + 1.0) + r0
    radial_weights = 0.5 * (1.0 - r0) * wr
    psi = math.pi * (xp + 1.0)
    azimuth_average_weights = 0.5 * wp

    rr, pp = np.meshgrid(r, psi, indexing="ij")
    weights = radial_weights[:, None] * azimuth_average_weights[None, :]
    u_t = rr + mu * np.sin(pp)
    u_r = mu * np.cos(pp)
    total_speed = np.sqrt(u_t * u_t + u_r * u_r + mu_z * mu_z)

    sigma_r = geometry.solidity.sigma(rr)
    factor = sigma_r * cd0 / 2.0

    # Johnson, Rotorcraft Aeromechanics, Eqs. 6.400-6.402.
    ct0 = np.sum(weights * factor * total_speed * (-mu_z))
    ch0 = np.sum(weights * factor * total_speed * (rr * np.sin(pp) + mu))
    cq0 = np.sum(weights * factor * total_speed * u_t * rr)
    return float(ct0), float(ch0), float(cq0)


def coefficients(
    mu,
    mu_z,
    pitch_input,
    geometry,
    model,
    profile_drag_model=PROFILE_DRAG_MODEL,
    induced_torque_model=INDUCED_TORQUE_MODEL,
    k_ind=K_IND,
    fx=FX_COLEMAN,
    fy=FY_COLEMAN,
):
    """Calculates rotor coefficients using zBET's hybrid BET formulation.

    CT, CHi, CY, CMx, and CMy use analytical weighted radial moments.
    PROFILE_DRAG_MODEL selects tangential analytical, vectorial analytical,
    or vectorial numerical profile drag for CT0/CH0/CQ0. CT is strictly the
    non-viscous thrust coefficient. CT0 is a separate viscous normal-force
    contribution and is never added to CT.
    INDUCED_TORQUE_MODEL selects direct analytical BET or the energy-balance
    closure for CQi.
    """
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
    ct0, ch0, cq0 = profile_drag_coefficients(mu, mu_z, geometry, profile_drag_model)

    # Induced longitudinal H-force CHi:
    chi = 0.25 * a * (lam * mu * t_mom[0] + lambda_1s * (t_mom[2] - 2.0 * lam * i_mom[1]))

    # Analytical BET induced torque:
    cqi_bet = 0.5 * a * (
        (lam + 0.5 * mu * lambda_1s) * t_mom[2]
        - lam * lam * i_mom[1]
        - 0.5 * (lambda_1c * lambda_1c + lambda_1s * lambda_1s) * i_mom[3]
    )

    if induced_torque_model == "energy_balance":
        if k_ind <= 0.0:
            raise ValueError("K_IND must be positive")
        # Rotor shaft torque energy balance.
        # CT is strictly the non-viscous thrust coefficient. CT0 is viscous and
        # does not enter CT or CQi.
        cqi = k_ind * lambda_i * ct + mu_z * ct - mu * chi
    elif induced_torque_model == "analytical_bet":
        cqi = cqi_bet
    else:
        raise ValueError("induced_torque_model must be 'analytical_bet' or 'energy_balance'")

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

    # 2. Aerodynamic power relative to the undisturbed air.
    # Shaft power is CQ. CPair adds the in-plane translational work mu*CH.
    # The climb contribution +mu_z*CT is already embedded in CQ through the
    # energy-balance torque closure. Equivalently:
    # CPair = CQ + mu*CH.
    # Equivalent energy-balance form:
    # CPair = K_IND*lambda_i*CT + mu_z*CT + mu_z*CT0 + CP0_air,
    # where CP0_shaft = CQ0 and
    # CP0_air = CQ0 + mu*CH0 - mu_z*CT0.
    # CT0 is used only in this air-power bookkeeping, never in CT.
    cp_air = cq_total + mu * ch_total

    # 3. Effective rotor L/D ratio in forward flight: (L/D)_eff = mu * CT / CPair
    l_d_eff = (mu * ct / cp_air) if (mu > 1e-6 and cp_air > 1e-12) else 0.0

    # 4. Dimensional outputs (Thrust [N] and Powers [kW])
    thrust_n = geometry.thrust_from_ct(ct)
    power_shaft_kw = cq_total * geometry.rho * geometry.disk_area * (geometry.vtip ** 3) / 1000.0
    power_air_kw = cp_air * geometry.rho * geometry.disk_area * (geometry.vtip ** 3) / 1000.0

    return {
        "CT": ct,
        "CT0": ct0,
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
    profile_drag_model=PROFILE_DRAG_MODEL,
    induced_torque_model=INDUCED_TORQUE_MODEL,
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
                    "profile_drag_model": profile_drag_model,
                    "induced_torque_model": induced_torque_model,
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
                    profile_drag_model=profile_drag_model,
                    induced_torque_model=induced_torque_model,
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
        "CPair": "C_Pair (Air Power: CQ + μ·CH)",
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
        profile_drag_model=PROFILE_DRAG_MODEL,
        induced_torque_model=INDUCED_TORQUE_MODEL,
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
