# zBET — Theory, Implementation, and Model Scope

## 1. Purpose

**zBET** is a fast, semi-empirical Blade Element Theory (BET) rotor solver intended for conceptual design, parametric sweeps, and preliminary performance studies. The implementation combines:

- analytical radial moments for the main blade-element loads;
- global momentum theory for mean induced velocity;
- first-harmonic inflow-gradient models;
- selectable profile-drag formulations;
- selectable direct-BET or energy-balance induced torque; and
- an energy balance for the reported air-power coefficient `CPair`.

The implementation is deliberately lighter than a comprehensive rotorcraft analysis. It does **not** solve blade dynamics, cyclic trim, nonlinear airfoil tables, dynamic stall, or a local momentum equation at each annulus.

### 1.1 Literature basis

The implementation has been checked against the formulation and interpretation used in:

- Wayne Johnson, *Rotorcraft Aeromechanics*, especially Chapter 6 (forward-flight section forces, rotor forces, and power), Chapter 7 (performance), and the profile-power discussion in Section 6.23.
- J. Gordon Leishman, *Principles of Helicopter Aerodynamics*, especially Chapter 3 (blade-element analysis) and Chapter 5 (helicopter performance).

Both references distinguish the **general blade-element force integration** from the **closed-form formulas obtained after simplifying assumptions**. Johnson also shows that force-balance and energy-balance methods are equivalent when they are built from the same assumptions and load model.

### 1.2 Aerodynamic model selectors

zBET exposes two independent selectors.

**Induced shaft torque**

- `analytical_bet`: compute $C_{Qi}$ from the direct BET torque integral;
- `energy_balance`: infer $C_{Qi}$ by reversing the energy balance,

$$
\boxed{
C_{Qi}
=
K_{\mathrm{ind}}\lambda_iC_T
+
\mu_zC_T
-
\mu C_{Hi}
}.
$$

The default is `energy_balance`.

**Profile drag**

- `analytical_tangential`: tangential-only closed-form $C_{H0}$ and $C_{Q0}$;
- `analytical_vectorial`: low-order vectorial closed form;
- `numerical_vectorial`: direct radial/azimuthal vector quadrature.

The default is `numerical_vectorial`.

$K_{\mathrm{ind}}$ is used only by energy-balance quantities: the `energy_balance` torque closure, $C_{Pair}$, and hover figure of merit. It is not inserted into the direct BET torque integral.

---

## 2. Coordinate System and Sign Conventions

### 2.1 Hub axes

The hub-centered Cartesian axes are:

- $+x$: forward;
- $+y$: right / starboard;
- $+z$: downward through the rotor disk.

Positive thrust acts upward, therefore along $-z$.

### 2.2 Rotor azimuth and rotation

Viewed from above, the rotor rotates counter-clockwise. The advancing side is the right side of the disk.

The local nondimensional tangential velocity is

$$
u_T = x + \mu \sin\psi
$$

with

$$
x=\frac{r}{R},
\qquad
\mu=\frac{V_h}{\Omega R}.
$$

At $\psi=90^\circ$,

$$
u_T=x+\mu,
$$

and at $\psi=270^\circ$,

$$
u_T=x-\mu.
$$

The in-plane radial velocity used by the profile quadrature is

$$
u_R=\mu\cos\psi.
$$

### 2.3 Axial flow

The total mean inflow used by the analytical lift model is

$$
\lambda=\mu_z+\lambda_i,
$$

where $\lambda_i\ge 0$ is induced downwash in the $+z$ direction.

The code convention is:

- $\mu_z>0$: imposed relative flow is downward through the disk;
- $\mu_z<0$: imposed relative flow is upward through the disk.

For `AXIAL_FLOW = "alpha"`,

$$
\mu_z=-\mu\tan\alpha.
$$

Therefore a positive rotor angle of attack, as defined in zBET, gives $\mu_z<0$.

### 2.4 Rotor loads

The nondimensional coefficients are

$$
C_T=\frac{T}{\rho A(\Omega R)^2},
\qquad
C_H=\frac{H}{\rho A(\Omega R)^2},
\qquad
C_Y=\frac{Y}{\rho A(\Omega R)^2},
$$

$$
C_Q=\frac{Q}{\rho A(\Omega R)^2R},
\qquad
C_{Mx}=\frac{M_x}{\rho A(\Omega R)^2R},
\qquad
C_{My}=\frac{M_y}{\rho A(\Omega R)^2R},
$$

with

$$
A=\pi R^2.
$$

Sign conventions:

- $C_T>0$: upward rotor thrust;
- $C_H>0$: aft rotor drag;
- $C_Y>0$: force toward starboard;
- $C_Q>0$: positive shaft-torque magnitude required to power the rotor;
- $C_{Mx}>0$: right wing down;
- $C_{My}>0$: nose up.

---

## 3. Blade Geometry, Solidity, and Pitch

The lifting span starts at

$$
x_0=\frac{r_0}{R}.
$$

For a linear chord distribution,

$$
c(x)=c_{\mathrm{root}}
+\left(c_{\mathrm{tip}}-c_{\mathrm{root}}\right)
\frac{x-x_0}{1-x_0}.
$$

The local rotor solidity is

$$
\sigma(x)=\frac{Nc(x)}{\pi R}=s_0+s_1x.
$$

### 3.1 Reference solidity

The reference solidity extrapolates the chord line to the hub:

$$
\sigma_{\mathrm{ref}}
=
\int_0^1 \sigma(x)\,dx
=
s_0+\frac{s_1}{2}.
$$

For a rectangular blade,

$$
\sigma_{\mathrm{ref}}=\frac{Nc}{\pi R}.
$$

### 3.2 Physical geometric solidity

The actual blade area begins at the root cutout:

$$
\sigma_{\mathrm{geom}}
=
\int_{x_0}^1\sigma(x)\,dx.
$$

For a rectangular blade,

$$
\sigma_{\mathrm{geom}}
=
(1-x_0)\sigma_{\mathrm{ref}}.
$$

### 3.3 Thrust-weighted solidity

zBET also reports

$$
\sigma_{\mathrm{thrust}}
=
3\int_{x_0}^1\sigma(x)x^2\,dx.
$$

For a rectangular blade,

$$
\sigma_{\mathrm{thrust}}
=
(1-x_0^3)\sigma_{\mathrm{ref}}.
$$

### 3.4 Why the BET equations use local solidity

The section force is proportional to local chord. Therefore the BET integrands are naturally scaled by

$$
\sigma(x)=\frac{Nc(x)}{\pi R}.
$$

The root cutout is already represented by the integration limits. Replacing local solidity by the physical-area solidity inside an integral that already begins at $x_0$ would apply the root-cutout penalty twice.

### 3.5 Pitch

For constant pitch,

$$
\theta(x)=\theta_0.
$$

For linear twist,

$$
\theta(x)=t_0+t_1x.
$$

Collective trim shifts the complete pitch distribution by a constant and preserves the specified twist.

---

## 4. Inflow Models

The local induced-flow model is

$$
\lambda_d(x,\psi)
=
\lambda
+
x\left(
\lambda_{1c}\cos\psi
+
\lambda_{1s}\sin\psi
\right),
$$

with

$$
\lambda_{1c}=K_x\lambda_i,
\qquad
\lambda_{1s}=K_y\lambda_i.
$$

The wake-skew helper used by the code is

$$
\tan\frac{\chi}{2}
=
\frac{\mu}
{\sqrt{\mu^2+\lambda^2}+|\lambda|}.
$$

### 4.1 Uniform

$$
K_x=0,
\qquad
K_y=0.
$$

### 4.2 Coleman simple

$$
K_x=\tan\frac{\chi}{2},
\qquad
K_y=0.
$$

### 4.3 Coleman-Feingold / NDARC form

$$
K_x
=
f_x\frac{15\pi}{32}
\tan\frac{\chi}{2},
\qquad
K_y=-2f_y\mu.
$$

### 4.4 Drees

$$
K_x
=
\frac{4}{3}
\left(1-1.8\mu^2\right)
\tan\frac{\chi}{2},
\qquad
K_y=-2\mu.
$$

### 4.5 Mean momentum closure

For a given blade pitch, zBET solves $\lambda_i$ from the intersection of analytical BET thrust and global momentum theory:

$$
C_T
=
2B^2\lambda_i
\sqrt{\mu^2+\lambda^2}.
$$

This is a **global** closure. zBET is therefore BET plus global momentum theory, not a multi-annulus BEMT solver.

---

## 5. Hover Trim

Hover trim is performed at

$$
\mu=0,
\qquad
\mu_z=0.
$$

### 5.1 Collective trim

With fixed RPM, the pitch distribution is shifted until the requested hover $C_T$ or dimensional thrust is reached.

### 5.2 RPM trim

With fixed pitch, RPM is adjusted using

$$
T
=
C_T\rho A(\Omega R)^2.
$$

Therefore

$$
\Omega
=
\frac{1}{R}
\sqrt{\frac{T}{\rho A C_T}}.
$$

### 5.3 No trim

With `HOVER_TRIM_MODE = "none"`, the specified RPM and pitch are used directly.

---

## 6. Aerodynamic Formulation

### 6.1 General blade-element reference model

The general section-level BET picture in Johnson and Leishman starts from the local relative velocity, section angle of attack, lift, and drag. In a two-dimensional blade section plane,

$$
U_{2D}
=
\sqrt{u_T^2+u_P^2},
$$

$$
\phi
=
\tan^{-1}\left(\frac{u_P}{u_T}\right),
\qquad
\alpha_s=\theta-\phi.
$$

The local lift and drag can then be resolved into normal and in-plane section forces:

$$
dF_z
=
dL\cos\phi-dD\sin\phi,
$$

$$
dF_x
=
dL\sin\phi+dD\cos\phi.
$$

A **fully integrated force-balance solver** would evaluate a common local aerodynamic state and consistently integrate the resulting forces for thrust, in-plane forces, torque, and hub moments over radius and azimuth.

zBET does not currently expose a fully integrated force-balance model.

### 6.2 Analytical weighted moments used by zBET

The principal zBET load model applies the standard small-inflow-angle analytical reduction and evaluates the radial dependence exactly through weighted moments.

Define

$$
J_n
=
\int_{x_0}^{B}x^n\,dx
=
\frac{B^{n+1}-x_0^{n+1}}{n+1},
$$

$$
I_m
=
\int_{x_0}^{B}\sigma(x)x^m\,dx,
$$

and

$$
T_m
=
\int_{x_0}^{B}\sigma(x)\theta(x)x^m\,dx.
$$

Because $\sigma(x)$ and $\theta(x)$ are linear in $x$, these moments are evaluated analytically.

The implemented thrust coefficient is

$$
C_T
=
\frac{a}{2}
\left[
T_2
+
\frac{\mu^2}{2}T_0
-
\left(
\lambda
+
\frac{\mu\lambda_{1s}}{2}
\right)I_1
\right].
$$

The induced longitudinal force is

$$
C_{Hi}
=
\frac{a}{4}
\left[
\lambda\mu T_0
+
\lambda_{1s}
\left(
T_2-2\lambda I_1
\right)
\right].
$$

The side force is

$$
C_Y
=
-\frac{a\lambda_{1c}}{4}
\left(
T_2-2\lambda I_1
\right).
$$

The rolling moment is

$$
C_{Mx}
=
-\frac{a\mu}{2}
\left(
T_2-\frac{\lambda I_1}{2}
\right)
+
\frac{a\lambda_{1s}}{4}I_3.
$$

The pitching moment is

$$
C_{My}
=
\frac{a\lambda_{1c}}{4}I_3.
$$

These expressions are independent of the profile-drag torque integration described below.

### 6.3 Profile-drag selector

All profile models return $C_{H0}$ and $C_{Q0}$.

#### Tangential analytical model

Neglecting the radial component of profile drag,

$$
\boxed{
C_{H0}^{\mathrm{tang}}
=
\frac{C_{d0}\mu}{2}I_1
}
$$

and

$$
\boxed{
C_{Q0}^{\mathrm{tang}}
=
\frac{C_{d0}}2
\left(
I_3+\frac{\mu^2}{2}I_1
\right)
}.
$$

For a rectangular blade without root cutout,

$$
C_{H0}^{\mathrm{tang}}=\frac{\sigma C_{d0}}4\mu,
\qquad
C_{Q0}^{\mathrm{tang}}=\frac{\sigma C_{d0}}8(1+\mu^2).
$$

#### Low-order vectorial analytical model

Retaining the leading vector corrections gives

$$
\boxed{
C_{H0}^{\mathrm{vec}}
=
\frac{3C_{d0}\mu}{4}I_1
}
$$

and

$$
\boxed{
C_{Q0}^{\mathrm{vec}}
=
\frac{C_{d0}}2
\left[
I_3+
\left(
\frac34\mu^2+\frac12\mu_z^2
\right)I_1
\right]
}.
$$

For a rectangular blade in edgewise flight without root cutout,

$$
C_{H0}^{\mathrm{vec}}\simeq\frac{3\sigma C_{d0}}8\mu,
\qquad
C_{Q0}^{\mathrm{vec}}\simeq\frac{\sigma C_{d0}}8(1+1.5\mu^2).
$$

#### Numerical vectorial model

The numerical model evaluates

$$
\boxed{
C_{H0}
=
\frac{C_{d0}}2
\int_{x_0}^{1}
\sigma(x)
\left\langle
W(x\sin\psi+\mu)
\right\rangle_\psi dx
}
$$

and

$$
\boxed{
C_{Q0}
=
\frac{C_{d0}}2
\int_{x_0}^{1}
\sigma(x)
\left\langle
Wu_Tx
\right\rangle_\psi dx
}
$$

with

$$
u_T=x+\mu\sin\psi,
\qquad
u_R=\mu\cos\psi,
\qquad
W=\sqrt{u_T^2+u_R^2+\mu_z^2}.
$$

The physical blade span $x_0\le x\le1$ is used for profile drag.
### 6.4 Induced shaft torque selector

For `analytical_bet`, $C_{Qi}$ is obtained directly from the moment of the lift-induced in-plane force about the shaft:

$$
\boxed{
C_{Qi}^{\mathrm{BET}}
=
\frac12
\int_{x_0}^{B}
\sigma(x)a
\left[
\left(\lambda+\frac{\mu\lambda_{1s}}2\right)
\theta(x)x^2
-
\lambda^2x
-
\frac{\lambda_{1c}^2+\lambda_{1s}^2}{2}x^3
\right]dx
}.
$$

For `energy_balance`, the same shaft-torque component is inferred from the power balance:

$$
\boxed{
C_{Qi}^{\mathrm{EB}}
=
K_{\mathrm{ind}}\lambda_iC_T
+
\mu_zC_T
-
\mu C_{Hi}
}.
$$

Thus $K_{\mathrm{ind}}$ appears in the energy-balance torque mode, but not in the direct BET torque mode.
### 6.5 Total shaft torque

$$
\boxed{C_Q=C_{Qi}+C_{Q0}}.
$$

$C_Q$ always remains $C_{Qi}+C_{Q0}$. Which expression supplies $C_{Qi}$ and $C_{Q0}$ is controlled by the two selectors above.

### 6.6 Scope of the mixed analytical/numerical formulation

Lift-induced loads remain analytical weighted-moment expressions, while $C_{H0}$ and $C_{Q0}$ use direct vectorial quadrature.

---

## 7. Performance Metrics

### 7.1 Shaft power

$$
P_{\mathrm{shaft}}=Q\Omega.
$$

In coefficient form,

$$
C_{P,\mathrm{shaft}}=C_Q.
$$

No $\mu C_H$ term belongs to shaft power: $C_Q$ is the mechanical torque coefficient about the rotor axis.

### 7.2 Air power (`CPair`)

Power is introduced here, after the torque calculation.

For uniform induced velocity,

$$
C_{Pi}=K_{\mathrm{ind}}\lambda_iC_T.
$$

The climb contribution is

$$
C_{Pc}=\mu_zC_T.
$$

The in-plane translational contribution associated with the lift-generated longitudinal force is

$$
C_{P,\mathrm{trans}}=\mu C_{Hi}.
$$

Keeping profile power unexpanded first,

$$
\boxed{
C_{Pair}
=
K_{\mathrm{ind}}\lambda_iC_T
+
\mu_zC_T
+
\mu C_{Hi}
+
C_{P0,\mathrm{air}}
}.
$$

Now expand the profile term:

$$
\boxed{
C_{P0,\mathrm{air}}=C_{Q0}+\mu C_{H0}
}.
$$

Therefore

$$
\boxed{
C_{Pair}
=
K_{\mathrm{ind}}\lambda_iC_T
+
\mu_zC_T
+
\mu C_{Hi}
+
C_{Q0}
+
\mu C_{H0}
}.
$$

Since $C_H=C_{Hi}+C_{H0}$,

$$
C_{Pair}
=
K_{\mathrm{ind}}\lambda_iC_T
+
\mu_zC_T
+
\mu C_H
+
C_{Q0}.
$$

For a rectangular blade in edgewise flight,

$$
C_{Q0}\simeq\frac{\sigma C_{d0}}8(1+1.5\mu^2),
\qquad
C_{H0}\simeq\frac{3\sigma C_{d0}}8\mu.
$$

Hence

$$
\begin{aligned}
C_{P0,\mathrm{air}}
&=C_{Q0}+\mu C_{H0}\\
&\simeq
\frac{\sigma C_{d0}}8
\left(1+1.5\mu^2+3\mu^2\right),
\end{aligned}
$$

and therefore

$$
\boxed{
C_{P0,\mathrm{air}}
\simeq
\frac{\sigma C_{d0}}8(1+4.5\mu^2)
}.
$$

The $4.5\mu^2$ factor belongs to profile power relative to the air, not to $C_{Q0}$.
### 7.3 Effective rotor lift-to-drag ratio

For $\mu>0$,

$$
\left(\frac{L}{D}\right)_{\mathrm{eff}}
=
\frac{\mu C_T}{C_{Pair}}.
$$

### 7.4 Hover figure of merit

The ideal hover power coefficient is

$$
C_{P,\mathrm{ideal}}
=
\frac{C_T^{3/2}}{\sqrt2}.
$$

The performance model uses

$$
C_{P,\mathrm{hover}}
=
K_{\mathrm{ind}}
\frac{C_T^{3/2}}{\sqrt2}
+
C_{Q0},
$$

so

$$
\boxed{
FoM
=
\frac{C_T^{3/2}/\sqrt2}
{K_{\mathrm{ind}}C_T^{3/2}/\sqrt2+C_{Q0}}
}.
$$

This definition keeps $K_{\mathrm{ind}}$ in the hover energy balance even when the direct BET torque selector is used.

---

## 8. Engineering Corrections

### 8.1 Tip loss

zBET can use an effective aerodynamic radius $B$ and integrate the **lift-induced analytical moments** from $x_0$ to $B$. Profile drag is integrated over the physical blade span $x_0$ to $1$, because the blade material still produces drag outside the effective lift radius.

The optional Sissingh-style relation is

$$
B
=
1-\frac{\sqrt{2C_T}}{N}.
$$

The momentum closure also uses the effective area factor $B^2$.

> This is an **effective-radius engineering approximation**. It should not be interpreted as a literal implementation of the full Prandtl finite-blade circulation/inflow correction.

### 8.2 Prandtl-Glauert lift-slope correction

If enabled, zBET modifies the linear section lift-curve slope using a representative subsonic Mach number:

$$
a(M)
=
\frac{a_0}
{\sqrt{1-M_{\mathrm{eff}}^2}}.
$$

The representative Mach construction and the numerical cap used in the code are engineering approximations, not a local compressible-airfoil solution.

---

## 9. Assumptions and Limitations

The current implementation is intentionally compact. Important limitations are:

- small-angle analytical treatment for the principal lift-induced rotor loads;
- linear section lift curve;
- constant `CD0` rather than an airfoil polar;
- no nonlinear stall or dynamic stall;
- no cyclic-pitch trim;
- no blade flapping solution in the aerodynamic load calculation;
- no elastic blade motion;
- no local annular momentum iteration;
- no full reverse-flow airfoil model;
- first-harmonic prescribed inflow gradients rather than a free wake;
- the vectorial profile integration uses constant `CD0`; it does not include a local airfoil polar;
- the vectorial profile path uses imposed $\mu_z$ in its local speed magnitude rather than local induced normal velocity;
- lift-induced loads remain analytical even when profile drag is numerical;

These limitations are compatible with the intended use of zBET as a rapid conceptual-analysis tool. They should be considered before applying the code to high advance ratio, severe descent, stalled conditions, or detailed loads work.

---

## 10. Configuration Guide

### Geometry and atmosphere

- `RHO`: air density
- `SPEED_OF_SOUND`: speed of sound
- `RPM`: rotor speed
- `R`: rotor radius
- `R0_BAR`: root cutout
- `A_LIFT`: linear lift-curve slope
- `CD0`: constant profile drag coefficient

### Solidity and planform

- `SOLIDITY_MODE = "sigma_ref"`
- `SOLIDITY_MODE = "sigma_geom"`
- `SOLIDITY_MODE = "chords"`

### Pitch and trim

- `PITCH_MODE = "constant"`
- `PITCH_MODE = "linear_twist"`
- `HOVER_TRIM_MODE = "collective"`
- `HOVER_TRIM_MODE = "rpm"`
- `HOVER_TRIM_MODE = "none"`

### Inflow

- `uniform`
- `coleman_simple`
- `coleman_feingold`
- `drees`

### Torque and profile paths

- `INDUCED_TORQUE_MODEL = "analytical_bet"`
- `INDUCED_TORQUE_MODEL = "energy_balance"` (default)
- `PROFILE_DRAG_MODEL = "analytical_tangential"`
- `PROFILE_DRAG_MODEL = "analytical_vectorial"`
- `PROFILE_DRAG_MODEL = "numerical_vectorial"` (default)
- `K_IND = 1.15` by default

`K_IND` affects the energy-balance torque mode, `CPair`, and hover FoM. It does not modify the direct BET torque integral.

---

## 11. References

1. Wayne Johnson, *Rotorcraft Aeromechanics*, Cambridge University Press, 2013. See especially Chapter 6, including the section-force relations around Eqs. 6.39–6.46, rotor-force decomposition around Eqs. 6.71–6.76, and the power/energy-balance development around Eqs. 6.105–6.115; also Chapter 7 and Section 6.23.
2. J. Gordon Leishman, *Principles of Helicopter Aerodynamics*, Cambridge University Press. See especially Chapter 3, “Blade Element Analysis,” and Chapter 5, “Basic Helicopter Performance.”
3. Wayne Johnson, *NDARC — NASA Design and Analysis of Rotorcraft: Theory*, NASA/TP-2009-215402.
4. R. P. Coleman, A. M. Feingold, and C. W. Stempin, *Evaluation of the Induced-Velocity Field of an Idealized Helicopter Rotor*, NACA ARR L5E10, 1945.
5. J. M. Drees, “A Theory of Airflow Through Rotors and Its Application to Some Helicopter Problems,” 1949.
