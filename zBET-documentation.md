# zBET — Theory, Implementation, and Model Scope

## 1. Purpose

**zBET** is a fast, semi-empirical Blade Element Theory (BET) rotor solver intended for conceptual design, parametric sweeps, and preliminary performance studies. The implementation combines:

- analytical radial moments for the main blade-element loads;
- global momentum theory for mean induced velocity;
- first-harmonic inflow-gradient models;
- optional numerical quadrature for profile drag; and
- an optional energy-balance closure for induced shaft torque.

The implementation is deliberately lighter than a comprehensive rotorcraft analysis. It does **not** solve blade dynamics, cyclic trim, nonlinear airfoil tables, dynamic stall, or a local momentum equation at each annulus.

### 1.1 Literature basis

The implementation has been checked against the formulation and interpretation used in:

- Wayne Johnson, *Rotorcraft Aeromechanics*, especially Chapter 6 (forward-flight section forces, rotor forces, and power), Chapter 7 (performance), and the profile-power discussion in Section 6.23.
- J. Gordon Leishman, *Principles of Helicopter Aerodynamics*, especially Chapter 3 (blade-element analysis) and Chapter 5 (helicopter performance).

Both references distinguish the **general blade-element force integration** from the **closed-form formulas obtained after simplifying assumptions**. Johnson also shows that force-balance and energy-balance methods are equivalent when they are built from the same assumptions and load model.

### 1.2 Aerodynamic model selectors

zBET uses two independent selectors whose names describe the physics directly.

| Selector | `"analytical_bet"` | Higher-fidelity alternative |
| --- | --- | --- |
| `PROFILE_DRAG_MODEL` | `"analytical_tangential"`: classical tangential-only closed form | `"analytical_vectorial"`: low-order vector closed form; `"numerical_vectorial"`: radial × azimuthal vector quadrature |
| `INDUCED_TORQUE_MODEL` | direct analytical BET induced torque | `"energy_balance"`: energy-balance shaft-torque closure with `K_IND` |

The defaults are

```python
PROFILE_DRAG_MODEL = "numerical_vectorial"
INDUCED_TORQUE_MODEL = "energy_balance"
```

These selectors are intentionally independent. zBET defines $C_T$ strictly as the **non-viscous thrust coefficient**. The profile selector separately supplies the viscous normal-force component $C_{T0}$ together with $C_{H0}$ and $C_{Q0}$. $C_{T0}$ is never added to $C_T$ and never enters the inflow or induced-torque calculation; it is used only in the $C_{Pair}$ air-power bookkeeping.

There is no generic “simple” or “complete” aerodynamic mode.

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

These expressions are independent of `PROFILE_DRAG_MODEL` and `INDUCED_TORQUE_MODEL`.

### 6.3 Profile-drag models

Johnson's profile-force formulation resolves the local drag vector into normal, longitudinal, and shaft-torque contributions. With

$$
u_T=x+\mu\sin\psi,
\qquad
u_R=\mu\cos\psi,
$$

and the imposed axial component $\mu_z$,

$$
W=\sqrt{u_T^2+u_R^2+\mu_z^2}.
$$

For constant section drag coefficient, the vectorial profile contributions are

$$
C_{T0}
=
\frac{C_{d0}}{2}
\int \sigma(x)
\left\langle W(-\mu_z)\right\rangle_\psi dx,
$$

$$
C_{H0}
=
\frac{C_{d0}}{2}
\int \sigma(x)
\left\langle W(x\sin\psi+\mu)\right\rangle_\psi dx,
$$

$$
C_{Q0}
=
\frac{C_{d0}}{2}
\int \sigma(x)
\left\langle W u_T x\right\rangle_\psi dx.
$$

The profile drag is integrated over the **physical blade span** $x_0\le x\le1$. An effective lift tip-loss radius $B<1$ does not truncate skin-friction/profile drag on the material blade.

#### Analytical tangential model

With `PROFILE_DRAG_MODEL = "analytical_tangential"`, radial and axial profile-drag projections are neglected. The implemented formulas are

$$
C_{T0}=0,
$$

$$
C_{H0}
=
\frac{C_{d0}\mu I_1}{2},
$$

$$
C_{Q0}
=
\frac{C_{d0}}{2}
\left(
I_3+\frac{\mu^2I_1}{2}
\right).
$$

For a rectangular blade with $x_0=0$,

$$
C_{H0}
=
\frac{\sigma C_{d0}}{4}\mu,
\qquad
C_{Q0}
=
\frac{\sigma C_{d0}}{8}(1+\mu^2).
$$

The corresponding edgewise profile power relative to the air is

$$
C_{P0,\mathrm{air}}
=
C_{Q0}+\mu C_{H0}
=
\frac{\sigma C_{d0}}{8}(1+3\mu^2).
$$

#### Analytical vectorial model

With `PROFILE_DRAG_MODEL = "analytical_vectorial"`, zBET uses the low-order expansion

$$
W
\simeq
x+\mu\sin\psi
+\frac{\mu^2\cos^2\psi+\mu_z^2}{2x}.
$$

After azimuthal averaging, the retained terms give

$$
C_{T0}
=
-\frac{C_{d0}\mu_z I_1}{2},
$$

$$
C_{H0}
=
\frac{3C_{d0}\mu I_1}{4},
$$

$$
C_{Q0}
=
\frac{C_{d0}}{2}
\left[
I_3+
\left(
\frac{3}{4}\mu^2+\frac{1}{2}\mu_z^2
\right)I_1
\right].
$$

For a rectangular blade with $x_0=0$ and $\mu_z=0$,

$$
C_{H0}
=
\frac{3\sigma C_{d0}}{8}\mu,
$$

$$
C_{Q0}
=
\frac{\sigma C_{d0}}{8}
(1+1.5\mu^2).
$$

Therefore

$$
C_{P0,\mathrm{air}}
=
C_{Q0}+\mu C_{H0}
=
\frac{\sigma C_{d0}}{8}
(1+4.5\mu^2).
$$

This distinction is essential: the familiar $1+4.5\mu^2$ (and Bennett's $1+4.65\mu^2$ approximation) is a **profile-power** factor, not a shaft-torque factor. Assigning it directly to $C_{Q0}$ double-counts the translational work $\mu C_{H0}$.

#### Numerical vectorial model

With `PROFILE_DRAG_MODEL = "numerical_vectorial"`, zBET evaluates the three vector projections above directly with radial × azimuthal Gauss-Legendre quadrature. It therefore captures radial-flow effects, reverse-flow sign changes in $u_T$, the axial profile-force contribution $C_{T0}$, and the exact local speed magnitude for the stated constant-$C_{d0}$ model.

The numerical model still uses the imposed axial component $\mu_z$ in the profile-drag kinematics rather than the local induced normal velocity. It is therefore a vectorial **profile-drag** model, not a full nonlinear section-force solver.

The reported thrust coefficient is $C_T$, defined exclusively by the non-viscous loading. The profile model additionally exposes $C_{T0}$ as a separate viscous normal-force component. It is **not** added to $C_T$. The mean-inflow solve and induced-torque model use $C_T$ directly.


### 6.4 Induced shaft-torque models

#### Direct analytical BET torque

With `INDUCED_TORQUE_MODEL = "analytical_bet"`,

$$
C_{Qi}
=
\frac{a}{2}
\left[
\left(
\lambda+\frac{\mu\lambda_{1s}}{2}
\right)T_2
-
\lambda^2I_1
-
\frac{
\lambda_{1c}^2+\lambda_{1s}^2
}{2}I_3
\right].
$$

This is the direct analytical small-angle BET torque expression used by the code.

#### Energy-balance torque

With `INDUCED_TORQUE_MODEL = "energy_balance"`,

$$
C_{Qi}
=
K_{\mathrm{ind}}\lambda_i C_T
+
\mu_z C_T
-
\mu C_{Hi}.
$$

This is an **energy-balance closure for shaft torque**. It is not a numerical integral of local lift and drag torque.

The factor

$$
K_{\mathrm{ind}}=\texttt{K_IND}
$$

allows an empirical induced-power correction. Because the default is

$$
K_{\mathrm{ind}}=1.15,
$$

the energy-balance result is intentionally different from ideal analytical BET.

The total shaft-torque coefficient is

$$
C_Q=C_{Qi}+C_{Q0}.
$$

### 6.5 Energy consistency of the `energy_balance` torque model

Shaft power is

$
C_{P,\mathrm{shaft}}=C_Q.
$

zBET defines the air-power bookkeeping quantity `CPair` by adding the in-plane translational work to shaft power:

$
\boxed{
C_{Pair}=C_Q+\mu C_H
}.
$

There is no separate $-\mu_zC_T$ correction in this torque route. The climb contribution is already embedded in $C_Q$ through the induced-torque balance.

Using

$
C_H=C_{Hi}+C_{H0},
$

and `INDUCED_TORQUE_MODEL = "energy_balance"`,

$
C_{Qi}
=
K_{\mathrm{ind}}\lambda_iC_T
+\mu_zC_T
-\mu C_{Hi}.
$

Therefore the torque route gives

$
\begin{aligned}
C_{Pair}
&=C_Q+\mu C_H \\
&=
K_{\mathrm{ind}}\lambda_iC_T
+\mu_zC_T
+C_{Q0}
+\mu C_{H0}.
\end{aligned}
$

The vector profile-power identity is

$
C_{P0,\mathrm{air}}
=
C_{Q0}+\mu C_{H0}-\mu_zC_{T0}.
$

The equivalent energy-balance route is

$
\boxed{
C_{Pair}
=
K_{\mathrm{ind}}\lambda_iC_T
+\mu_zC_T+\mu_zC_{T0}
+C_{P0,\mathrm{air}}
}.
$

This is the equivalent **energy-balance route**. Here $C_T$ remains purely non-viscous. The separate $+\mu_zC_{T0}$ in the axial-work term cancels the $-\mu_zC_{T0}$ contained in $C_{P0,\mathrm{air}}$, so $C_{T0}$ affects only the $C_{Pair}$ bookkeeping and never redefines $C_T$. The two routes are algebraically identical:

$
\boxed{
C_Q+\mu C_H
=
K_{\mathrm{ind}}\lambda_iC_T
+\mu_zC_T+\mu_zC_{T0}
+C_{P0,\mathrm{air}}
}.
$


### 6.6 Relationship between the selector choices

The two selectors are independent. $C_T$, $C_{Hi}$, $C_Y$, $C_{Mx}$, and $C_{My}$ retain the same analytical weighted-moment formulation, while `PROFILE_DRAG_MODEL` changes the profile-drag calculation and `INDUCED_TORQUE_MODEL` changes the induced shaft-torque calculation.

The tangential, analytical-vectorial, and numerical-vectorial profile models are not required to agree numerically. The analytical-vectorial model is the low-order expansion of the same vector kinematics used by the numerical-vectorial model, so they approach one another as $\mu$ and $|\mu_z|$ become small. The tangential model intentionally omits radial and axial profile-force components.

### 6.7 Fully integrated force balance

A fully integrated force-balance formulation would use the same local velocity, angle of attack, lift/drag model, reverse-flow treatment, and induced-flow state for every force and moment channel, then integrate all of them over radius and azimuth.

Doing this robustly requires additional modeling choices that zBET currently avoids, especially:

- reverse-flow lift and drag;
- nonlinear $C_l(\alpha,M,Re)$ and $C_d(\alpha,M,Re)$;
- stall;
- blade flapping and cyclic pitch;
- the distinction between spanwise flow effects on lift and drag;
- consistent local induced velocity in the section kinematics.

zBET therefore keeps its current mixed analytical/numerical architecture explicit instead of presenting it as a full force-balance solver.

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

zBET reports air power by either of two exactly equivalent routes.

**Torque route**

$
\boxed{
C_{Pair}=C_Q+\mu C_H
}.
$

Here $C_Q$ is the shaft-power coefficient. The term $\mu C_H$ is the in-plane translational work and is not contained in $C_Q$.

**Energy-balance route**

$
\boxed{
C_{Pair}
=
K_{\mathrm{ind}}\lambda_iC_T
+\mu_zC_T+\mu_zC_{T0}
+C_{P0,\mathrm{air}}
},
$

with

$
C_{P0,\mathrm{air}}=C_{Q0}+\mu C_{H0}-\mu_zC_{T0}.
$

For $\mu_z>0$ (climb in the zBET convention), the axial-work term is $+\mu_zC_T+\mu_zC_{T0}$. The $C_{T0}$ part exists only inside the $C_{Pair}$ energy bookkeeping and cancels the $-\mu_zC_{T0}$ contained in $C_{P0,\mathrm{air}}$.


### 7.3 Effective rotor lift-to-drag ratio

For $\mu>0$,

$$
\left(\frac{L}{D}\right)_{\mathrm{eff}}
=
\frac{\mu C_T}
{C_Q+\mu C_H}.
$$

### 7.4 Hover figure of merit

$$
FoM
=
\frac{C_T^{3/2}}
{\sqrt{2}\,C_Q}.
$$

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
- the vectorial profile models use constant `CD0`; they do not include a local airfoil polar;
- the numerical-vectorial path uses imposed $\mu_z$ in profile kinematics rather than local induced normal velocity;
- lift-induced loads remain analytical even when profile drag is numerical;
- `INDUCED_TORQUE_MODEL = "energy_balance"` is an energy-balance closure, not direct torque quadrature.

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

- `INDUCED_TORQUE_MODEL = "analytical_bet"`: direct analytical induced torque
- `INDUCED_TORQUE_MODEL = "energy_balance"`: energy-balance induced torque
- `PROFILE_DRAG_MODEL = "analytical_tangential"`: classical tangential-only profile drag
- `PROFILE_DRAG_MODEL = "analytical_vectorial"`: low-order vectorial profile drag
- `PROFILE_DRAG_MODEL = "numerical_vectorial"`: numerical vector profile-drag quadrature with $C_{T0}$

---

## 11. References

1. Wayne Johnson, *Rotorcraft Aeromechanics*, Cambridge University Press, 2013. See especially Chapter 6, including the section-force relations around Eqs. 6.39–6.46, rotor-force decomposition around Eqs. 6.71–6.76, and the power/energy-balance development around Eqs. 6.105–6.115; also Chapter 7 and Section 6.23.
2. J. Gordon Leishman, *Principles of Helicopter Aerodynamics*, Cambridge University Press. See especially Chapter 3, “Blade Element Analysis,” and Chapter 5, “Basic Helicopter Performance.”
3. Wayne Johnson, *NDARC — NASA Design and Analysis of Rotorcraft: Theory*, NASA/TP-2009-215402.
4. R. P. Coleman, A. M. Feingold, and C. W. Stempin, *Evaluation of the Induced-Velocity Field of an Idealized Helicopter Rotor*, NACA ARR L5E10, 1945.
5. J. M. Drees, “A Theory of Airflow Through Rotors and Its Application to Some Helicopter Problems,” 1949.
