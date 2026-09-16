# zBET — Advanced Rotor Blade Element Momentum Theory Solver

## 1. Scope and Overview

`zBET.py` is an aerodynamic analysis and simulation tool for helicopter and rotary-wing blades based on Blade Element Theory (BET) coupled with Momentum Theory (BEMT). It calculates rotor forces, moments, power requirements, and overall efficiency metrics across advance ratio sweeps ($\mu$) and axial flight conditions ($\alpha$, $\mu_z$, or vertical climb/descent velocity $w$).

The primary aerodynamic coefficients and performance metrics evaluated are:

$$
C_T,\ C_Q,\ C_{Qi},\ C_{Q0},\ C_H,\ C_{Hi},\ C_{H0},\ C_Y,\ C_{My},\ C_{Mx},\ \lambda,\ \lambda_i,\ C_{P,\mathrm{air}},\ (L/D)_{\mathrm{eff}},\ FoM.
$$

Key solver capabilities include:
1. **Four Inflow Models**:
   - Classic Uniform Inflow (Glauert, 1926);
   - Simple Coleman Inflow (Coleman, Feingold & Stempin, 1945);
   - NDARC Coleman-Feingold Inflow (Wayne Johnson, NASA/TP-2009-215402) with longitudinal gradient $K_x$ and lateral gradient $K_y$;
   - Drees Inflow (Drees, 1949) with longitudinal and lateral inflow gradients.
2. **Three Solidity Definitions and Blade Taper**:
   - $\sigma_{\mathrm{ref}}$: Reference solidity extrapolated to the hub ($r=0$);
   - $\sigma_{\mathrm{geom}}$: True geometric physical solidity (actual blade area divided by disk area);
   - $\sigma_{\mathrm{thrust}}$: Thrust-weighted equivalent solidity ($r^2$-weighted).
3. **Blade Pitch, Linear Twist, and Hover Trim**:
   - Constant collective pitch $\theta_0$;
   - Linear twist defined by $\theta_{\mathrm{root}}$ (at root cutout $x_0$) and $\theta_{\mathrm{tip}}$ (at tip $x=1.0$);
   - Hover trim by **collective** (adjusts collective while keeping total twist $\Delta\theta = \theta_{\mathrm{tip}} - \theta_{\mathrm{root}}$ strictly constant);
   - Hover trim by **RPM** (adjusts rotational speed to match target thrust in Newtons).
4. **Advanced Aerodynamic Corrections**:
   - Blade tip loss factor $B$ (fixed or Sissingh self-adjusting);
   - Prandtl-Glauert compressibility correction $a(M)$ on section lift curve slope.
5. **Global Performance Metrics**:
   - Total aerodynamic power against the air $C_{P,\mathrm{air}} = C_Q + \mu C_H$;
   - Effective rotor lift-to-drag ratio $(L/D)_{\mathrm{eff}} = \mu C_T / C_{P,\mathrm{air}}$;
   - Hover Figure of Merit ($FoM$).

---

## 2. Coordinate System and Sign Conventions

### 2.1 Disk Angles and Axial Velocity
Following standard helicopter flight mechanics conventions:
* $\alpha$: rotor disk angle of attack (positive for wind from below, propulsive tilt forward);
* $\mu_z$: non-dimensional axial velocity component through the rotor hub:
  $$\boxed{\mu_z = -\mu\tan\alpha}$$
* For vertical velocity input $w$ [m/s]:
  $$\mu_z = \frac{w}{V_{\mathrm{tip}}} = \frac{w}{\Omega R}$$

The total non-dimensional axial inflow velocity through the rotor disk is:
$$\boxed{\lambda = \mu_z + \lambda_i}$$
where $\lambda_i \ge 0$ is the mean induced downwash velocity.

### 2.2 Moments and Axes
Moments follow standard aircraft body axes:
* $+C_{Mx}$: rolling moment, **right wing down** (advancing blade down);
* $+C_{My}$: pitching moment, **nose up**.

### 2.3 Non-Dimensional Coefficients
$$\mu = \frac{V_h}{\Omega R}, \qquad \mu_z = \frac{V_z}{\Omega R}, \qquad x = \bar{r} = \frac{r}{R}$$
$$C_T = \frac{T}{\rho A (\Omega R)^2}, \qquad C_Q = \frac{Q}{\rho A (\Omega R)^2 R}, \qquad C_H = \frac{H}{\rho A (\Omega R)^2}$$
$$C_Y = \frac{Y}{\rho A (\Omega R)^2}, \qquad C_{Mx} = \frac{M_x}{\rho A (\Omega R)^2 R}, \qquad C_{My} = \frac{M_y}{\rho A (\Omega R)^2 R}$$
where $A = \pi R^2$ is the rotor disk area.

---

## 3. Rotor Solidities

The first aerodynamic lifting station begins at $x_0 = \bar{r}_0 = r_0 / R$ (blade root cutout).
For a tapered blade, the chord distribution between $x_0$ and $1.0$ is:
$$c(x) = c_{\mathrm{root}} + (c_{\mathrm{tip}} - c_{\mathrm{root}})\frac{x - x_0}{1 - x_0}$$

The local non-dimensional solidity is:
$$\sigma(x) = \frac{N c(x)}{\pi R} = s_0 + s_1 x$$
where:
$$\sigma_{\mathrm{root}} = \frac{N c_{\mathrm{root}}}{\pi R}, \qquad \sigma_{\mathrm{tip}} = \frac{N c_{\mathrm{tip}}}{\pi R}$$
$$s_1 = \frac{\sigma_{\mathrm{tip}} - \sigma_{\mathrm{root}}}{1 - x_0}, \qquad s_0 = \sigma_{\mathrm{root}} - s_1 x_0$$

The solver evaluates three complementary solidities:

### 3.1 Reference Solidity ($\sigma_{\mathrm{ref}}$)
Extrapolates the blade chord line to the hub center ($x=0$):
$$c(0) = c_{\mathrm{root}} - \left(\frac{c_{\mathrm{tip}} - c_{\mathrm{root}}}{1 - x_0}\right) x_0$$
$$\boxed{\sigma_{\mathrm{ref}} = \frac{A_{\mathrm{extrapolated}}}{\pi R^2} = \frac{N}{\pi R} \int_0^1 c(x)\,dx = s_0 + \frac{s_1}{2}}$$
*For rectangular blades of constant chord $c$, $\sigma_{\mathrm{ref}} = \frac{N c}{\pi R}$.*

### 3.2 True Geometric Physical Solidity ($\sigma_{\mathrm{geom}}$)
Ratio of actual lifting surface area (integrated strictly from $x_0$ to $1$) to disk area:
$$\boxed{\sigma_{\mathrm{geom}} = \frac{A_{\mathrm{physical}}}{\pi R^2} = \frac{N (1 - x_0)}{\pi R} \left[\frac{c_{\mathrm{root}} + c_{\mathrm{tip}}}{2}\right] = (1 - x_0)\left[s_0 + \frac{s_1 (1 + x_0)}{2}\right]}$$
*For rectangular blades of constant chord $c$, $\sigma_{\mathrm{geom}} = (1 - x_0)\sigma_{\mathrm{ref}}$. With $x_0 = 0.15$, $\sigma_{\mathrm{geom}} = 0.85\,\sigma_{\mathrm{ref}}$.*

### 3.3 Thrust-Weighted Solidity ($\sigma_{\mathrm{thrust}}$)
Equivalent solidity weighted by $x^2$ (Wayne Johnson, *Rotorcraft Aeromechanics*, Eq. 3.122):
$$\boxed{\sigma_{\mathrm{thrust}} = 3 \int_{x_0}^1 \sigma(x) x^2\,dx = 3 \left(s_0 J_2 + s_1 J_3\right) = s_0(1 - x_0^3) + \frac{3}{4}s_1(1 - x_0^4)}$$

### 3.4 Why $\sigma_{\mathrm{ref}}$ Enters the BET Equations Rather than $\sigma_{\mathrm{geom}}$
In Blade Element Theory, the differential lift on a section $dr$ is:
$$dL = \frac{1}{2}\rho V^2 c(r) C_l\,dr$$
Dividing by $\rho A (\Omega R)^2$ and summing over $N$ blades, the factor that scales section loading is the **local solidity**:
$$\sigma(x) = \frac{N c(x)}{\pi R}$$
For constant-chord blades, the local solidity along the entire lifting blade span (from $x_0$ to $1$) is constant and identically equal to the reference solidity:
$$\sigma(x) = \frac{N c}{\pi R} = \sigma_{\mathrm{ref}}$$
The effect of the root cutout $x_0$ is already accounted for **in the integration limits**:
$$\int_{x_0}^1 x^m\,dx = J_m = \frac{1 - x_0^{m+1}}{m+1}$$
If $\sigma_{\mathrm{geom}} = (1 - x_0)\sigma_{\mathrm{ref}}$ were substituted into the canonical integral, the root cutout penalty would be **counted twice**: once in the solidity and once in the integration limits. Therefore:
* **In BET equations**: local solidity $\sigma(x)$ is used (whose level is $\sigma_{\mathrm{ref}}$ for rectangular blades);
* **$\sigma_{\mathrm{geom}}$**: represents physical blade area accounting, but should not multiply integrals that already start from $x_0$.

---

## 4. Inflow Models and Spatial Gradients

The local induced velocity over the rotor disk is modeled with the first-order linear harmonic distribution:
$$\boxed{\lambda_d(x, \psi) = \lambda + x \left(\lambda_{1c}\cos\psi + \lambda_{1s}\sin\psi\right)}$$
where:
$$\lambda = \mu_z + \lambda_i, \qquad \lambda_{1c} = K_x \lambda_i, \qquad \lambda_{1s} = Ky \lambda_i$$
and the wake skew angle $\chi$ is given by:
$$\tan\left(\frac{\chi}{2}\right) = \frac{\mu}{\sqrt{\mu^2 + \lambda^2} + |\lambda|}$$

### 4.1 Uniform Inflow (Glauert, 1926)
$$K_x = 0, \qquad K_y = 0 \implies \lambda_d(x, \psi) = \lambda$$

### 4.2 Simple Coleman (Coleman, Feingold & Stempin, 1945)
$$K_x = \tan\left(\frac{\chi}{2}\right), \qquad K_y = 0$$

### 4.3 NDARC Coleman-Feingold (Wayne Johnson, NASA/TP-2009-215402)
Incorporates the $15\pi/32$ vortex cylinder harmonic factor longitudinally and the lateral inflow gradient:
$$\boxed{K_x = f_x \frac{15\pi}{32} \tan\left(\frac{\chi}{2}\right), \qquad K_y = -f_y 2\mu}$$
In standard NDARC configuration, $f_x = 1.0$ and $f_y = 1.0$.

### 4.4 Drees (Drees, 1949; Johnson, 2013)
$$\boxed{K_x = \frac{4}{3}\left(1 - 1.8\mu^2\right) \tan\left(\frac{\chi}{2}\right), \qquad K_y = -2\mu}$$

### 4.5 Blade Tip Loss Factor ($B$)
To model vortex shedding and vanishing circulation at the blade tip, the effective blade radius is reduced to $B R$, where $B \in (x_0, 1.0]$:
* **Disabled (`TIP_LOSS_MODE = "none"`):** $B = 1.0$.
* **Prescribed (`TIP_LOSS_MODE = "fixed"`):** $B = \text{TIP\_LOSS\_B}$ (typically $0.96$ to $0.98$).
* **Sissingh Formula (`TIP_LOSS_MODE = "sissingh"`):**
  $$B = 1 - \frac{\sqrt{2 C_T}}{N}$$
The active blade span is integrated from $x_0$ to $B$, and the effective momentum disk area is scaled by $B^2$:
$$C_T = 2 B^2 \lambda_i \sqrt{\mu^2 + \lambda^2}$$

### 4.6 Prandtl-Glauert Compressibility Correction
Captures the subsonic compressibility effect on section lift curve slope $a = dC_l / d\alpha$:
* **Disabled (`USE_PRANDTL_GLAUERT = False`):** $a = a_0 = \text{A\_LIFT}$.
* **Enabled (`USE_PRANDTL_GLAUERT = True`):**
  $$M_{\mathrm{tip}} = \frac{\Omega R}{a_{\mathrm{sound}}}, \qquad M_{\mathrm{eff}} = M_{\mathrm{tip}} \sqrt{0.75^2 + \frac{\mu^2}{2}}$$
  $$\boxed{a(M) = \frac{a_0}{\sqrt{1 - \min(M_{\mathrm{eff}}^2, 0.85^2)}}}$$

---

## 5. Hover Trim Modes (`HOVER_TRIM_MODE`)

Hover trim is performed strictly at hover ($\mu = 0, \alpha = 0$). The trimmed parameters are then frozen across the forward flight sweep:

### 5.1 Collective Trim (`HOVER_TRIM_MODE = "collective"`)
* Rotor speed (RPM) is held constant.
* Target can be non-dimensional $C_T$ (`CT_HOVER_TARGET`) or dimensional thrust in Newtons (`THRUST_HOVER_N`).
* **Constant pitch:** solves for exact collective $\theta_0$.
* **Linear twist ($\theta_{\mathrm{root}}, \theta_{\mathrm{tip}}$):** total twist $\Delta\theta = \theta_{\mathrm{tip}} - \theta_{\mathrm{root}}$ is kept **strictly constant**. The solver adds a uniform collective offset $\Delta\theta_0$:
  $$\theta(x) = \left(\theta_{\mathrm{root}} + \Delta\theta_0\right) + \left(\frac{\theta_{\mathrm{tip}} - \theta_{\mathrm{root}}}{1 - x_0}\right)(x - x_0)$$

### 5.2 RPM Trim (`HOVER_TRIM_MODE = "rpm"`)
* Blade pitch angles and twist are held fixed.
* Required target is thrust in Newtons (`THRUST_HOVER_N`, e.g. aircraft gross weight $W$).
* The solver calculates the hover $C_{T,\mathrm{hover}}$ produced by the prescribed pitch and adjusts rotor RPM:
  $$V_{\mathrm{tip}} = \sqrt{\frac{T}{\rho A C_{T,\mathrm{hover}}}} \implies \Omega = \frac{V_{\mathrm{tip}}}{R} \implies \mathrm{RPM} = \frac{60\Omega}{2\pi}$$

---

## 6. Analytical Formulation via Weighted Radial Moments

To unify the analytical formulation with or without tip loss $B$:

$$J_n = \int_{x_0}^B x^n\,dx = \frac{B^{n+1} - x_0^{n+1}}{n+1}$$
$$I_m = \int_{x_0}^B \sigma(x) x^m\,dx = s_0 J_m + s_1 J_{m+1}$$
$$T_m = \int_{x_0}^B \sigma(x)\theta(x) x^m\,dx = p_0 J_m + p_1 J_{m+1} + p_2 J_{m+2}$$
where $p_0 = s_0 t_0$, $p_1 = s_0 t_1 + s_1 t_0$, and $p_2 = s_1 t_1$.

### 6.1 Rotor Forces and Moments
* **Thrust Coefficient $C_T$:**
  $$\boxed{C_T = \frac{a}{2}\left[ T_2 + \frac{\mu^2}{2}T_0 - \left(\lambda + \frac{1}{2}\mu\lambda_{1s}\right)I_1 \right]}$$
* **Longitudinal H-Force $C_H = C_{Hi} + C_{H0}$:**
  $$C_{Hi} = \frac{a}{4}\left[\lambda\mu T_0 + \lambda_{1s} T_2\right]$$
  $$C_{H0}^{\mathrm{comp}} = \frac{C_{d0}}{2}\int_{x_0}^1 \sigma(x) \langle W (x\sin\psi + \mu) \rangle_\psi\,dx, \qquad C_{H0}^{\mathrm{simp}} = \frac{C_{d0}\mu I_1}{2}$$
* **Side Force $C_Y$:**
  $$\boxed{C_Y = -\frac{a\lambda_{1c}}{4}\left[T_2 - 2\lambda I_1\right]}$$
* **Rolling Moment $C_{Mx}$** (+Mx right wing down):
  $$\boxed{C_{Mx} = -\frac{a\mu}{2}\left[T_2 - \frac{1}{2}\lambda I_1\right] + \frac{a\lambda_{1s}}{4}I_3}$$
* **Pitching Moment $C_{My}$** (+My nose up):
  $$\boxed{C_{My} = +\frac{a\lambda_{1c}}{4}I_3}$$
* **Rotor Torque $C_Q = C_{Qi} + C_{Q0}$:**
  $$C_{Qi}^{\mathrm{comp}} = K_{\mathrm{ind}}\lambda_i C_T + \mu_z C_T - \mu C_{Hi}$$
  $$C_{Q0}^{\mathrm{comp}} = \frac{C_{d0}}{2}\int_{x_0}^1 \sigma(x) \langle W u_T x \rangle_\psi\,dx, \qquad C_{Q0}^{\mathrm{simp}} = \frac{C_{d0}}{2}\left(I_3 + \frac{\mu^2 I_1}{2}\right)$$

---

## 7. Global Rotor Performance Metrics

### 7.1 Hover Figure of Merit ($FoM$)
Compares actual hover power to ideal Froude momentum power:
$$\boxed{FoM = \frac{C_T^{3/2}}{\sqrt{2}\,C_Q}}$$

### 7.2 Total Aerodynamic Power Against the Air ($C_{P,\mathrm{air}}$)
Total aerodynamic power imparted to the free stream, combining rotor shaft power ($P_{\mathrm{shaft}} = Q\Omega$) and translational power overcoming longitudinal rotor drag ($P_{\mathrm{trans}} = H V_x$ with $V_x = \mu\Omega R$):
$$P_{\mathrm{total}} = Q\Omega + H V_x = \rho A (\Omega R)^3 \left(C_Q + \mu C_H\right)$$
In non-dimensional form:
$$\boxed{C_{P,\mathrm{air}} = C_Q + \mu C_H}$$

### 7.3 Effective Rotor Lift-to-Drag Ratio ($(L/D)_{\mathrm{eff}}$)
In forward flight ($\mu > 0$), the ratio of forward transport work rate $T V_x$ to total aerodynamic power $P_{\mathrm{total}}$ defines rotor transport efficiency:
$$\boxed{\left(\frac{L}{D}\right)_{\mathrm{eff}} = \frac{\mu C_T}{C_{P,\mathrm{air}}} = \frac{\mu C_T}{C_Q + \mu C_H}}$$

---

## 8. Summary Table of Analytical Formulas

| Output | Generalized Formulation | Canonical Baseline ($\sigma=\text{const}$, $\theta=\theta_0$, $K_y=0$) |
| :--- | :--- | :--- |
| **$C_T$** | $\displaystyle \frac{a}{2}\left[T_2 + \frac{\mu^2}{2}T_0 - \left(\lambda + \frac{1}{2}\mu\lambda_{1s}\right)I_1\right]$ | $\displaystyle \frac{\sigma a}{2}\left[\theta_0\left(J_2 + \frac{\mu^2}{2}J_0\right) - \lambda J_1\right]$ |
| **$C_{Qi}^{\mathrm{comp}}$** | $\displaystyle K_{\mathrm{ind}}\lambda_i C_T + \mu_z C_T - \mu C_{Hi}$ | $\displaystyle K_{\mathrm{ind}}\lambda_i C_T + \mu_z C_T - \mu C_{Hi}$ |
| **$C_{Qi}^{\mathrm{simp}}$** | $\displaystyle \frac{a}{2}\left[\left(\lambda + \frac{\mu\lambda_{1s}}{2}\right)T_2 - \lambda^2 I_1 - \frac{\lambda_{1c}^2+\lambda_{1s}^2}{2}I_3\right]$ | $\displaystyle \frac{\sigma a}{2}\left[\theta_0\lambda J_2 - \lambda^2 J_1 - \frac{\lambda_{1c}^2}{2}J_3\right]$ |
| **$C_{Q0}^{\mathrm{simp}}$** | $\displaystyle \frac{C_{d0}}{2}\left(I_3 + \frac{\mu^2 I_1}{2}\right)$ | $\displaystyle \frac{\sigma C_{d0}}{2}\left(J_3 + \frac{\mu^2 J_1}{2}\right)$ |
| **$C_{Hi}$** | $\displaystyle \frac{a}{4}\left[\lambda\mu T_0 + \lambda_{1s} T_2\right]$ | $\displaystyle \frac{\sigma a\theta_0\lambda\mu J_0}{4}$ |
| **$C_{H0}^{\mathrm{simp}}$** | $\displaystyle \frac{C_{d0}\mu I_1}{2}$ | $\displaystyle \frac{\sigma C_{d0}\mu J_1}{2}$ |
| **$C_Y$** | $\displaystyle -\frac{a\lambda_{1c}}{4}\left[T_2 - 2\lambda I_1\right]$ | $\displaystyle -\frac{\sigma a\lambda_{1c}}{4}\left(\theta_0 J_2 - 2\lambda J_1\right)$ |
| **$C_{Mx}$** | $\displaystyle -\frac{a\mu}{2}\left[T_2 - \frac{\lambda I_1}{2}\right] + \frac{a\lambda_{1s}}{4}I_3$ | $\displaystyle -\frac{\sigma a\mu}{2}\left(\theta_0 J_2 - \frac{\lambda J_1}{2}\right)$ |
| **$C_{My}$** | $\displaystyle +\frac{a\lambda_{1c}}{4}I_3$ | $\displaystyle +\frac{\sigma a\lambda_{1c}J_3}{4}$ |
| **$C_{P,\mathrm{air}}$** | $\displaystyle C_Q + \mu C_H$ | $\displaystyle C_Q + \mu C_H$ |
| **$(L/D)_{\mathrm{eff}}$** | $\displaystyle \frac{\mu C_T}{C_{P,\mathrm{air}}} = \frac{\mu C_T}{C_Q + \mu C_H}$ | $\displaystyle \frac{\mu C_T}{C_Q + \mu C_H}$ |
| **$FoM$** | $\displaystyle \frac{C_T^{3/2}}{\sqrt{2}\,C_Q}$ | $\displaystyle \frac{C_T^{3/2}}{\sqrt{2}\,C_Q}$ |

---

## 9. References

1. **Glauert, H.**, *The Analysis of Experimental Investigations of the Aerodynamics of Flying Models*, R&M 1026, 1926.
2. **Coleman, R. P., Feingold, A. M., Stempin, C. W.**, *Evaluation of the Induced-Velocity Field of an Idealized Helicopter Rotor*, NACA ARR L5E10 / WR L-101, 1945.
3. **Drees, J. M.**, *A Theory of Airflow Through Rotors and Its Application to Some Helicopter Problems*, Journal of the Helicopter Association of Great Britain, Vol. 3, No. 2, 1949.
4. **Johnson, Wayne**, *Helicopter Theory*, Princeton University Press, 1980.
5. **Johnson, Wayne**, *Rotorcraft Aeromechanics*, Cambridge University Press, 2013.
6. **Johnson, Wayne**, *NDARC — NASA Design and Analysis of Rotorcraft: Theory*, NASA/TP-2009-215402, 2009.
7. **Leishman, J. Gordon**, *Principles of Helicopter Aerodynamics*, 2nd Edition, Cambridge University Press, 2006.
