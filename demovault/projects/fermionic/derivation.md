---
title: Wegner velocity derivation for the one-site flow
tags:
  - fermionic
  - hubbard
  - physics-informed-kernels
  - wegner-equation
---

# Wegner Velocity Derivation

This note derives the one-dimensional velocity formula used to transport a
time-dependent density

$$
p_t(\phi) \propto e^{-S_t(\phi)}
$$

along a flow

$$
\frac{d\phi_t}{dt} = \dot\phi_t(\phi_t).
$$

The same derivation applies to real or holomorphic complexified fields, as
long as the integration contour and boundary terms are well defined.

## Setup

Let the unnormalized density be

$$
w_t(\phi) = e^{-S_t(\phi)}
$$

and the normalized density be

$$
\rho_t(\phi) = \frac{w_t(\phi)}{Z_t},
\qquad
Z_t = \int d\phi\, e^{-S_t(\phi)}.
$$

Expectation values at flow time $t$ are

$$
\langle A\rangle_t
=
\frac{\int d\phi\, A(\phi)e^{-S_t(\phi)}}{\int d\phi\, e^{-S_t(\phi)}}.
$$

The goal is to find a velocity field

$$
v_t(\phi) = \dot\phi_t(\phi)
$$

such that moving samples by this velocity preserves the probability density
specified by $S_t$.

## Continuity Equation

Probability conservation under the deterministic flow gives the one-dimensional
continuity equation

$$
\partial_t \rho_t(\phi)
+
\partial_\phi\!\left[v_t(\phi)\rho_t(\phi)\right]
= 0.
$$

Equivalently,

$$
\partial_\phi\!\left[v_t(\phi)\rho_t(\phi)\right]
=
-\partial_t \rho_t(\phi).
$$

## Time Derivative of the Normalized Density

Start from

$$
\rho_t(\phi)=\frac{e^{-S_t(\phi)}}{Z_t}.
$$

Taking a time derivative gives two terms:

$$
\partial_t\rho_t(\phi)
=
-\partial_t S_t(\phi)\rho_t(\phi)
-
\rho_t(\phi)\,\partial_t\log Z_t.
$$

Now compute the normalization term:

$$
\partial_t Z_t
=
\int d\phi\,[-\partial_t S_t(\phi)]e^{-S_t(\phi)}
=
-Z_t\langle\partial_t S_t\rangle_t.
$$

Therefore

$$
\partial_t\log Z_t
=
\frac{\partial_t Z_t}{Z_t}
=
-\langle\partial_t S_t\rangle_t.
$$

Substitute this back into the derivative of $\rho_t$:

$$
\partial_t\rho_t(\phi)
=
\left[
-\partial_t S_t(\phi)
+
\langle\partial_t S_t\rangle_t
\right]\rho_t(\phi).
$$

This is where the centering term comes from. It is not optional: it is the
normalization correction from differentiating $Z_t$.

## First-Order Equation for the Velocity

Insert the result into the continuity equation:

$$
\partial_\phi\!\left[v_t(\phi)\rho_t(\phi)\right]
=
\left[
\partial_t S_t(\phi)
-
\langle\partial_t S_t\rangle_t
\right]\rho_t(\phi).
$$

Since $\rho_t=e^{-S_t}/Z_t$, the constant $Z_t$ cancels:

$$
\partial_\phi\!\left[v_t(\phi)e^{-S_t(\phi)}\right]
=
\left[
\partial_t S_t(\phi)
-
\langle\partial_t S_t\rangle_t
\right]e^{-S_t(\phi)}.
$$

This is the one-dimensional Wegner equation in divergence form.

Expanding the derivative gives the equivalent local PDE

$$
\partial_\phi v_t(\phi)
-
v_t(\phi)\partial_\phi S_t(\phi)
=
\partial_t S_t(\phi)
-
\langle\partial_t S_t\rangle_t.
$$

Depending on convention, this may also be written with all terms moved to the
other side:

$$
\partial_t S_t
+
v_t\partial_\phi S_t
-
\partial_\phi v_t
-
\langle\partial_t S_t\rangle_t
= 0.
$$

If $S_t$ is defined only up to an additive time-dependent constant, the
centering term removes that constant from the velocity equation.

## Integral Solution

Define

$$
g_t(\phi)
=
\left[
\partial_t S_t(\phi)
-
\langle\partial_t S_t\rangle_t
\right]e^{-S_t(\phi)}.
$$

The divergence equation is

$$
\partial_\phi\!\left[v_t(\phi)e^{-S_t(\phi)}\right]
=
g_t(\phi).
$$

Integrate from a left boundary $-L$ to $\phi$:

$$
v_t(\phi)e^{-S_t(\phi)}
-
v_t(-L)e^{-S_t(-L)}
=
\int_{-L}^{\phi} g_t(\phi')\,d\phi'.
$$

If the boundary term vanishes,

$$
v_t(-L)e^{-S_t(-L)} \approx 0,
$$

then

$$
v_t(\phi)
=
e^{S_t(\phi)}
\int_{-L}^{\phi}
\left[
\partial_t S_t(\phi')
-
\langle\partial_t S_t\rangle_t
\right]
e^{-S_t(\phi')}\,d\phi'.
$$

This is the left-boundary form.

Alternatively, integrate from $\phi$ to a right boundary $L$:

$$
v_t(L)e^{-S_t(L)}
-
v_t(\phi)e^{-S_t(\phi)}
=
\int_{\phi}^{L} g_t(\phi')\,d\phi'.
$$

If the right boundary term vanishes,

$$
v_t(L)e^{-S_t(L)} \approx 0,
$$

then

$$
v_t(\phi)
=
-e^{S_t(\phi)}
\int_{\phi}^{L}
\left[
\partial_t S_t(\phi')
-
\langle\partial_t S_t\rangle_t
\right]
e^{-S_t(\phi')}\,d\phi'.
$$

This is the right-boundary form.

> [!warning] Sign convention
> The upper-limit integral carries a minus sign with the continuity equation
> convention used above. A formula with $\int_\phi^L$ and no minus sign uses
> the opposite sign convention for the velocity or has a sign error.

## Why the Centering Term Is Necessary

The source term must integrate to zero:

$$
\int d\phi\,
\left[
\partial_t S_t(\phi)
-
\langle\partial_t S_t\rangle_t
\right]e^{-S_t(\phi)}
=0.
$$

This follows immediately from the definition of the average. Without this
condition, the divergence equation would inject or remove total probability,
which cannot be represented by a pure transport flow with vanishing boundary
flux.

The centering term also makes the velocity insensitive to time-dependent
additive constants in the action. If

$$
S_t(\phi)\mapsto S_t(\phi)+c(t),
$$

then

$$
\partial_t S_t(\phi)\mapsto \partial_t S_t(\phi)+c'(t),
\qquad
\langle\partial_t S_t\rangle_t
\mapsto
\langle\partial_t S_t\rangle_t+c'(t),
$$

so the difference

$$
\partial_t S_t-\langle\partial_t S_t\rangle_t
$$

is unchanged.

This is why a constant offset such as the Path 1 value

$$
S_{t=0}^{(1)}(\phi)
=
\frac{\beta}{2U}\phi^2 - 2\log 2
$$

does not change the normalized Gaussian measure or the velocity field.

## Boundary Conditions

The integral solution assumes that the probability flux vanishes at the chosen
boundary:

$$
v_t(\phi)e^{-S_t(\phi)} \to 0
\qquad
\text{as } \phi\to\pm\infty.
$$

Numerically, the infinite boundary is replaced by a finite $L$ chosen so that
$e^{-\operatorname{Re} S_t(L)}$ is negligible. In the fermionic one-site flow
notebook, this role is played by `L_BOUND`.

For complex actions, this condition should be interpreted along the integration
contour. The derivation is algebraically the same, but convergence depends on
the real part of the action along that contour.

## Connection to the CDF Formula

The same result can be written using the cumulative distribution

$$
F_t(\phi)=\int_{-L}^{\phi}\rho_t(\phi')\,d\phi'.
$$

Differentiate with respect to $t$:

$$
\partial_t F_t(\phi)
=
\int_{-L}^{\phi}\partial_t\rho_t(\phi')\,d\phi'.
$$

Using the continuity equation,

$$
\partial_t F_t(\phi)
=
-\int_{-L}^{\phi}\partial_{\phi'}
\left[v_t(\phi')\rho_t(\phi')\right]d\phi'.
$$

With vanishing left boundary flux,

$$
\partial_t F_t(\phi)
=
-v_t(\phi)\rho_t(\phi).
$$

Therefore

$$
v_t(\phi)
=
-\frac{\partial_t F_t(\phi)}{\rho_t(\phi)}.
$$

This is the CDF form used in the PIK literature. Substituting the expression
for $\partial_t\rho_t$ reproduces the left-boundary integral solution above.

## Practical Formula Summary

With

$$
\bar S_t'
=
\langle\partial_t S_t\rangle_t,
$$

the two equivalent boundary forms are:

$$
\boxed{
\dot\phi_t(\phi)
=
e^{S_t(\phi)}
\int_{-L}^{\phi}
\left[
\partial_t S_t(\phi')-\bar S_t'
\right]
e^{-S_t(\phi')}\,d\phi'
}
$$

and

$$
\boxed{
\dot\phi_t(\phi)
=
-e^{S_t(\phi)}
\int_{\phi}^{L}
\left[
\partial_t S_t(\phi')-\bar S_t'
\right]
e^{-S_t(\phi')}\,d\phi'
}
$$

These two forms agree when the total centered source integrates to zero and the
boundary fluxes vanish.

