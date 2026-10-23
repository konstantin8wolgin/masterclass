---
title: "phi4imaginary code review"
date: 2026-05-20
tags:
  - code-review
  - phi4imaginary
  - sign-problem
  - pik
status: draft
---

# phi4imaginary Code Review

Scope: current `projects/phi4imaginary` implementation, generated checkpoint data, and the background notes `[[ihssen-2026-solving-sign-problems-piks]]`, `[[wegner-equation]]`, and `[[sign-problem]]`.

The implementation is not yet a reliable implementation of the complex PIK flow for the real-time proxy phi4 model. The biggest issues are that the velocity field is numerically unstable and symmetry-breaking, the complex flow is evaluated using only the real coordinate, and the observable pipeline is mostly a regularised real-axis reweighting calculation rather than a transported-sample PIK estimator.

## Critical Findings

1. **The saved velocity field violates the expected odd symmetry by enormous amounts.**

   For the Z2-even action in `01_action_check.jl:8` and `reg.jl:24`, the Wegner velocity should satisfy `phi_dot(-phi,t) = -phi_dot(phi,t)`. The saved `data/wegner.jld2` does not. Measured from the current checkpoint:

   | t | max abs velocity | max odd-symmetry error |
   |---|---:|---:|
   | 0.0 | `6.506684250584152e61` | `6.506684250584152e61` |
   | 0.5 | `3.2532178281418445e38` | `3.2532178281418445e38` |
   | 1.0 | `4.831101163730034e6` | `5.361688627254067e6` |

   The root cause is the one-sided tail formula in `02_wegner.jl:28-32`: it computes `exp(S_hat_reg(phi,t)) * integral(phi,L,...)`. In the tails, the integral is tiny and cancellation-sensitive while `exp(S)` can be huge, so roundoff becomes a fake velocity. Negative `phi` is especially fragile because the integral crosses nearly the whole symmetric domain and depends on cancellation.

2. **The velocity field is computed from a regularised action but still uses the unregularised Wegner prefactor.**

   `02_wegner.jl:21-32` computes the centered integrand with `S_hat_reg`, but then multiplies by the fixed prefactor `(im*pi/2)` from the bare path derivative. If the action is changed to `S_reg = S_bare + eps*phi^2/2`, then `dS_reg/dt` is not simply `(-im*pi/2) * S_reg`. The Tikhonov term has no `t` derivative. This means the velocity no longer solves the Wegner equation for the action actually being integrated.

3. **The complex RK4 flow discards the imaginary coordinate when evaluating the velocity.**

   `03_flow.jl:15-18` clamps and interpolates only `real(phi_complex)`. Once RK4 moves a point into the complex plane, `k2`, `k3`, and `k4` should evaluate the kernel at the complex field value. Instead, the code projects every intermediate point back onto the real axis. The resulting trajectory is not the complex PIK flow described in `[[wegner-equation]]`.

4. **The flow output is numerically blown up at `t=1`.**

   The current saved `data/flow.jld2` has endpoint `-611.6743903824411 + 39499.76059698362im` for the first trajectory. This is not a plausible transported sample for the intended zero-dimensional phi4 benchmark; it is consistent with integrating the corrupted velocity field from finding 1.

5. **The code labels the trajectory plot as a PIKfold / Picard-Lefschetz thimble without computing the PIKfold condition.**

   `03_flow.jl:58-87` labels endpoints as `PIKfold` and says the flow goes to a `Picard-Lefschetz thimble`. The background note defines the PIKfold by the ODE `d phi_I / d phi_R = -Im[p_t]/Re[p_t]` in `[[ihssen-2026-solving-sign-problems-piks]]`. This ODE is not solved anywhere in `03_flow.jl`. The plotted trajectories are RK4 integrations of the approximate velocity field, not a verified PIKfold.

6. **The thimble overlay uses the wrong angle for the free-theory PIK trajectory.**

   `03_flow.jl:74-79` uses `phi0 * exp(im*pi/2*t_end)`. The wiki reference for the free real-time proxy says the free solution is `phi_t = phi_0 * exp(i*pi*t/4)`, not `exp(i*pi*t/2)`. The overlay is therefore rotated twice as far as the analytic PIK flow.

7. **Observables are not computed from transported samples.**

   `04_observables.jl:19-38` performs direct real-axis quadrature with a regularised complex weight. The reference algorithm in `[[ihssen-2026-solving-sign-problems-piks]]` says observables at `t=1` should be evaluated on transported samples. The current `flow.jld2` is not used by `04_observables.jl`, so the observable pipeline does not test the PIK flow.

8. **The observable run stops at `t=0.99`, but summary output reports it as `t=1.0`.**

   `04_observables.jl:46` sets `t_grid = range(0.0, 0.99, length=N_t)`. `05_summary.jl:26-32` asks for `t_val = 1.0`, uses the nearest point, and prints `0.99` in the table. This is easy to misread as a target-time result. The target theory in the wiki is the endpoint `t=1`.

9. **The real-axis regulator changes the problem being solved.**

   `reg.jl:24-31` adds a real positive damping term or a shifted mass. This may be a useful numerical experiment, but it is not the bare S3.4 action. Since the pipeline saves only one `wegner.jld2` and downstream scripts consume it as the main velocity field, the project mixes the regulated auxiliary problem with the claimed target problem.

10. **The regulator comparison calls the free-theory contour a thimble for the interacting phi4 action.**

   `07_regulator_comparison.jl:13-32` uses the contour `phi_c = phi0 * exp(i*pi*t/4)` and labels it `:thimble` / `free-theory thimble contour`. This is analytic for the free theory, but with `LAM = 1` the quartic term is still present. It is not shown to be the actual Lefschetz thimble or PIKfold for the interacting model.

11. **The thimble contour comment has a wrong quartic phase.**

   `07_regulator_comparison.jl:14-18` says the quartic becomes `LAM * phi0^4 * exp(i*pi*t/2) / 4`. For `S_t(phi) = exp(-i*pi*t/2) * (m2 phi^2/2 + lam phi^4/4)` and `phi = phi0 * exp(i*pi*t/4)`, the quartic phase is `exp(+i*pi*t/2)`, but multiplied by the action prefactor `exp(-i*pi*t/2)`, so the quartic term becomes real positive, not `exp(i*pi*t/2)`. The comment's reasoning is inconsistent with the formula.

12. **The `epsmass` implementation does not match its documented meaning.**

   `reg.jl:15-20` documents `:epsmass` as a Feynman-style prescription; the current code at `reg.jl:26-29` shifts `m2 -> m2 + i*eps` only. The quartic coupling remains unregularised. At `t=1`, the quartic part is still purely oscillatory on the real axis. This may be intentional, but the documentation and numerical interpretation are not precise enough for a benchmark.

13. **Velocity-field plots hide the symmetry failure.**

   `02_wegner.jl:90-118` creates log plots with `abs` and `clamp_pos`, and the log-log panels show only `phi > 0`. Those plots are useful for magnitude, but they cannot diagnose odd symmetry. The linear plot also clamps `ylims=(-3,3)` in `02_wegner.jl:74-77`, hiding the enormous corrupted values present in `wegner.jld2`.

14. **The phi grid does not contain `phi=0`, but tests and comments treat it as if it does.**

   `02_wegner.jl:36-38` uses `N_phi = 120` on `[-5,5]`, so zero is not a grid point. The closest saved value is `-0.04201680672268908`. `test_wegner_quick.jl:48-54` labels the closest grid point as `phi_dot(0,0)`, which creates a false symmetry failure at the center and hides the real tail/cancellation failure.

15. **No automated assertion protects the core analytic invariants.**

   The project has scripts named `test_*`, but they are not part of a Julia `Test` suite and do not fail the process on invariant violations. There is no enforced test for:

   - `phi_dot(-phi,t) + phi_dot(phi,t) ~= 0`
   - the free-theory analytic solution `phi_dot = i*pi*phi/4`
   - finite velocities over the plotted grid
   - convergence of the velocity field with respect to `L`, `rtol`, `N_phi`, and `N_t`
   - agreement between transported-sample observables and a trusted quadrature benchmark

16. **Several scripts depend on the current working directory.**

   `03_flow.jl:4`, `04_observables.jl:6`, `05_summary.jl:4-6`, `06_diagnostics.jl:6`, and `07_regulator_comparison.jl:6` load `data/...` using relative paths. `01_action_check.jl` and `02_wegner.jl` use `@__DIR__`. Running later scripts from outside the project directory can load or write the wrong path.

17. **Generated data is treated as authoritative even when upstream code changes.**

   `05_summary.jl` and `06_diagnostics.jl` load whatever is in `data/*.jld2`. They do not check that `reg_method`, `eps_reg`, grid sizes, or source script settings are mutually consistent. A stale `flow.jld2` can be combined with a fresh `wegner.jld2` or `observables.jld2` without warning.

18. **The diagnostics use the bare action even when the pipeline uses a regularised action.**

   `06_diagnostics.jl:30-43` plots landscapes for `S_hat`, not `S_hat_reg`, while the active velocity and observables use `S_hat_reg`. These figures can make the user inspect the wrong effective problem.

19. **The quiver diagnostic is dominated by corrupted outliers.**

   `06_diagnostics.jl:59-66` scales all arrows by `maximum(abs.(phi_dot_grid))`. Because the saved velocity grid contains extreme tail artifacts, most meaningful arrows are likely flattened.

20. **The free-theory comparison in `04_observables.jl` is not a validation of the interacting run.**

   `04_observables.jl:72-78` compares the interacting observable with the free-theory formula `exp(i*pi*t/2)/M2`. That formula is only the `LAM = 0` check. With `LAM = 1`, disagreement does not directly diagnose correctness, and apparent similarity would not validate the interacting PIK flow.

21. **`abs(phi)` is used as a complex observable in direct reweighting.**

   `04_observables.jl:30` and `07_regulator_comparison.jl:30,42` compute `abs(phi)` or `abs(phi_c)` under complex weights. This is not holomorphic. If the purpose is contour deformation / Cauchy-equivalent observables, non-holomorphic observables need explicit justification and should not be mixed with holomorphic checks like `phi^2`.

22. **The project mixes several interpretations without naming separate outputs.**

   The current folder contains:

   - bare action checks in `01_action_check.jl`
   - regularised real-axis Wegner data in `02_wegner.jl`
   - projected-complex trajectory visualisation in `03_flow.jl`
   - regularised real-axis reweighting in `04_observables.jl`
   - regulator comparison in `07_regulator_comparison.jl`

   These are different experiments, but they all write into the same `data/` namespace and summary figures. The result looks like one coherent pipeline when it is not.

## Recommended Fix Order

1. Separate the bare S3.4 target, each regulator, and each contour experiment into explicitly named outputs.
2. Replace the unstable one-sided tail formula or compute it in a symmetry-preserving way before trusting any velocity plot.
3. Add hard tests for odd symmetry, finite velocities, and the free-theory analytic result.
4. Decide whether the velocity is for the bare action or a regulated action; then use the matching `dS/dt` in the Wegner equation.
5. For complex trajectories, evaluate the kernel on complex `phi`, not only `real(phi)`.
6. Compute the PIKfold ODE explicitly before labeling any plot as a PIKfold.
7. Keep real-axis reweighting, contour quadrature, and transported-sample estimators as separate scripts with separate filenames and legends.

## Verification Commands Run During Review

```bash
julia --project=. -e 'using JLD2; d=load("data/wegner.jld2"); ...'
julia --project=. -e 'using JLD2; d=load("data/observables.jld2"); ...'
julia --project=. -e 'using JLD2; d=load("data/flow.jld2"); ...'
julia --project=. -e 'using JLD2; d=load("data/regulator_comparison.jld2"); ...'
```

Observed checkpoint facts:

- `data/wegner.jld2`: `reg_method=tikhonov`, `eps_reg=0.001`, no exact `phi=0` grid point.
- `data/observables.jld2`: `t` range is `0.0..0.99`, not `0.0..1.0`.
- `data/flow.jld2`: `t` range is `0.0..1.0`, but trajectories are corrupted by the unstable velocity grid.
- `data/regulator_comparison.jld2`: at `t=1`, Tikhonov and epsmass produce `phi2 ~= -13.2805 - 26.5930im`, while the contour experiment gives `phi2 ~= 0.2598 + 0.4376im`; these are different auxiliary calculations, not mutually validating results.
