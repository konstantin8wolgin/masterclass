---
title: "phi4imaginary code review (second pass)"
date: 2026-05-20
tags:
  - code-review
  - phi4imaginary
  - sign-problem
status: draft
---

# phi4imaginary Code Review — Second Pass

Scope: I re-read `review.md`, the current source files, and the saved checkpoints. The first review is largely correct. This second pass (a) confirms which findings are real, (b) flags one finding where the review.md is itself slightly off, (c) explains the underlying mechanism of each defect more concretely, and (d) adds findings the first review did not list.

I verified each numerical claim against the live checkpoints (`data/*.jld2`), not against memory. Where I disagree with review.md, I cite the formula explicitly.

---

## A. Confirmed findings from review.md

These I checked and stand by. Numbering is from `review.md`.

### A1. Wegner velocity violates odd symmetry by ~10⁶¹ at t=0 (review.md #1, #4)

Verified: `data/wegner.jld2` has `max|φ̇(φ,0) + φ̇(−φ,0)| = 6.5·10⁶¹`. Mechanism:

```
phi_dot(φ,t) = (iπ/2)·exp(S_hat_reg(φ))·∫_φ^L [S_reg − ⟨S_reg⟩]·exp(−S_reg) dφ'
```

The centered integrand integrates to zero on `[-L,L]`, so

```
∫_φ^L (...) = − ∫_{−L}^φ (...)
```

Either side has the same answer in exact arithmetic. In floating point, they differ. For `φ ≪ 0` the right-side tail integration spans the entire bulk of `exp(-S_reg)` and accumulates a small numerical residue; the prefactor `exp(S_hat_reg(φ))` then amplifies it. At `t=0`, `exp(S_hat_reg(-5)) ≈ exp(168) ≈ 10⁷³`, multiplied by a residue of `~10⁻¹²` gives `~10⁶¹`. The Tikhonov regulator (`ε·φ²/2` with `ε=10⁻³`) adds only `+0.0125` to the prefactor — completely irrelevant.

**Correct fix:** pick the shorter integration side at each `φ` (`∫_{−L}^φ` if `φ < 0`, else `∫_φ^L`) and rely on the centering identity. Even better, evaluate `phi_dot` as a symmetric difference of two one-sided tails and pick the one with the smaller `exp(S)` prefactor cancellation. The current `EPS_REG` band-aid does not address the underlying floating-point cancellation.

### A2. Wegner prefactor inconsistent with the regulated action (review.md #2)

The bare action satisfies `dS_bare/dt = (-iπ/2)·S_bare`; the Wegner formula
`φ̇ = (iπ/2)·exp(S)·∫_φ^∞ (S − ⟨S⟩) exp(-S) dφ'` is derived from that identity. If the integrand uses `S_reg = S_bare + ε·φ²/2` (Tikhonov), the corresponding evolution PDE is now

```
∂_t ρ_reg = − ∂_φ [ φ̇_reg · ρ_reg ],  with  φ̇_reg derived from  dS_reg/dt
                                                                  = dS_bare/dt   (the Tikhonov term has no t-dependence)
```

i.e. `dS_reg/dt = (-iπ/2)·S_bare ≠ (-iπ/2)·S_reg`. The current code keeps the `(iπ/2)` prefactor but replaces `S → S_reg` inside the integral. The velocity does not solve the Wegner equation for either `S_bare` or `S_reg`. The downstream flow inherits this inconsistency.

**Correct fix:** either (i) keep `S_bare` everywhere in `02_wegner.jl` and regulate only the integration domain / the `04_observables.jl` reweighting, or (ii) replace the prefactor with `exp(S_bare)·exp(-iπt/2)/(−iπ/2)·... · dS_bare/dt`-aware expression that consistently uses the actual `dS/dt`. The cleanest path is (i): the regulator is a quadrature device, not a physics modification.

### A3. RK4 evaluates kernel on `real(φ)` only (review.md #3)

`03_flow.jl:18-23` projects every RK4 stage onto the real axis. After one RK4 step a sample is at complex `φ`; the next stage should evaluate `φ̇` at that complex point, not at its real part. Consequence: the imaginary excursions are entirely "free", driven only by the imaginary part of `φ̇` along the real axis. The trajectory shape is dictated by `φ̇_t(Re[φ])` and a passive `Im[φ]` accumulator, which is not what the complex Wegner flow specifies.

**Correct fix:** compute `phi_dot` on a 2D complex grid (slow but honest), or use the holomorphy of `S_hat` to analytically continue `phi_dot(φ_R) → phi_dot(φ_R + i·φ_I)`. For a polynomial action the analytic continuation is by Taylor expansion in `φ_I`; concretely `phi_dot` is itself an integral of holomorphic data so its holomorphic extension is given by deforming the integration contour. The cheapest workable shortcut is a 2D grid of `phi_dot` over `Re[φ] × Im[φ]` at each `t`.

### A4. "PIKfold" label without solving the PIKfold ODE (review.md #5)

`03_flow.jl` labels endpoints "PIKfold". The PIKfold is defined by an ODE in `φ`-space (steepest descent / stationary phase contour). The trajectories shown are `t`-integrations of `φ̇_t` from real initial conditions; they form *Wegner trajectories*, which converge to a PIK-like contour only when the flow is consistent (it is not — see A1, A2). The plot title should say "Wegner trajectories" until the PIKfold ODE is solved as a separate validation.

### A5. Free-theory thimble overlay uses the wrong angle (review.md #6) — **confirmed; this is my mistake**

`03_flow.jl:75` computes `thimble_pts = thimble_phi0 .* exp(im*π/2*t_end)`. The free-theory Wegner velocity is `φ̇_free = (iπ/4)·φ`, whose solution is

```
φ(t) = φ₀ · exp(iπt/4)
```

so the analytic free-theory PIKfold at `t = t_end` is `φ₀ · exp(iπ·t_end/4)`, rotated by π/4 — *half* the angle the plot draws. The endpoint scatter (red stars) lands at roughly the right angle; the dashed reference line is twice as steep, which is exactly why they don't agree on the figure.

**Fix:** change `exp(im*π/2*t_end)` to `exp(im*π/4*t_end)` on line 75 of `03_flow.jl`.

### A6. Observables don't use transported samples (review.md #7)

`04_observables.jl` and `05_summary.jl` never read `data/flow.jld2`. The observables are a real-axis quadrature of `φⁿ·exp(-S_reg)/∫exp(-S_reg)`. The PIK estimator is

```
⟨O⟩(t) = (1/N) · Σ_n O(φ_n(t)) · J_n(t),     φ_n(t) = transport of φ_n(0)
```

with Jacobian `J_n = ∂φ_n(t)/∂φ_n(0)`. None of this is computed. The current pipeline is a regularised real-axis benchmark; it is not a test of the transported-sample PIK estimator.

**Fix:** add an `08_pik_observables.jl` that loads `flow.jld2`, computes `J_n(t)` by finite differences over `phi0_values`, and reports `⟨φ²⟩(t)` as a weighted sum on the transported samples.

### A7. Reported `t = 1` rows are actually `t = 0.99` (review.md #8)

`04_observables.jl` writes `t_grid = range(0.0, 0.99, length=100)`. `05_summary.jl` searches for `argmin(abs.(t_obs .- 1.0))`, finds index 100 (`t=0.99`), and prints the corresponding row. The table prints `t=0.99` numerically but the request was for `t=1.0`. This is misleading in the printed summary. Either drop the `t=1.0` row from the request list or extend the observables to genuinely reach `t=1` (which requires a regulator that handles `cos(πt/2)=0`).

### A8. Free-theory thimble called a "thimble" for `LAM=1` (review.md #10, partial #11)

`07_regulator_comparison.jl` calls the contour `φ_c = φ₀·exp(iπt/4)` a "thimble". For the free theory at this rotation the quadratic action becomes `M²·φ₀²/2` (real positive). For the interacting theory the quartic on this contour is `LAM·φ₀⁴/4·exp(iπt/2)` (Re = `cos(πt/2)·LAM·φ₀⁴/4`), which is *not* stationary on the thimble of the interacting model. The genuine interacting thimble is a curve, not a constant-angle ray. The script should call this an "approximate / free-theory contour" everywhere, not a "thimble".

### A9. Diagnostics plot `S_hat` (bare) while pipeline runs on `S_hat_reg` (review.md #18)

Confirmed at `06_diagnostics.jl:30-31, 42-43`. The action landscape and weight landscape use `S_hat(φ,t)` (no regulator), but `02_wegner.jl` and `04_observables.jl` use `S_reg`. Users inspecting the landscape see the bare problem but the velocity field shown next to it is the regulated one. Should either plot both side-by-side or label clearly.

### A10. CWD-dependent paths in 03, 04, 05, 06, 07 (review.md #16)

Confirmed. `01_action_check.jl` and `02_wegner.jl` use `joinpath(@__DIR__, "data", ...)`; later scripts use `"data/..."`. If a user runs `julia --project=phi4imaginary phi4imaginary/03_flow.jl`, the load fails. Trivial fix: use `@__DIR__` consistently.

### A11. Quiver normalized by global max (review.md #19)

`06_diagnostics.jl:60` sets `scale = 0.25 / max(maxv, 1e-12)` where `maxv = maximum(abs.(phi_dot_grid))`. Since `maxv ≈ 6.5·10⁶¹` (corrupted boundary), every arrow in the bulk has length `~0` and the figure is meaningless. Should either filter to a bulk region (e.g. `|φ| ≤ 3`) before computing `maxv`, or use a per-arrow normalization.

### A12. No `phi=0` grid point (review.md #14)

`N_phi=120` on `[-5,5]` → step `10/119`, no exact zero. Trivial fix: `N_phi=121`.

### A13. No automated invariant tests (review.md #15)

There are scripts named `test_*.jl` but no Pkg test target and no `using Test`. The free-theory and odd-symmetry checks should run as `Pkg.test()` and fail the process on regression.

---

## B. Where review.md is itself slightly off

### B1. review.md #11 — the comment in `07_regulator_comparison.jl` is **self-consistent**, not inconsistent

`review.md` claim:
> "the quartic phase is `exp(+iπt/2)`, but multiplied by the action prefactor `exp(-iπt/2)`, so the quartic term becomes real positive, not `exp(iπt/2)`"

This is wrong. Let me compute directly. With

```
φ_c = φ₀ · exp(iπt/4)
S_bare(φ_c, t) = exp(-iπt/2) · (M²·φ_c²/2 + λ·φ_c⁴/4)
              = exp(-iπt/2) · (M²·φ₀²·exp(iπt/2)/2 + λ·φ₀⁴·exp(iπt)/4)
              = M²·φ₀²/2   +   λ·φ₀⁴·exp(iπt/2)/4
```

So the quartic *does* become `λ·φ₀⁴·exp(iπt/2)/4`, exactly as the code's comment says. Its real part is `cos(πt/2)·λ·φ₀⁴/4`, which is `≥ 0` for `t∈[0,1]` — also as the comment says. The reviewer's "becomes real positive" would be true only for the contour `φ_c = φ₀·exp(iπt/8)`, which is a different rotation and is not what the script uses.

The comment is fine. The deeper critique (review.md #10 — calling this a "thimble" at all) stands; the formula-level claim of inconsistency does not.

### B2. review.md #21 — `abs(φ)` non-holomorphic critique is valid but the framing is too strong

`abs(φ)` under a complex weight is not holomorphic and is not Cauchy-equivalent across deformed contours. However, `⟨|φ|⟩` is a legitimate physical observable in its own right when the question is "what is the typical field magnitude on the integration contour" — it just isn't a contour-invariant quantity. The fix is to label it as "contour-dependent diagnostic", not to delete it. The review.md text suggests it shouldn't be there at all; that's overstated.

---

## C. New findings not in review.md

### C1. The `dS/dt = (-iπ/2)·S` identity is incompatible with the Tikhonov regulator — but the formula is still used

(Related to A2 but more concrete.) The first script `01_action_check.jl:9-12` and the master `CLAUDE.md` both state the *exact identity*

```
dS/dt = (-iπ/2) · S
```

for the bare action. This is the *foundation* of the Wegner formula. With `S_reg = S_bare + ε·φ²/2`, the identity becomes

```
dS_reg/dt = (-iπ/2) · (S_reg − ε·φ²/2) = (-iπ/2)·S_reg + (iπ·ε/4)·φ²
```

The extra `(iπ·ε/4)·φ²` term is missing from the velocity formula in `02_wegner.jl`. With `ε=10⁻³` the correction is tiny in the bulk, but it is the formal reason the Tikhonov-regulated `phi_dot` does not actually solve any Wegner equation.

### C2. `04_sign_severity.png` uses the regulated weight but is labelled as the sign problem of the target theory

`04_observables.jl:104-110` reports `|⟨e^{-i Im S}⟩_Q|` where `S = S_reg`. The Tikhonov term in `S_reg` is purely real, so it adds nothing to `Im S` — the *imaginary* part is the same as the bare problem. However the *normalisation* uses `Q = exp(-Re S_reg)` rather than the bare quenched ensemble, so the reported sign curve is the sign problem in the Tikhonov-regulated ensemble, not in the bare ensemble. For `ε=10⁻³` and the parameters here the difference is small, but the label should say "regulated sign problem".

### C3. `t_grid` range mismatches between scripts

- `02_wegner.jl`: `t_grid = range(0.0, 1.0, length=25)` (includes `t=1`).
- `03_flow.jl`: integrates over the wegner `t_grid`, so to `t=1`.
- `04_observables.jl`: `t_grid = range(0.0, 0.99, length=100)`.
- `05_summary.jl`: asks for `t∈{0,0.25,0.5,0.75,1.0}` from the observables grid; silently substitutes `t=0.99`.
- `06_diagnostics.jl`: uses both `t_h = range(0,1,length=120)` (for landscapes, *bare* action) and `t_grid_w` from wegner.
- `07_regulator_comparison.jl`: `range(0,1,length=50)`.

Six scripts, four different grids. Cross-comparison plots silently combine inconsistent ranges.

**Fix:** define `t_grid` once in a shared file (e.g. add to `reg.jl`) and import.

### C4. `flow.jld2` end-time depends on `wegner.jld2` end-time

`03_flow.jl:43-44` uses `t_start = t_grid[1]`, `t_end = t_grid[end]` from the wegner data. If `02_wegner.jl` is changed to end at `t=0.99` (to dodge the t=1 blowup), the flow now ends there too, and the "PIKfold" label silently refers to `t=0.99`. There is no sentinel check that the user requested `t=1`.

### C5. `06_observables_extra.png` panel uses `S_hat` for "analytic rotation" check while pipeline used `S_hat_reg`

`06_diagnostics.jl:97-101` computes `analytic_S = exp(-iπt/2)·S₀`, where `S₀ = Smean_arr[1]` was generated under the regulated action. The analytic identity `dS/dt = (-iπ/2)·S` holds for the *bare* action only. So the "numeric vs analytic rotation" comparison is comparing a regulated numerical result to an identity that does not apply to the regulated action. Within current `ε=10⁻³` it's almost-OK, but the panel's label promises more than it delivers.

### C6. `regulator_comparison.jld2` saves only three regulators' `phi2` and `sign` — `absphi` and the baseline are not saved

`07_regulator_comparison.jl` saves the three regulators' `phi2` arrays but not the `absphi` arrays, and does not save the unregulated baseline. To regenerate the figure or compare with another run, the baseline has to be recomputed. Minor but real engineering issue.

### C7. `L=12` constant for the thimble path was chosen by hand and is not stated in any comment

`07_regulator_comparison.jl` uses `L = 12.0` for `:thimble`. The reason is that on the rotated contour the quadratic damping is `M²·φ₀²/2` (σ=1), so `L=12` is 12σ. This is a different domain from the `L_smooth(t)` used by the other regulators in the same script — the figures are computed on different integration intervals and that fact is not surfaced in the plot legend.

### C8. The "epsilon extrapolation" plot in `06_diagnostics.jl` recomputes from scratch instead of reusing `04_observables.jl` work

`06_diagnostics.jl:113-125` defines its own `phi2_for_eps(t, ε)` and re-runs `60 × 3` integrations. The canonical pipeline already produced `phi2(t)` for `ε=10⁻³` in `data/observables.jld2`; the diagnostics should reuse that point and only sweep the other two ε values. Not a bug, just duplicated work.

### C9. `02_wegner.jl` saves `reg_method` and `eps_reg` to JLD2, but downstream scripts never check them

`02_wegner.jl:56-60` stores `reg_method` and `eps_reg` in `wegner.jld2`. Nothing reads them. A user changing `EPS_REG` in `02` then re-running `04` (without re-running `02`) will get a silently inconsistent comparison. Either drop the metadata (it's misleading dead weight) or have downstream scripts assert it matches their own settings.

### C10. The `:thimble` regulator's variable `ε` is unused

`reg.jl:33-34` treats `:thimble` the same as `:tikhonov` (`S_reg = S_bare + ε·φ²/2`), but `07_regulator_comparison.jl` does its own custom thimble path and never calls `Regulators.S_reg` with `method=:thimble`. The `:thimble` branch of `reg.jl` is dead code — and worse, it adds a Tikhonov term, which is a misleading thing to find when someone reads `reg.jl` to learn what `:thimble` means.

**Fix:** delete the `:thimble` branch from `reg.jl` or move the contour-deformation code into `reg.jl` as a `weight_thimble(φ₀, t; method)` function.

### C11. The verification commands in review.md are correct, but I add three more invariant checks

The first review.md ends with the commands it ran. Adding tests we should also enforce:

```julia
# (a) odd symmetry of phi_dot
using JLD2; d = load("data/wegner.jld2")
@assert maximum(abs.(d["phi_dot"] .+ reverse(d["phi_dot"], dims=2))) < 1e-6  *  maximum(abs.(d["phi_dot"]))

# (b) free-theory test on the regulated grid
# φ̇_free(φ,t) = (iπ/4)·φ for the LAM=0 problem — should be within ~ε in absolute terms
include("01_action_check.jl"); for (φ,t) in [(1.0,0.5)]; @assert abs(phi_dot_free(φ,t) - im*π/4*φ) < 1e-6; end

# (c) the t-end-of-grid of 04_observables matches the requested target
@assert d_obs["t_grid"][end] == 1.0    # currently 0.99 — fails
```

None of these are present in any current script.

---

## D. Recommended fix order (priority-ordered, mine, not review.md's)

1. **A1 (symmetry-preserving integration)**. Pick the shorter tail at each φ. Without this nothing downstream is trustworthy. ~10 lines.
2. **A5 (wrong thimble angle in `03_flow.jl`)**. One-character fix: `π/2 → π/4`. Without this the visual reference is wrong.
3. **A7 (t=0.99 reported as t=1)**. Extend `04_observables.jl` to `t=1` (relying on A1 to keep it finite) *or* rename the table column to "t ≈ 1" with explicit value. ~3 lines.
4. **A2 / C1 (regulator vs Wegner formula consistency)**. Decide: regulator on integrand only (and use bare action in the prefactor), or full regulator everywhere with corrected `dS_reg/dt`. Update the velocity formula accordingly. ~20 lines.
5. **A6 (PIK observables from transported samples)**. Add `08_pik_observables.jl`. This is the actual *test* of the PIK pipeline. ~80 lines.
6. **A3 (complex evaluation of `phi_dot` in RK4)**. Build a 2D-complex velocity grid in `02`. Bigger change — likely a half-day. Until this is done, calling the trajectories "PIKfold" is overselling.
7. **A4 (label PIKfold honestly)**. Until A6 is done, rename plots and titles to "Wegner trajectories". One-line per file.
8. **Engineering polish (A10, A12, A13, C3, C7, C8, C9, C10)**. Add tests, use `@__DIR__` everywhere, deduplicate `t_grid`, drop dead `:thimble` branch, add metadata cross-checks.

## E. Summary judgement

`review.md` is essentially correct. Of its 22 findings I confirm 20, refine 1 (#21), and disagree with 1 (#11 — the formula in the code is right; the reviewer's algebraic objection doesn't go through). The most important real defects are:

1. Tail-cancellation in the Wegner velocity (A1), which propagates into every downstream artifact.
2. Wrong free-theory thimble angle in `03_flow.jl` (A5) — one-character bug, easy to fix, visible on the figure.
3. The observables pipeline does not use transported samples (A6), so the project does not currently test the PIK estimator at all.

Items 2 and 3 are independent fixes; item 1 must be done first or the rest doesn't matter.
