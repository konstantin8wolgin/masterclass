# AGENTS.md

This file provides guidance to Codex (Codex.ai/code) when working with code in this repository.

## What this repository is

A dual-workspace vault for research in lattice field theory and physics-informed kernels:

| Workspace | When it applies | Root rule |
|-----------|----------------|-----------|
| **Code** (`projects/`) | implementing simulations | all code and data stay inside the relevant project directory |
| **Wiki** (`wiki/`) | ingesting papers, answering queries, maintaining knowledge | LLM creates/updates everything under `wiki/`; never touches `raw/` |

---

## Workspace A: Julia code projects

### Running scripts

Scripts are numbered pipelines; run them in order. Each script saves output to `data/` (JLD2 format) and the next script loads from it.

```bash
# Run from the project directory
julia --project=. 01_action_check.jl
julia --project=. -t auto 02_wegner.jl   # -t auto enables multithreading for Threads.@threads
julia --project=. 03_flow.jl
# etc.
```

First-time setup or after `Project.toml` changes:
```bash
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

Dependencies (both projects): `QuadGK`, `JLD2`, `Plots`, `Interpolations`. No `DifferentialEquations.jl` — RK4 is implemented manually in script 03.

### Project layout

```
projects/<name>/
├── Project.toml / Manifest.toml
├── data/                   JLD2 checkpoints + PNG figures (one file per script)
├── 01_action_check.jl      define S_t, validate analytic identities (free-theory check)
├── 02_wegner.jl            compute φ̇_t on (N_t × N_phi) grid via quadrature, plot
├── 03_flow.jl              RK4-integrate 10 000 samples t=0→1, save snapshots
├── 04_*.jl                 observables / histogram comparison against theory
└── 05_summary.jl           load all data, print statistics table, combined figure
```

Active projects:
- **`projects/phi4imaginary/`** — zero-dimensional phi4 with a complex action (Minkowski rotation); primary reference: Ihssen 2026 Supplement S3.4. Uses complex arithmetic throughout.
- **`projects/linearphi4/`** — phi4 with a real linear action path. No complex numbers.

### Physics reference (phi4imaginary)

Action rotating from Euclidean (t=0) to Minkowski (t=1):
```
S_hat(φ, t) = exp(-iπt/2)/2 · m²·φ² + exp(-iπt/2)/4 · λ·φ⁴
dS_dt       = (-iπ/2) · S_hat    (exact — both couplings carry the same exp factor)
```

Velocity field from Wegner equation:
```
φ̇_t(φ) = (iπ/2) · exp(S_hat(φ,t)) · ∫_φ^∞ [S_hat(φ') - ⟨S_hat⟩] · exp(-S_hat(φ')) dφ'
```

The centering term `⟨S_hat⟩` is essential — without it the free-theory check (`φ̇ = (iπ/4)φ` when λ=0) fails. Analytic validation targets: identity check err < 1e-14, free-theory φ̇ err < 1e-6.

### Physics reference (linearphi4)

Action interpolating from free Gaussian (t=0) to bimodal phi4 (t=1):
```
S_hat(φ, t) = m²(t)/2 · φ² + λ(t)/4 · φ⁴
m²(t) = 1 − 3t,  λ(t) = t/2

dŜ/dt(φ) = −3/2 · φ² + 1/8 · φ⁴      (analytic)

φ̇_t(φ) = exp(Ŝ_t(φ)) · ∫_φ^∞ [dŜ/dt(φ') − ⟨dŜ/dt⟩_t] · exp(−Ŝ_t(φ')) dφ'
```

Free-theory analytic check: `φ̇_t^free = −3/(2(1−3t)) · φ`. Singularity at t=1/3 (m²=0) is regularised by the quartic — quadrature still converges.

### Integration bounds

`integration_bound(t)` computes a finite L that approximates +∞ by ensuring `exp(-Re[S(L,t)]) < 1e-12`. Scripts call `quadgk(..., φ, L; rtol=1e-10)` rather than `quadgk(..., φ, Inf)` for numerical stability. Typical bounds: 8–30 for phi4imaginary (driven by `cos(πt/2)`), 5–50 for linearphi4 (driven by quartic or quadratic term).

### Interpolant dimension order

The velocity field matrix is shaped `(N_t, N_phi)`. Build the interpolant as `(t_grid, phi_grid)` and access as `interp(t, phi)`.

---

## Workspace B: Wiki

### Operational discipline

**On every ingest:**
1. Write or update the source page (`wiki/sources/<firstauthor-year-slug>.md`)
2. Update or create concept / method / equation pages the source touches
3. Update `wiki/index.md` — add or update rows for every page touched
4. Append to `wiki/log.md`: `## [YYYY-MM-DD] ingest | <title> | <N pages updated>`
5. Report which pages were created or updated

**On every query:**
1. Read `wiki/index.md` first
2. Use `tools/search.sh "<query>"` when the index alone is insufficient
3. Read relevant pages in full; answer with citations as `[[page-name]]`

**Never:** modify any file under `raw/`; create duplicate pages; leave `index.md` or `log.md` out of date after an ingest.

### Page frontmatter templates

**Source** (`wiki/sources/<firstauthor-year-slug>.md`):
```yaml
---
title: ""
authors: []
year:
type: paper  # paper | article | note
tags: []
source_file: raw/papers/<filename>
---
```
Sections: Summary, Key claims, Methods used, Notation introduced, Open questions, Links.

**Concept / Method** (`wiki/concepts/` or `wiki/methods/`):
```yaml
---
title: ""
tags: []
sources: []
---
```
Sections: Definition, Intuition, Key results, Connections, Open questions, Sources. Method pages add an Algorithm section.

**Equation** (`wiki/equations/`):
```yaml
---
title: ""
tags: []
sources: []
---
```
Sections: Statement, Variables, Derivation sketch, Limits and special cases, Sources.

**Problem** (`wiki/problems/`):
```yaml
---
title: ""
status: open  # open | in-progress | resolved
date_opened: YYYY-MM-DD
date_resolved:
---
```
Sections: Problem statement, Context, Attempts, Resolution.

### Search tool

```bash
tools/search.sh "<query>"
```

Shells out to `qmd search` (Node.js `@tobilu/qmd`) against `wiki/`. Use when `index.md` alone is insufficient.

### Style rules (wiki pages, log entries, all responses)

1. No emojis anywhere.
2. Clear order: numbered steps, bulleted lists, tables with headers.
3. Concise: one sentence where one sentence is enough.
4. No commentary on own actions — just do it and report the result.
