# Demovault — Session Instructions

Read this file at the start of every session before doing anything else.
This vault has two distinct workspaces. Check which one applies to your session.

---

## Workspace A: phi4 project (code)

**Applies when**: working in `projects/phi4/` or asked to implement the phi4 simulation.
**Rule**: all code and data files must stay inside `projects/phi4/`. Do not write files elsewhere.

### What to read first

Read these wiki files in full before writing any code:

- `wiki/sources/ihssen-2026-solving-sign-problems-piks.md` — primary reference;
  Supplement S3.4 is the exact model being implemented
- `wiki/equations/wegner-equation.md` — derivation and closed-form solution
- `wiki/concepts/sign-problem.md` — physics background and reweighting formula
- `wiki/sources/lawrence-2023-deep-learning-fermion-sign.md` — Lawrence & Yamauchi
  context on fermion sign problem (background only)

### Project overview

Zero-dimensional phi4 theory with a complex action (real-time path integral).
The action is rotated from Euclidean (t=0) to Minkowski (t=1) along a one-parameter path.
The sign problem is addressed via the PIK flow (Ihssen 2026, Supplement S3.4).

Full design spec: `../docs/superpowers/specs/2026-05-05-phi4-sign-problem-design.md`

Read that spec before implementing. It contains: physics equations, algorithm, script
responsibilities, numerical parameters, dependencies, and validation checks.

### Quick-reference equations

Action path:
    S_t(phi) = m_t^2/2 * phi^2 + lambda_t/4 * phi^4
    m_t^2 = exp(-i*pi*t/2),  lambda_t = exp(-i*pi*t/2),  t in [0, 1]

Key identity:
    dS_t/dt = (-i*pi/2) * S_t(phi)

Closed-form velocity field (from Wegner equation, normalized):
    phi_dot_t(phi) = (i*pi/2) * exp(S_hat_t(phi)) * integral_{phi}^{inf} [S_hat_t(phi') - <S_hat_t>] * exp(-S_hat_t(phi')) dphi'

where S_hat_t is the field-dependent part (no log Z), and <S_hat_t> = integral S_hat_t e^{-S_hat_t} / integral e^{-S_hat_t}.
The centering term is essential: without it the free-theory check (lambda=0, phi_dot = i*pi/4 * phi) fails.

Observable (reweighting):
    <O>_t = <O * exp(-i * Im[S_t])>_Q / <exp(-i * Im[S_t])>_Q
    where Q_t proportional to exp(-Re[S_t(phi)])

### File layout

    projects/phi4/
    ├── Project.toml
    ├── data/                     (JLD2 checkpoints, one per script)
    ├── 01_action_check.jl        define S_t, validate free-theory analytic result
    ├── 02_wegner.jl              solve Wegner eq on grid, save phi_dot_t, plot
    ├── 03_flow.jl                integrate flow ODE (RK4), save trajectories, plot PIKfold
    ├── 04_observables.jl         quadrature + reweighting, save <O>_t, plot
    └── 05_summary.jl             load all data, print table, combined 2x2 figure

### Language and dependencies

Julia. Dependencies: QuadGK.jl, JLD2.jl, Plots.jl, Interpolations.jl.
No DifferentialEquations.jl — use manual RK4 in script 03.

---

## Workspace B: wiki (knowledge base)

**Applies when**: ingesting papers, answering queries, maintaining wiki pages.
**Rule**: wiki pages live in `wiki/`. Never modify files under `raw/`.

### What this wiki is

A persistent, LLM-maintained knowledge base for research in lattice field theory and
physics-informed kernels. The human curates sources and asks questions. The LLM does all
the filing, cross-referencing, and maintenance.

### Layer boundaries

| Directory | Owner  | Rule                              |
|-----------|--------|-----------------------------------|
| raw/      | Human  | LLM reads, never modifies         |
| wiki/     | LLM    | LLM creates and updates everything|
| tools/    | Shared | Human installs; LLM invokes       |
| CLAUDE.md | Shared | Evolves by mutual agreement       |

### File naming conventions

- Source pages: `wiki/sources/<firstauthor-year-slug>.md`
- All other wiki pages: `wiki/<dir>/<kebab-case-title>.md`
- Raw papers: keep original filename; note it in source page frontmatter

### Page format templates

#### Source page

```yaml
---
title: ""
authors: []
year:
type: paper
tags: []
source_file: raw/papers/<filename>
---
```

Sections (in order): Summary, Key claims, Methods used, Notation introduced,
Open questions, Links.

#### Concept page

```yaml
---
title: ""
tags: []
sources: []
---
```

Sections: Definition, Intuition, Key results, Connections, Open questions, Sources.

#### Method page

Same as concept page plus an Algorithm section.

#### Equation page

```yaml
---
title: ""
tags: []
sources: []
---
```

Sections: Statement, Variables, Derivation sketch, Limits and special cases, Sources.

#### Problem page

```yaml
---
title: ""
status: open
date_opened: YYYY-MM-DD
date_resolved:
---
```

Sections: Problem statement, Context, Attempts, Resolution.

### Style rules

1. No emojis anywhere.
2. Clear order: numbered steps, bulleted lists, tables with headers.
3. Concise: one sentence where one sentence is enough.
4. No commentary on your own actions — just do it and report the result.

### Operational discipline

**On every ingest:**
1. Write or update the source page
2. Update or create any concept / method / equation pages the source touches
3. Update `wiki/index.md` — add or update rows for every page touched
4. Append to `wiki/log.md`: `## [YYYY-MM-DD] ingest | <title> | <N pages updated>`
5. Report which pages were created or updated

**On every query:**
1. Read `wiki/index.md` first to find relevant pages
2. Use `tools/search.sh <query>` when the index alone is insufficient
3. Read the relevant pages in full
4. Answer with citations in the form `[[page-name]]`
5. If the answer is substantive, offer to file it as a new concept or creative page

**On lint:**
1. Scan all pages in `wiki/` for: orphans, contradictions, missing pages, stale claims
2. Report a bullet list of issues with severity (fix now / fix later / flag only)
3. Fix only what the user approves

**Never:**
- Modify any file under `raw/`
- Create duplicate pages
- Leave `index.md` or `log.md` out of date after an ingest
- Add entries to `log.md` in any format other than the prefix format above

### Tooling

```bash
tools/search.sh "<query>"
```

Shells out to `qmd search` against `wiki/`. Use when `index.md` alone is not finding
relevant pages.
