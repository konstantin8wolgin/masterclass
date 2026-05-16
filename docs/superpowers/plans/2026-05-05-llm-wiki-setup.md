# LLM Wiki Setup Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Scaffold a two-layer LLM-maintained physics wiki in the existing Obsidian vault at `demovault/`, complete with CLAUDE.md schema, initial wiki pages, and qmd search tooling.

**Architecture:** Raw sources live in `raw/` (immutable); the LLM-owned wiki lives in `wiki/` with a source layer and a synthesis layer; `CLAUDE.md` is the schema governing all LLM behavior across sessions. A `tools/search.sh` wrapper shells out to `qmd` for fast full-wiki search once the index file alone is insufficient.

**Tech Stack:** Markdown, YAML frontmatter, Obsidian, qmd (local hybrid search), bash

---

### Task 1: Initialize git and directory scaffold

**Files:**
- Create: `demovault/.gitignore`
- Create: `demovault/raw/papers/.gitkeep`
- Create: `demovault/raw/articles/.gitkeep`
- Create: `demovault/raw/assets/.gitkeep`
- Create: `demovault/wiki/sources/.gitkeep`
- Create: `demovault/wiki/concepts/.gitkeep`
- Create: `demovault/wiki/methods/.gitkeep`
- Create: `demovault/wiki/equations/.gitkeep`
- Create: `demovault/wiki/people/.gitkeep`
- Create: `demovault/wiki/creative/.gitkeep`
- Create: `demovault/wiki/problems/.gitkeep`
- Create: `demovault/tools/.gitkeep`

- [ ] **Step 1: Initialize git repo**

```bash
cd /home/konstantin/empty/demovault
git init
```

Expected: `Initialized empty Git repository in .../demovault/.git/`

- [ ] **Step 2: Create .gitignore**

Create `demovault/.gitignore`:

```
.obsidian/workspace.json
.obsidian/workspace-mobile.json
.trash/
```

- [ ] **Step 3: Create all directories with placeholder files**

```bash
mkdir -p raw/papers raw/articles raw/assets
mkdir -p wiki/sources wiki/concepts wiki/methods wiki/equations wiki/people wiki/creative wiki/problems
mkdir -p tools
touch raw/papers/.gitkeep raw/articles/.gitkeep raw/assets/.gitkeep
touch wiki/sources/.gitkeep wiki/concepts/.gitkeep wiki/methods/.gitkeep
touch wiki/equations/.gitkeep wiki/people/.gitkeep wiki/creative/.gitkeep wiki/problems/.gitkeep
touch tools/.gitkeep
```

- [ ] **Step 4: Verify structure**

```bash
find . -not -path './.git/*' -not -path './.obsidian/*' | sort
```

Expected output includes:
```
./raw/papers/.gitkeep
./raw/articles/.gitkeep
./raw/assets/.gitkeep
./wiki/sources/.gitkeep
./wiki/concepts/.gitkeep
./wiki/methods/.gitkeep
./wiki/equations/.gitkeep
./wiki/people/.gitkeep
./wiki/creative/.gitkeep
./wiki/problems/.gitkeep
./tools/.gitkeep
./.gitignore
```

- [ ] **Step 5: Commit**

```bash
git add .gitignore raw/ wiki/ tools/
git commit -m "feat: scaffold llm wiki directory structure"
```

---

### Task 2: Create initial wiki pages

**Files:**
- Create: `demovault/wiki/index.md`
- Create: `demovault/wiki/log.md`
- Create: `demovault/wiki/overview.md`

- [ ] **Step 1: Create index.md**

Create `demovault/wiki/index.md`:

```markdown
# Wiki Index

A catalog of all pages in this wiki. Updated by the LLM on every ingest.

## Sources

| Page | Summary | Date Added |
|------|---------|------------|
| *(none yet)* | | |

## Concepts

| Page | Summary | Sources | Date Added |
|------|---------|---------|------------|
| *(none yet)* | | | |

## Methods

| Page | Summary | Sources | Date Added |
|------|---------|---------|------------|
| *(none yet)* | | | |

## Equations

| Page | Summary | Sources | Date Added |
|------|---------|---------|------------|
| *(none yet)* | | | |

## People

| Page | Summary | Date Added |
|------|---------|------------|
| *(none yet)* | | |

## Creative

| Page | Summary | Date Added |
|------|---------|------------|
| *(none yet)* | | |

## Problems

| Page | Status | Date Added |
|------|--------|------------|
| *(none yet)* | | |
```

- [ ] **Step 2: Create log.md**

Create `demovault/wiki/log.md`:

```markdown
# Operation Log

Append-only record of all wiki operations. Each entry begins with a parseable prefix.

Format:
- `## [YYYY-MM-DD] ingest | <source title> | <N pages updated>`
- `## [YYYY-MM-DD] query  | <question summary> | <pages consulted>`
- `## [YYYY-MM-DD] lint   | <summary of issues found/fixed>`

---

## [2026-05-05] init | Wiki created | 3 pages created
```

- [ ] **Step 3: Create overview.md**

Create `demovault/wiki/overview.md`:

```markdown
# Overview

*This page is the evolving synthesis of how topics in this wiki connect. The LLM updates it after significant ingests or when a new cross-cutting theme emerges.*

## Domain

This wiki covers **lattice field theory** and **physics-informed kernels** — two areas that intersect around discretization of continuous field theories and machine-learning methods constrained by physical laws.

## Current understanding

*(To be built up as sources are ingested.)*

## Open threads

*(Questions and themes that cut across multiple sources, to be resolved over time.)*

## Key tensions

*(Places where sources disagree or where the theory is unsettled.)*
```

- [ ] **Step 4: Verify files**

```bash
ls -1 wiki/
```

Expected:
```
concepts/
creative/
equations/
index.md
log.md
methods/
overview.md
people/
problems/
sources/
```

- [ ] **Step 5: Commit**

```bash
git add wiki/index.md wiki/log.md wiki/overview.md
git commit -m "feat: add initial wiki pages (index, log, overview)"
```

---

### Task 3: Write CLAUDE.md schema

**Files:**
- Create: `demovault/CLAUDE.md`

- [ ] **Step 1: Write CLAUDE.md**

Create `demovault/CLAUDE.md`:

```markdown
# LLM Wiki — Schema

This file governs how the LLM behaves in this vault across all sessions. Read it at the start of every session before doing anything else.

## What this wiki is

A persistent, LLM-maintained knowledge base for research in **lattice field theory** and **physics-informed kernels**. The human curates sources and asks questions. The LLM does all the filing, cross-referencing, and maintenance.

Code lives in a separate repository. This wiki is for understanding theory, methods, and design decisions — not for storing code.

## Layer boundaries

| Directory | Owner | Rule |
|-----------|-------|------|
| `raw/` | Human | LLM reads, never modifies |
| `wiki/` | LLM | LLM creates and updates everything |
| `tools/` | Shared | Human installs tools; LLM invokes them |
| `CLAUDE.md` | Shared | Evolves by mutual agreement |

## File naming conventions

- Source pages: `wiki/sources/<firstauthor-year-slug>.md` (e.g., `wilson-1974-confinement.md`)
- All other wiki pages: `wiki/<dir>/<kebab-case-title>.md`
- Raw papers: keep original filename; note it in the source page frontmatter

## Page format templates

### Source page

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

Sections (in order):
1. **Summary** — 3–5 sentences covering the paper's contribution
2. **Key claims** — bullet list of the main claims or results
3. **Methods used** — what techniques, discretizations, algorithms appear
4. **Notation introduced** — any non-standard notation to remember
5. **Open questions** — things this paper leaves unresolved
6. **Links** — wikilinks to relevant concept/method/equation pages

### Concept page

```yaml
---
title: ""
tags: []
sources: []  # list of wikilinks to source pages
---
```

Sections:
1. **Definition** — precise statement of what this concept is
2. **Intuition** — physical or geometric intuition
3. **Key results** — important theorems, limits, known properties
4. **Connections** — wikilinks to related concepts, methods, equations
5. **Open questions** — what is still unsettled
6. **Sources** — inline citations to source pages

### Method page

Same frontmatter and section structure as concept page, with an additional section:
- **Algorithm** — step-by-step pseudocode or prose description of the procedure

### Equation page

```yaml
---
title: ""
tags: []
sources: []
---
```

Sections:
1. **Statement** — the equation in LaTeX
2. **Variables** — definition of each symbol
3. **Derivation sketch** — key steps or reference to where derivation lives
4. **Limits and special cases**
5. **Sources**

### Problem page

```yaml
---
title: ""
status: open  # open | in-progress | resolved
date_opened: YYYY-MM-DD
date_resolved:  # leave blank if not resolved
---
```

Sections:
1. **Problem statement**
2. **Context** — why this matters, where it arose
3. **Attempts** — what has been tried, with results
4. **Resolution** — fill in when solved

### Creative page

```yaml
---
title: ""
date: YYYY-MM-DD
---
```

Free-form body. Use whatever structure fits the idea.

## Style rules

These apply to every response, every wiki page, every log entry — no exceptions.

1. **No emojis.** Never use emojis in wiki pages, log entries, responses, or any other output.
2. **Clear order.** Always present information in a logical, numbered or bulleted sequence. Avoid unordered prose dumps. Steps are numbered, lists are bulleted, tables have headers.
3. **Concise.** Prefer shorter over longer. One sentence where one sentence is enough. No padding, no filler, no restating what the user can read in the page.
4. **No commentary on your own actions.** Don't explain that you are "now updating index.md" — just do it and report the result at the end.

## Operational discipline

**On every ingest:**
1. Write or update the source page
2. Update or create any concept / method / equation / people pages that the source touches
3. Update `wiki/index.md` — add or update rows for every page touched
4. Append an entry to `wiki/log.md` with the format: `## [YYYY-MM-DD] ingest | <title> | <N pages updated>`
5. Report to the user which pages were created or updated

**On every query:**
1. Read `wiki/index.md` first to find relevant pages
2. Use `tools/search.sh <query>` when the index alone is insufficient (wiki has grown large)
3. Read the relevant pages in full
4. Answer with citations in the form `[[page-name]]`
5. If the answer is substantive, offer to file it as a new concept or creative page

**On lint:**
1. Scan all pages in `wiki/` for: orphans (no inbound wikilinks), contradictions, concepts mentioned without their own page, stale claims
2. Report a bullet list of issues with severity (fix now / fix later / flag only)
3. Fix only what the user approves

**Never:**
- Modify any file under `raw/`
- Create duplicate pages (search index.md and use `tools/search.sh` before creating)
- Leave `index.md` or `log.md` out of date after an ingest
- Add entries to `log.md` in any format other than the prefix format above

## Tooling

### qmd search

```bash
tools/search.sh "<query>"
```

This shells out to `qmd search` against the `wiki/` directory. Use it when `index.md` alone is not finding the relevant pages (wiki has more than ~50 pages, or the query is subtle).

### Obsidian plugins (enable manually in Obsidian settings)

- **Dataview** — allows querying frontmatter. Useful for e.g. listing all open problems.
- **Marp** (optional) — render markdown as slide decks. Useful for generating presentations from wiki content.

### Obsidian Web Clipper

Browser extension that converts web articles to markdown. After installing:
1. Configure output directory to `raw/articles/`
2. Clip article → it appears in `raw/articles/` → run ingest workflow

### Image attachments

In Obsidian Settings → Files and links:
- Set "Attachment folder path" to `raw/assets`
- Bind "Download attachments for current file" to a hotkey (e.g. Ctrl+Shift+D)

Note: the LLM cannot read markdown with inline images in one pass. Workflow: LLM reads the markdown text first, then views referenced images separately if needed.

## Ingest workflow (detailed)

1. User says: "ingest raw/papers/filename.pdf" (or pastes text, or gives a URL)
2. LLM reads the source
3. LLM gives a 3–5 sentence summary and asks for confirmation or redirect
4. User confirms or says what to emphasize differently
5. LLM writes/updates source page at `wiki/sources/<slug>.md`
6. LLM identifies which concepts, methods, equations, people are relevant and updates those pages (or creates them if they don't exist — check index.md first)
7. LLM updates `wiki/index.md`
8. LLM appends to `wiki/log.md`
9. LLM reports: "Filed as [[source-slug]]. Updated: [[concept-a]], [[method-b]]. Created: [[equation-c]]."
```

- [ ] **Step 2: Verify CLAUDE.md exists and is non-empty**

```bash
wc -l CLAUDE.md
```

Expected: more than 100 lines.

- [ ] **Step 3: Commit**

```bash
git add CLAUDE.md
git commit -m "feat: add CLAUDE.md schema with page templates and workflows"
```

---

### Task 4: Install qmd and create search wrapper

**Files:**
- Create: `demovault/tools/search.sh`
- Modify: `demovault/tools/.gitkeep` (delete it once search.sh exists)

- [ ] **Step 1: Install qmd**

Check if Go is available:

```bash
go version
```

If Go is available, install qmd:

```bash
go install github.com/tobi/qmd@latest
```

If Go is not available, download the prebuilt binary from https://github.com/tobi/qmd/releases and place it somewhere on `$PATH` (e.g., `~/.local/bin/qmd`).

Verify:

```bash
qmd --help
```

Expected: qmd help output with `search`, `index`, and other subcommands.

- [ ] **Step 2: Create search wrapper**

Create `demovault/tools/search.sh`:

```bash
#!/usr/bin/env bash
# Search the wiki using qmd hybrid search (BM25 + vector + LLM rerank).
# Usage: tools/search.sh "<query>"
# Run from the demovault/ root.
set -euo pipefail

WIKI_DIR="$(dirname "$0")/../wiki"
QUERY="${1:?Usage: search.sh '<query>'}"

qmd search "$QUERY" --dir "$WIKI_DIR"
```

```bash
chmod +x tools/search.sh
```

- [ ] **Step 3: Index the wiki**

```bash
cd /home/konstantin/empty/demovault
qmd index wiki/
```

Expected: qmd indexes the wiki directory (may be fast since the wiki is currently empty).

- [ ] **Step 4: Run a test search**

```bash
tools/search.sh "lattice field theory"
```

Expected: either no results (wiki is empty) or a list of matching pages — no error. If qmd outputs "no results found" that is acceptable; a crash or "command not found" is not.

- [ ] **Step 5: Remove placeholder and commit**

```bash
rm tools/.gitkeep
git add tools/search.sh
git rm tools/.gitkeep
git commit -m "feat: add qmd search wrapper at tools/search.sh"
```

---

### Task 5: Smoke-test the full setup

This task has no new files — it verifies the whole system is coherent.

- [ ] **Step 1: Check final directory structure**

```bash
cd /home/konstantin/empty/demovault
find . -not -path './.git/*' -not -path './.obsidian/*' | sort
```

Expected structure:
```
.
./.gitignore
./CLAUDE.md
./raw/
./raw/articles/.gitkeep
./raw/assets/.gitkeep
./raw/papers/.gitkeep
./tools/
./tools/search.sh
./wiki/
./wiki/concepts/.gitkeep
./wiki/creative/.gitkeep
./wiki/equations/.gitkeep
./wiki/index.md
./wiki/log.md
./wiki/methods/.gitkeep
./wiki/overview.md
./wiki/people/.gitkeep
./wiki/problems/.gitkeep
./wiki/sources/.gitkeep
```

- [ ] **Step 2: Verify git log shows all commits**

```bash
git log --oneline
```

Expected (4 commits in order, most recent first):
```
<hash> feat: add qmd search wrapper at tools/search.sh
<hash> feat: add CLAUDE.md schema with page templates and workflows
<hash> feat: add initial wiki pages (index, log, overview)
<hash> feat: scaffold llm wiki directory structure
```

- [ ] **Step 3: Confirm Obsidian can open the vault**

Open Obsidian and open `demovault/` as a vault. Navigate to `wiki/index.md`. Verify it renders without errors. Enable Dataview plugin if desired (Settings → Community plugins → Browse → Dataview).

- [ ] **Step 4: Final commit if any loose changes**

```bash
git status
# If clean: nothing to do.
# If there are changes: review and commit with an appropriate message.
```

---

## Appendix: First ingest

Once setup is complete, test the whole pipeline with your first paper:

1. Copy or save the paper to `raw/papers/`
2. In a new Claude Code session (so CLAUDE.md is loaded fresh), say:
   `"Ingest raw/papers/<filename>"`
3. Confirm the summary
4. Check that `wiki/sources/`, `wiki/index.md`, and `wiki/log.md` were all updated
5. Run `tools/search.sh "<concept from the paper>"` and verify it returns the new page
