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
2. Update or create any concept / method / equation  pages that the source touches
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
6. LLM identifies which concepts, methods, equations are relevant and updates those pages (or creates them if they don't exist — check index.md first)
7. LLM updates `wiki/index.md`
8. LLM appends to `wiki/log.md`
9. LLM reports: "Filed as [[source-slug]]. Updated: [[concept-a]], [[method-b]]. Created: [[equation-c]]."
