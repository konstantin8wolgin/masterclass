# LLM Wiki — Design

**Date:** 2026-05-05
**Domain:** Lattice field theory and physics-informed kernels (research + learning)
**Vault location:** `/home/konstantin/empty/demovault/`

## Purpose

A persistent, LLM-maintained knowledge base for the user's research in lattice field theory and physics-informed kernels. Sources (papers, articles, notes) are read once and integrated into a structured wiki that compounds over time. The user curates sources and asks questions; the LLM does the summarizing, cross-referencing, filing, and bookkeeping. Code lives in a separate repository — this wiki is for understanding theory, methods, and design decisions.

## Architecture

Three layers:

1. **Raw sources** (`raw/`) — immutable user-curated documents. The LLM reads but never modifies.
2. **Wiki** (`wiki/`) — LLM-generated and LLM-maintained markdown. Two sub-layers:
   - *Source layer* (`wiki/sources/`) — one page per ingested document, stable once filed.
   - *Synthesis layer* (`wiki/concepts/`, `methods/`, `equations/`, `people/`, `creative/`, `problems/`) — pages that integrate knowledge across multiple sources and evolve with the wiki.
3. **Schema** (`CLAUDE.md`) — conventions, page formats, and workflows that govern LLM behavior across sessions.

## Directory Structure

```
demovault/
├── raw/
│   ├── papers/                 # PDFs / markdown of papers
│   ├── articles/               # web clips, blog posts
│   └── assets/                 # downloaded images
│
├── wiki/
│   ├── sources/                # one page per ingested source
│   ├── concepts/               # synthesized topic pages
│   ├── methods/                # algorithms, numerical techniques
│   ├── equations/              # key equations with derivation notes
│   ├── people/                 # recurring authors / figures
│   ├── creative/               # ideas, brainstorming, speculative directions
│   ├── problems/               # open problems and approaches to solve them
│   ├── overview.md             # evolving synthesis of how topics connect
│   ├── index.md                # catalog of all wiki pages
│   └── log.md                  # append-only operation log
│
├── tools/                      # CLI helpers (search, etc.)
│
└── CLAUDE.md                   # schema + conventions + workflows
```

## Page Formats

### Source page (`wiki/sources/<author-year-slug>.md`)

```yaml
---
title:
authors: []
year:
type: paper | article
tags: []
source_file: raw/papers/<filename>
---
```

Sections: Summary · Key claims · Methods used · Notation introduced · Open questions · Links to concept pages.

### Concept / Method / Equation page

```yaml
---
title:
tags: []
sources: []   # links to source pages
---
```

Sections: Definition · Intuition · Key results · Connections to other concepts · Open questions · Source citations.

### Problem page (`wiki/problems/<slug>.md`)

```yaml
---
title:
status: open | in-progress | resolved
---
```

Sections: Problem statement · Context · Attempts · Resolution.

### Creative page (`wiki/creative/<slug>.md`)

Free-form. Frontmatter with title and date; body follows the shape of the idea.

### `index.md`

A table grouped by directory: page link · type · one-line summary · source count · date added.

### `log.md`

Append-only. Each entry begins with a uniform prefix so it is parseable with simple shell tools:

```
## [YYYY-MM-DD] ingest | <source title> | <N pages updated>
## [YYYY-MM-DD] query  | <question>     | <pages consulted>
## [YYYY-MM-DD] lint   | <summary>
```

## Workflows

### Ingest

1. User points to a file in `raw/` (or pastes text).
2. LLM reads it and offers a 3–5 sentence summary.
3. User confirms or redirects emphasis (lightweight back-and-forth — not deep discussion every time).
4. LLM writes a source page; updates / creates relevant concept, method, equation, people pages; updates `index.md` and `log.md`.
5. LLM reports which pages were touched.

### Query

1. LLM consults `index.md` (and `qmd search` once the wiki grows large) to find relevant pages.
2. LLM reads the relevant pages and synthesizes an answer with citations.
3. If the answer is substantial, LLM offers to file it as a new page in `concepts/` or `creative/` so explorations compound back into the wiki.

### Lint

Run on request. LLM scans for:

- Orphan pages with no inbound links
- Contradictions between pages
- Concepts referenced repeatedly but lacking their own page
- Stale claims superseded by newer sources
- Missing cross-references

Reports a short punch list and fixes what the user approves.

## Tooling

- **qmd** — local hybrid search (BM25 + vector + LLM re-rank) over markdown. Installed during setup; LLM shells out to `qmd search <query>` against `wiki/` once the index file alone is no longer sufficient.
- **Obsidian plugins** (user enables manually): Dataview for frontmatter queries; optionally Marp for slide decks. Listed in `CLAUDE.md`.
- **Obsidian Web Clipper** browser extension — clips articles to `raw/articles/` as markdown. Setup notes in `CLAUDE.md`.
- **Image handling** — Obsidian configured to download attachments to `raw/assets/`; LLM views images out-of-band when needed since markdown with inline images cannot be read in one pass.

## CLAUDE.md (the schema)

`CLAUDE.md` is the contract that lets the LLM behave consistently across sessions. It will contain:

- Project description (this wiki, this domain)
- Layer boundaries (`raw/` is immutable; `wiki/` is LLM-owned)
- Page format templates (verbatim from this design)
- Workflow descriptions for ingest / query / lint
- Conventions: file naming (`<author-year-slug>` for sources, kebab-case for concepts), wikilink style, frontmatter keys
- Tooling pointers (`tools/qmd`, where to invoke search)
- A short "operational discipline" section: always update `index.md` and `log.md` on every ingest; never modify `raw/`; prefer updating existing pages over creating duplicates

## Out of scope

- Code repository — lives separately, brought in by the user in a different session if needed.
- Automated batch ingestion — workflow is user-supervised, one source at a time.
- Public sharing / publishing — vault is private and local.

## Success criteria

- The user can drop a paper into `raw/papers/` and a single LLM session produces a filed source page plus updates to all genuinely affected concept / method / equation pages.
- A query against the wiki returns an answer grounded in cited wiki pages within one round-trip.
- After ten sources have been ingested, the wiki visibly synthesizes across them — concept pages cite multiple sources, the overview reflects the evolving thesis, and `qmd search` finds material reliably.
