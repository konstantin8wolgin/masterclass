#!/usr/bin/env bash
# Search the wiki using qmd hybrid search (BM25 + vector + LLM rerank).
# Usage: tools/search.sh "<query>"
# Run from the demovault/ root.
set -euo pipefail

# Ensure the Node.js qmd (@tobilu/qmd) is used, not the Python qmd package.
export PATH="$HOME/.bun/bin:$HOME/.nvm/versions/node/v20.20.1/bin:$PATH"

WIKI_DIR="$(cd "$(dirname "$0")/../wiki" && pwd)"
QUERY="${1:?Usage: search.sh '<query>'}"

qmd search "$QUERY" -c wiki
