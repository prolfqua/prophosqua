#!/usr/bin/env bash
set -euo pipefail

if [ $# -ne 1 ]; then
  echo "Usage: $0 <diagram.puml>"
  exit 1
fi

PUML="$1"
BASE="${PUML%.*}"

# 1) Render SVG with PlantUML
plantuml -tsvg "$PUML"
plantuml -tpng "$PUML"

# 2) Convert SVG → EMF using Inkscape
if command -v inkscape >/dev/null 2>&1; then
  inkscape "${BASE}.svg" --export-filename="${BASE}.emf"
  echo "✅ Generated ${BASE}.emf via Inkscape"
else
  echo "❌ Inkscape not found. Install it with:"
  echo "     brew install --cask inkscape"
  exit 1
fi

echo "✅ Outputs:"
echo "   • ${BASE}.svg"
echo "   • ${BASE}.png"
