# CLAUDE.md

Project-specific instructions for Claude Code.

## Rules

- Do NOT modify DESCRIPTION or NAMESPACE unless explicitly asked
- One function, one task! Keep functions minimal and focused.

## Vignettes

- ONLY edit vignettes in `vignettes/`, NEVER in `doc/`
- `doc/` contains built output — it is regenerated from `vignettes/` via
  `make build-vignettes` or
  [`devtools::build_vignettes()`](https://devtools.r-lib.org/reference/build_vignettes.html)
- Workflow: edit `vignettes/*.Rmd` -\> build vignettes -\> `doc/` is
  updated automatically

## Development

Use the Makefile for all package development tasks:

- `make test` - run testthat tests
- `make check-fast` - R CMD check without vignettes (quick validation)
- `make check` - full R CMD check
- `make document` - generate roxygen2 docs
- `make build` - build tarball
- `make install` - install package locally
- `make lint` - run lintr
- `make format` - format with air
- `make clean` - remove build artifacts
