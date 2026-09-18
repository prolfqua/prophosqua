.PHONY: sync-quarto-assets

build build-vignettes test check-fast vignette precommit-prepare: sync-quarto-assets

sync-quarto-assets:
	Rscript data-raw/sync_quarto_assets.R

help-package:
	@echo ""
	@echo "Package-specific:"
	@echo "  make sync-quarto-assets - synchronize FGCZ Quarto assets"
