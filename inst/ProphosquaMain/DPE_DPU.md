# Differential PTM-feature Expression (DPE)
## What it tests
- tests raw PTM-feature change between conditions, without any correction for its parent protein.

## When to use
To flag PTM-features whose signal changes, even if the protein itself is also up- or down-regulated.

## How it works
- Log_2-transform PTM-feature intensities and apply normalization (prolfqua::LFQTransformer)
- Fit your usual differential-expression model (prolfqua::build_model)

# Differential PTM-feature Usage (DPU)

## What it tests
- tests **protein normalized** changes of PTM-features (see Figure below)

## When to use
To home in whether the **fraction** of protein thats modified shifts.

## How it works

- Summarize protein intensities per sample
- For each PTM compute log2(PTM) - log2(protein) and fit your usual differential-expression model
- or fit your usual DE model and compute log2FC(PTM) - log2FC(Protein) and recalculate p-values.
 