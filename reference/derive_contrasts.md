# Derive Contrasts from a Sample Annotation

Two annotation layouts are in use and both have to work:

## Usage

``` r
derive_contrasts(annot, annot_label = "the annotation")
```

## Arguments

- annot:

  Data frame read from the annotation file.

- annot_label:

  Name used in error messages, normally the file name.

## Value

A named character vector of contrast expressions.

## Details

\* \*\*explicit\*\* — \`Contrast\` holds the linear-model expression and
\`ContrastName\` its label. Taken as given. \* \*\*dataset\*\* —
\`Group\` names the experimental group and \`Control\` marks it as
control (\`C\`) or treatment (\`T\`). Every treatment-versus-control
pair is generated, named \`\<treatment\>\_vs\_\<control\>\` and
expressed with the \`G\_\` prefix prolfqua gives group coefficients.

Column names of the dataset layout are matched case-insensitively:
annotations in the wild spell the second column both \`Control\` and
\`CONTROL\`.
