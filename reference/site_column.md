# Name of the Site Identifier Column

A prolfquapp DEA run names the site identifier either \`site\` or
\`protein_Id_site\`, depending on its version. Every caller that joins
site annotation onto a result table has to ask which one it got.

## Usage

``` r
site_column(x)
```

## Arguments

- x:

  Data frame from a DEA result.

## Value

The name of the site column.
