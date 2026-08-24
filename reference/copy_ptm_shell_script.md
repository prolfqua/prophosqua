# Copy the PTM Command Wrapper into a Working Directory

Places \`ptm.sh\` from \`inst/application/bin\` in \`workdir\`. It is
one script taking the step to run as its first argument – \`ptm.sh
dpa_dpu\`, \`ptm.sh render\`, \`ptm.sh help\` – so a work directory
holds a single wrapper however many steps the package grows.

## Usage

``` r
copy_ptm_shell_script(workdir = getwd())
```

## Arguments

- workdir:

  directory where to copy the wrapper - default is the current working
  directory.

## Value

The path copied, from \[prolfqua::script_copy_helper_vec()\].

## Details

The copy carries no analysis logic: it resolves both the list of
commands and the script each one runs from the installed package, and so
cannot drift from it. It exists so that a person at a prompt runs the
pipeline's steps through the same entry point the Snakefile uses.

## Examples

``` r
basename(copy_ptm_shell_script(tempdir()))
#> copy /home/runner/work/_temp/Library/prophosqua/application/bin/ptm.sh to /tmp/RtmpO4XXFX/ptm.sh
#> your working directory now should contain: 1 new files:
#> [1] "ptm.sh"
```
