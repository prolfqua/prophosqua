#!/bin/bash
#
# ptm.sh -- the entry point to prophosqua's PTM analysis steps.
#
#   ptm.sh <command> [options]
#   ptm.sh help
#
# Resolves the command's R script from the installed prophosqua, so that the
# Snakefile and a person at a prompt run the same code and neither carries a
# copy of it.
#
# One wrapper naming its command as an argument, rather than one wrapper per
# command: a work directory then holds a single file regardless of how many
# steps the package grows, and the list of steps cannot fall out of step with
# the package -- `help` reads it from what is installed.

set -euo pipefail

# Use the R that invoked us when there is one. R_HOME is set by any R process,
# so a wrapper spawned from R runs that same installation rather than whichever
# Rscript happens to be first on PATH -- which is also what "Writing R
# Extensions" 1.6 asks for, and what R CMD check enforces by putting a stub
# Rscript on PATH that refuses to run. In a plain shell R_HOME is unset and this
# is the ordinary PATH lookup.
RSCRIPT="${R_HOME:+${R_HOME}/bin/}Rscript"

APPLICATION_PATH=$("$RSCRIPT" --vanilla -e "cat(system.file('application', package = 'prophosqua'))")

if [[ -z "$APPLICATION_PATH" || ! -d "$APPLICATION_PATH" ]]; then
    echo "Error: prophosqua is not installed, or its installation is incomplete." >&2
    echo "Install it with: make -C <prophosqua checkout> install" >&2
    exit 1
fi

# The installed command scripts, as `name<TAB>summary`. The summary is the
# script's own first comment line, so it is written where the command is.
command_table() {
    local script name
    for script in "$APPLICATION_PATH"/CMD_*.R; do
        [[ -f "$script" ]] || continue
        name=${script##*/CMD_}
        name=${name%.R}
        printf '%s\t%s\n' \
            "$(printf '%s' "$name" | tr '[:upper:]' '[:lower:]')" \
            "$(sed -n '2s/^# *//p' "$script")"
    done
}

usage() {
    cat <<'EOF'
Usage: ptm.sh <command> [options]

Runs one step of the PTM analysis. Each step has its own options; ask it:

    ptm.sh <command> --help

Commands:
EOF
    command_table | awk -F'\t' '{printf "  %-16s %s\n", $1, $2}'
    printf '\nThe commands are prophosqua'"'"'s CMD_<COMMAND>.R, installed in\n%s\n' \
        "$APPLICATION_PATH"
}

if [[ $# -eq 0 ]]; then
    usage >&2
    exit 2
fi

COMMAND=$1
shift

case "$COMMAND" in
    help | --help | -h)
        usage
        exit 0
        ;;
esac

if [[ ! "$COMMAND" =~ ^[A-Za-z0-9_]+$ ]]; then
    echo "Error: '$COMMAND' is not a command name." >&2
    usage >&2
    exit 2
fi

R_SCRIPT_PATH="${APPLICATION_PATH}/CMD_$(printf '%s' "$COMMAND" | tr '[:lower:]' '[:upper:]').R"

if [[ ! -f "$R_SCRIPT_PATH" ]]; then
    echo "Error: no such command: '$COMMAND'" >&2
    usage >&2
    exit 2
fi

exec "$RSCRIPT" --vanilla "$R_SCRIPT_PATH" "$@"
