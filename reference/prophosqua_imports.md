# Imports the Package Needs but its own R Code Does Not Name

The \`CMD\_\*.R\` front ends under \`inst/application\` parse their
arguments with optparse. \`R CMD check\` does not read \`inst/\`, so
without this the dependency looks unused and the shipped applications
would break on an installation that resolved it as merely suggested.
\`prolfquapp\` declares the same import the same way, for the same
scripts.
