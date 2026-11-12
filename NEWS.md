# fdmcorrect 0.1.1

- Base R only: the count tables are read with `scan()` instead of `data.table::fread()`, so the package has no dependencies. Outputs are unchanged (byte-identical).

# fdmcorrect 0.1.0

- First release: the false-double-mutant error correction of nf-core/deepmutscan as an R package (`fdmcorrect`) with the command-line tool `dms-fdmcorrect`.
