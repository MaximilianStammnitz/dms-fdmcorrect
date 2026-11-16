# fdmcorrect 0.1.2

- Libraries with no double mutants, no usable anchors or an empty count table no longer stop `evidence`, `estimate` or `correct`: they give tables without rows and are left uncorrected (or corrected with the rates of the other libraries of their pool).
- The libraries of a pool are combined in the order of their names (the `library` column of their evidence), not in the order the files were given: floating-point sums depend on the order of their terms, so the rates could differ in the last digits.
- `--method mle` needs at least 1,000 false-double reads in a pool (`eb`: 100, as before).
- The reference is read in upper case, so a soft-masked FASTA gives the same results.

# fdmcorrect 0.1.1

- Base R only: the count tables are read with `scan()` instead of `data.table::fread()`, so the package has no dependencies. Outputs are unchanged (byte-identical).

# fdmcorrect 0.1.0

- First release: the false-double-mutant error correction of nf-core/deepmutscan as an R package (`fdmcorrect`) with the command-line tool `dms-fdmcorrect`.
