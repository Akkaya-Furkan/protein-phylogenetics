# Local validation

Validation date: 24 October 2026. All observations below are from the frozen
source and executed analysis/tests, not planned outcomes.

| Acceptance item | Expected | Observed |
|---|---:|---:|
| Distinct accessions, taxon IDs, binomial species | 12 each | 12, 12, 12 |
| Taxa per sampled order | 4 each | 4, 4, 4 |
| Alignment rows / ML tips | 12 / 12 | 12 / 12 |
| Bootstrap trees with complete tip sets | 1,000 | 1000 |
| Supported internal edges | 9 | 9 |
| Prespecified order assessments | 3 | 3 |
| Tree figures (PNG + SVG) | 2 | 2 |
| Source files (response, FASTA, metadata, provenance) | 4 | 4 |
| Result files (including checksum manifest) | 19 | 19 |
| Checksummed source/result files | 22 | 22 |
| pytest cases passed | 16 | 16 |

The original alignment/inference run took 16.361 seconds
as measured by the workflow. The final local pytest invocation reported
16 passed in 21.68 seconds. Its integration test ran the complete workflow again
in a temporary directory with the real MAFFT and IQ-TREE binaries.

Tests independently verify the full downloaded response against the metadata
and FASTA, all 1,000 bootstrap tip sets, and each ML support label against the
frequency of its split in the bootstrap trees (within rounding tolerance).
All 22 SHA-256 manifest entries passed verification. Both figures were opened
and visually inspected; labels, support values and branch-length axes are visible.

Taxonomy agreement is 2/3 order splits; agreement was not a
pass condition. The absent Carnivora split is reported in the README.

Remote commit identity and CI are verified after pushing; the corresponding
GitHub Actions run is the authoritative remote test record. Local checks alone
do not establish a remote CI conclusion.
