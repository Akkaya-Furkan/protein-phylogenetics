# One-day scope — 24 October 2026

## Question
Does the topology of one protein family agree with known taxonomy?

## Bounded experiment
Analyze mitochondrial cytochrome b (CYTB) from 12 preselected mammal species:
four Primates (human, chimpanzee, gorilla, orangutan), four Carnivora (cat,
dog, brown bear, harbor seal), and four Artiodactyla (cow, sheep, goat, pig).
Retrieve real UniProtKB sequences and UniProt taxonomy lineages. Prefer reviewed
records, reject fragments, and select deterministically by accession if multiple
eligible records exist. Freeze the downloaded source responses and provenance.

Align full-length proteins with MAFFT, then estimate a maximum-likelihood tree
with IQ-TREE, model selection, and 1,000 ultrafast bootstrap replicates with NNI
correction. Use seed 22 and one thread. Evaluate each of the three orders as an
unrooted edge split; do not infer an ancestral root from this experiment.
Produce a labeled tree figure and machine-readable metrics.

## Numerical acceptance criteria
1. Exactly 12 distinct UniProt accessions and 12 distinct species, with 4 taxa
   per order; each input sequence is 350–400 amino acids and contains only the
   20 standard amino acids. Every species lineage contains its expected order.
2. One alignment with exactly 12 unique tips, equal aligned lengths, and
   unchanged ungapped sequences; one ML tree containing exactly the same 12 tips.
3. Exactly 1,000 bootstrap trees, each containing all 12 tips; numeric support
   in [0, 100] on all 9 internal edges of the unrooted binary ML tree.
4. Exactly 3 order-split assessments with support when a split is present.
   Agreement is a measured result, not a pass condition: disagreement must be
   reported rather than filtered away.
5. One PNG and one SVG tree figure; a JSON summary, per-order TSV, provenance
   and checksums; all pytest tests pass, including a real-tool integration test.
6. A conda environment file, download script, reproducible offline analysis
   command and GitHub Actions workflow; commit and push to main, verify the
   remote commit, and read the CI conclusion for that commit.

## Exclusions
No species-tree claim, molecular clock, ligand work, or ancestral reconstruction.
No additional families, exhaustive taxon sampling, trimming optimization, or
independent alignment/model sensitivity analysis. Single-gene mitochondrial
history can differ from species history. Ultrafast bootstrap is an approximation
to conventional bootstrap and is not a posterior probability.
