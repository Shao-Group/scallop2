# GNN transcript-to-bundle module

Updated: 2026-10-01

## Goal

On the `gnn` branch, load annotation transcripts from `-b <annotation.gtf>` and alignments from `-i <alignments.bam>`, assign transcripts from direct read-bundle evidence, and guarantee that every assigned transcript is represented by a source-to-sink path in the resulting splice graph.

## Current state

- Active branch: `gnn`; transcript-to-bundle assignment was committed as `2725ed0`.
- The `main.cc` assembly call and `assembler::assemble()` BAM loop are active.
- There are unrelated untracked files above `src/`; leave them alone.
- Commit `800f43d` added the GNN processing skeleton and commit `c74ea76` added `-b` plus `genome gm` loading.
- No `AGENTS.md` was found from `src/` through the checked parent directories.

## Relevant flow

1. `main.cc` parses options, previews the BAM, constructs `assembler`, and normally calls `assembler::assemble()`.
2. `assembler::assemble()` reads BAM records, creates `hit` objects, splits them between plus/minus `bundle_base` instances, and flushes a pool of bundles.
3. `assembler::process_gnn()` constructs a `transcript_matcher` directly from each `bundle_base`, queries annotation candidates, evaluates them, and updates per-transcript counts.
4. The accepted transcript models are then passed into `bundle`. Their exons, junctions, and outer boundaries become structural graph evidence, and `bundle::build(1, true)` builds and revises the splice graph while preserving every assigned source-to-sink transcript path.
5. `genome` loads genes/transcripts from GTF. Exons are sorted and adjacent exon records are merged by `gene::shrink()` / `transcript::shrink()`.

All relevant genomic intervals are 0-based and right-open. A GTF exon `start..end` becomes `[start-1,end)`. Bundle bounds, partial exons, graph vertices, and transcript exons can therefore be compared directly.

## Proposed definitions

For transcript `T` and read bundle `B`:

- `same_locus`: `T.seqname == B.chrm` and `T.strand == B.strand`; for an ambiguous bundle strand, the matcher applies the same majority-of-`XS` inference previously used by `bundle::compute_strand()`.
- `exon_ratio`: length of `union(T.exons) intersect union(B.mmap)` divided by `T.length()`.
- `splicing_position_ratio`: number of unique transcript donor/acceptor coordinates found among the bundle's observed donor/acceptor coordinates divided by the number of unique transcript donor/acceptor coordinates.
- Assign when `same_locus`, `exon_ratio >= 0.5`, and `splicing_position_ratio >= 0.5` all hold by default.

At the default 0.5 setting, a four-exon transcript has six splice positions and needs at least three of them in the bundle. A single-exon transcript has no splice positions, so the splice-position rule is vacuously true and exon overlap decides it. Boundary-only contact contributes zero exon bases. Assignment is many-to-many.

The exact-strand rule means a `.` bundle matches only a `.` transcript. This follows the current wording but may yield no labels for ambiguous bundles in unstranded libraries; confirm or revise this policy when real data exposes the case.

## Implemented design

The source-level `transcript_index` component replaces use of `genome::locate_gene()` for candidate retrieval:

- Flatten all nonempty loaded transcripts into stable records containing a numeric ID, gene/transcript identity, chromosome, strand, bounds, exons, and precomputed introns/exonic length.
- Partition records by `(chromosome, strand)`.
- For each partition, build a Boost split interval map whose codomain is the set of transcript IDs active over a genomic segment. This reuses the repository's `ROI` right-open interval type and does not scan unrelated chromosomes or strands.
- Query with the bundle half-open span `[lpos,rpos)`. Return numeric record IDs; do exact scoring separately.

This index is built once after the GTF is loaded and remains immutable while BAM bundles stream through `assembler`.

Files:

- `transcript_index.h/.cc`: immutable transcript records and chromosome/strand interval maps.
- `transcript_match.h/.cc`: `bundle_base` exon/splice-position extraction and exact criterion scoring.
- `assembler.cc`: one-time index construction and per-bundle candidate assignment.
- `config.h/.cc`: `-b` validation and configurable exon/splice-position thresholds.
- `tests/test_transcript_assignment.cc`: focused index, criterion, and `bundle_base::mmap` extraction tests.
- `tests/data/transcript_assignment.{sam,gtf}`: end-to-end fixture.

## Completed implementation

1. Repaired the branch baseline.
   - Restored `extern string ref_file;` and removed the duplicate `ref_file1` declaration.
   - Removed the temporary early return and restored `asmb.assemble()`.
   - Added `-i`/`-b` value and file-open checks, documented `-b`, and require it outside preview mode.

2. Added `transcript_index.h/.cc` and registered them in `src/Makefile.am`.
   - Builds immutable flattened transcript records from `genome`.
   - Builds one split interval map of transcript-ID sets per chromosome/strand.
   - Exposes `query(chromosome, strand, left, right)` plus a const record accessor.

3. Added an independently testable matcher.
   - Joins the covered intervals in `bundle_base::mmap` into an exon union.
   - Extracts unique donor and acceptor positions from the junction pairs in bundle hits.
   - Computes exon overlap without double counting and counts transcript splice-position membership.
   - Returns both the boolean decision and raw counts for debugging and later GNN features.

4. Integrated it at the GNN processing boundary.
   - Builds the index once in `assembler` after loading the GTF.
   - Keeps streaming/batching behavior unchanged, but routes every pool flush through `process_gnn()`, including the final pool.
   - Queries candidates using bundle chromosome, inferred strand, and bounds, then runs the matcher directly against `bundle_base`.
   - Clears the processed pool in `process_gnn()` to prevent repeated batches.

5. Added deterministic diagnostic output.
   - Verbosity 1 prints candidate and assignment counts per bundle.
   - Verbosity 2 prints assigned `gene_id`/`transcript_id` and raw exon/splice-position match counts.
   - Assignment order is stable by transcript numeric ID.

6. Verified with focused synthetic fixtures and an end-to-end run.
   - Index: chromosomes/strands, boundary-touch cases, and multiple returned candidates.
   - Matcher: exact 50% exon coverage, below-threshold coverage, one-of-three versus two-of-three junctions, graph splice-edge extraction, and single-exon transcripts.
   - Pipeline: a tiny coordinate-sorted BAM plus GTF proving batch and final flush behavior and deterministic assignments.
   - Force or cleanly rebuild all affected sources; do not accept a no-op `make` against stale objects as validation.

7. Fixed a historical mixed-build crash in the former `bundle::assigned_transcripts` field.
   - The existing dependency file for `bundle.cc` was a dummy, so changing `bundle.h` rebuilt `assembler.o` but left an older `bundle.o`; that constructor did not construct the newly added vector.
   - The field was explicitly initialized at the time; the direct-`bundle_base` refactor later removed it entirely.

8. Added per-transcript bundle assignment counts.
   - `assembler` maintains one 64-bit count per stable transcript-index ID and increments it once for each bundle that accepts that transcript.
   - At successful completion, every annotation transcript, including zero-count transcripts, is printed and saved as deterministic TSV.
   - The default output is `<gtf-file>.bundle_counts.tsv`; `--transcript_bundle_count_file <filename>` overrides it.
   - TSV columns are `transcript_index`, `transcript_id`, `gene_id`, `chromosome`, `strand`, and `bundle_count`.

9. Refactored matching to operate directly on `bundle_base`.
   - Removed splice-graph construction from `process_gnn()` and removed the unused `bundle::assigned_transcripts` field.
   - Builds the joined exon interval map directly from `bundle_base::mmap`, following `region::build_join_interval_map()` semantics.
   - Replaced exact `int64_t` junction storage with a set of `int32_t` donor/acceptor positions collected from bundle hits.
   - Scores the fraction of unique query splice positions present in the bundle; the default threshold is 0.5.
   - Renamed the option to `--min_bundle_transcript_splicing_position_overlap`; the former junction-named option remains accepted as a compatibility alias.

10. Added assigned transcripts to splice-graph construction.
   - `bundle` stores copies of its assigned annotation transcripts and exposes them to `bundle_bridge` as reference transcripts.
   - Annotation-only junctions are inserted into `bundle::junctions` with minimum structural support; junctions already supported by reads are not duplicated.
   - Transcript starts and ends are included among region partition coordinates, and bundle bounds expand when an assigned transcript extends beyond the read-supported span.
   - The union of assigned exon intervals is compared with the structural interval map; only uncovered subintervals are added, so existing read-coverage values are unchanged.
   - The newly added interval union is retained on `bundle`, and every partial exon and corresponding graph vertex records the number of newly added bases it contains.
   - After graph revision, missing source, contiguous-exon, junction, and sink edges are restored at minimum weight, vertices on assigned paths are retained, and every assigned transcript is explicitly validated as a source-to-sink path.

## Remaining decisions

- The phrase “share half of the exon regions” could mean half of exon count rather than half of exonic bases. The proposed definition uses exonic bases because it handles partial overlaps and unequal exon lengths predictably.
- Exact same-strand matching for `.` bundles may be too strict for unstranded libraries. Keep it strict initially as requested and expose counts of skipped ambiguous bundles.
- Matching describes direct read-bundle evidence, while assigned transcripts are subsequently guaranteed as paths in the revised graph. The per-vertex `newly_added_length` field now explicitly identifies annotation-supplied exon bases for downstream GNN serialization.
- The final GNN tensor/graph serialization is not yet specified. Matching remains independent of serialization; bundle indices and deterministic transcript counts are available for the next layer.

## Validation performed

- Focused tests pass for chromosome/strand partitions, boundary-touch non-overlap, multiple candidates, exact 50% exon coverage, below-threshold coverage, splice-position threshold behavior, single-exon handling, strand rejection, and exon extraction from a `bundle_base` coverage map.
- The assigned-transcript graph test covers expanded transcript bounds, annotation-only junction insertion, boundary-based region partitioning, adding only the union gaps absent from `fmap` without changing existing coverage, retained missing intervals, per-partial-exon/per-vertex newly-added lengths, post-revision path restoration, and final source-to-sink path validation.
- Every `src/*.cc` translation unit compiles and the full executable links successfully in a clean temporary directory.
- End-to-end fixture output: two annotations loaded; the plus-strand bundle fetched one candidate and assigned `t1` with exon overlap `100/100` and junction overlap `1/1`; the minus-strand annotation was excluded.
- `git diff --check` passes.
- Reproduced the reported large-input crash under GDB at `assembler.cc` while pushing into an unconstructed `assigned_transcripts` vector; object timestamps confirmed the mixed build.
- After rebuilding `bundle.o`, the same `tests/star.sort.bam` and `tests/scallop2.gtf` run passed the former crash point and ran for 180 seconds without a segmentation fault before an intentional timeout.
- The two-transcript end-to-end fixture prints and saves `t1` with one assigned bundle and the opposite-strand `t2` with zero assigned bundles.
- Automake regeneration remains unavailable in this restricted workspace because it tries to write the parent `autom4te.cache`; canonical source registration is complete in tracked `src/Makefile.am`.

## Next session

Start with `git status --short --branch` and read this file plus `skills/scallop2-gnn/SKILL.md`. Direct `bundle_base` assignment and assigned-transcript graph preservation are implemented and validated. The next functional task is to define the GNN example serialization and how it should encode assigned transcript paths and annotation-derived structural edges.
