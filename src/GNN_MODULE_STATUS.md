# GNN transcript-to-bundle module

Updated: 2026-09-30

## Goal

On the `gnn` branch, load annotation transcripts from `-b <annotation.gtf>` and alignments from `-i <alignments.bam>`, build one splice graph for each accepted read bundle, and attach every annotation transcript that satisfies the initial assignment criterion.

## Current state

- Active branch: `gnn` at `c74ea76` (`load gtf file`), tracking `origin/gnn`.
- The `main.cc` assembly call and `assembler::assemble()` BAM loop are active.
- There are unrelated untracked files above `src/`; leave them alone.
- Commit `800f43d` added the GNN processing skeleton and commit `c74ea76` added `-b` plus `genome gm` loading.
- No `AGENTS.md` was found from `src/` through the checked parent directories.

## Relevant flow

1. `main.cc` parses options, previews the BAM, constructs `assembler`, and normally calls `assembler::assemble()`.
2. `assembler::assemble()` reads BAM records, creates `hit` objects, splits them between plus/minus `bundle_base` instances, and flushes a pool of bundles.
3. `assembler::process_gnn()` constructs `bundle bd(bb)`, calls `bd.build(1, true)`, queries annotation candidates, assigns matching transcript IDs, and prints diagnostics.
4. The `bundle` constructor prepares aligned intervals, observed junctions, regions, partial exons, and mappings. `build(1, true)` creates and revises `bd.gr` and builds `bd.hs`.
5. `genome` loads genes/transcripts from GTF. Exons are sorted and adjacent exon records are merged by `gene::shrink()` / `transcript::shrink()`.

All relevant genomic intervals are 0-based and right-open. A GTF exon `start..end` becomes `[start-1,end)`. Bundle bounds, partial exons, graph vertices, and transcript exons can therefore be compared directly.

## Proposed definitions

For transcript `T` and final revised graph `G`:

- `same_locus`: `T.seqname == G.chrm` and `T.strand == G.strand`.
- `exon_ratio`: length of `union(T.exons) intersect union(G interior vertex intervals)` divided by `T.length()`.
- `junction_ratio`: number of exact transcript introns `(exon[i].second, exon[i+1].first)` present among graph splice edges divided by `T.exons.size() - 1`.
- Assign when `same_locus`, `exon_ratio >= 0.5`, and `junction_ratio >= 0.5` all hold.

At the default 0.5 setting, the comparisons are equivalent to `2 * numerator >= denominator`. Thus one of three junctions is insufficient and two of three pass. A zero-junction transcript passes the junction rule and is decided by exon overlap. Boundary-only contact contributes zero exon bases. Assignment is many-to-many.

The exact-strand rule means a `.` bundle matches only a `.` transcript. This follows the current wording but may yield no labels for ambiguous bundles in unstranded libraries; confirm or revise this policy when real data exposes the case.

## Implemented design

The source-level `transcript_index` component replaces use of `genome::locate_gene()` for candidate retrieval:

- Flatten all nonempty loaded transcripts into stable records containing a numeric ID, gene/transcript identity, chromosome, strand, bounds, exons, and precomputed introns/exonic length.
- Partition records by `(chromosome, strand)`.
- For each partition, build a Boost split interval map whose codomain is the set of transcript IDs active over a genomic segment. This reuses the repository's `ROI` right-open interval type and does not scan unrelated chromosomes or strands.
- Query with the graph/bundle half-open span `[lpos,rpos)`. Return numeric record IDs; do exact scoring separately.

This index is built once after the GTF is loaded and remains immutable while BAM bundles stream through `assembler`.

Files:

- `transcript_index.h/.cc`: immutable transcript records and chromosome/strand interval maps.
- `transcript_match.h/.cc`: final-graph exon/junction extraction and exact criterion scoring.
- `assembler.cc`: one-time index construction and per-bundle candidate assignment.
- `bundle.h`: public `assigned_transcripts` numeric IDs.
- `config.h/.cc`: `-b` validation and configurable exon/junction thresholds.
- `tests/test_transcript_assignment.cc`: focused index, criterion, and graph extraction tests.
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
   - Extracts the union of interior graph vertex intervals, excluding source and sink.
   - Extracts exact splice-junction pairs from non-contiguous interior graph edges.
   - Computes exon overlap without double counting and counts exact transcript-junction membership.
   - Returns both the boolean decision and raw counts for debugging and later GNN features.

4. Integrated it at the GNN processing boundary.
   - Builds the index once in `assembler` after loading the GTF.
   - Keeps streaming/batching behavior unchanged, but routes every pool flush through `process_gnn()`, including the final pool.
   - After `bd.build(1, true)`, queries candidates using `bd.gr.chrm`, `bd.gr.strand`, and graph bounds, runs the matcher, and stores accepted IDs in `bd.assigned_transcripts`.
   - Clears the processed pool in `process_gnn()` to prevent repeated batches.

5. Added deterministic diagnostic output.
   - Verbosity 1 prints candidate and assignment counts per bundle.
   - Verbosity 2 prints assigned `gene_id`/`transcript_id` and raw exon/junction match counts.
   - Assignment order is stable by transcript numeric ID.

6. Verified with focused synthetic fixtures and an end-to-end run.
   - Index: chromosomes/strands, boundary-touch cases, and multiple returned candidates.
   - Matcher: exact 50% exon coverage, below-threshold coverage, one-of-three versus two-of-three junctions, graph splice-edge extraction, and single-exon transcripts.
   - Pipeline: a tiny coordinate-sorted BAM plus GTF proving batch and final flush behavior and deterministic assignments.
   - Force or cleanly rebuild all affected sources; do not accept a no-op `make` against stale objects as validation.

7. Fixed a mixed-build crash in `bundle::assigned_transcripts`.
   - The existing dependency file for `bundle.cc` was a dummy, so changing `bundle.h` rebuilt `assembler.o` but left an older `bundle.o`; that constructor did not construct the newly added vector.
   - Explicitly initialized `assigned_transcripts` in `bundle::bundle()`. This both documents the required construction and forces `bundle.cc` to rebuild when the fix is applied.

## Remaining decisions

- The phrase “share half of the exon regions” could mean half of exon count rather than half of exonic bases. The proposed definition uses exonic bases because it handles partial overlaps and unequal exon lengths predictably.
- Exact same-strand matching for `.` bundles may be too strict for unstranded libraries. Keep it strict initially as requested and expose counts of skipped ambiguous bundles.
- Matching against the revised graph makes labels consistent with the actual GNN input, but a true annotation junction removed during revision will count as absent. This is intentional in the proposed design and should be measured.
- The final GNN tensor/graph serialization is not yet specified. Matching remains independent of serialization; bundle IDs and diagnostic counts are ready for that next layer.

## Validation performed

- Focused tests pass for chromosome/strand partitions, boundary-touch non-overlap, multiple candidates, exact 50% exon coverage, below-threshold coverage, one-of-three versus two-of-three junctions, single-exon handling, strand rejection, and extraction from a real `splice_graph` object.
- Every `src/*.cc` translation unit compiles and the full executable links successfully in a clean temporary directory.
- End-to-end fixture output: two annotations loaded; the plus-strand bundle fetched one candidate and assigned `t1` with exon overlap `100/100` and junction overlap `1/1`; the minus-strand annotation was excluded.
- `git diff --check` passes.
- Reproduced the reported large-input crash under GDB at `assembler.cc` while pushing into an unconstructed `assigned_transcripts` vector; object timestamps confirmed the mixed build.
- After rebuilding `bundle.o`, the same `tests/star.sort.bam` and `tests/scallop2.gtf` run passed the former crash point and ran for 180 seconds without a segmentation fault before an intentional timeout.
- Automake regeneration remains unavailable in this restricted workspace because it tries to write the parent `autom4te.cache`; canonical source registration is complete in tracked `src/Makefile.am`.

## Next session

Start with `git status --short --branch` and read this file plus `skills/scallop2-gnn/SKILL.md`. The transcript assignment module is implemented and validated. The next functional task is to define and add the GNN example serialization that consumes `bundle::assigned_transcripts`.
