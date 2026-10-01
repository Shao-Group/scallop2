---
name: scallop2-gnn
description: Continue design, implementation, or review of transcript-to-bundle assignment and GNN data generation on the Scallop2 gnn branch.
---

# Scallop2 GNN module

Work from `src/` on the `gnn` branch. Read [`../../GNN_MODULE_STATUS.md`](../../GNN_MODULE_STATUS.md) before changing code; update its status and validation sections after meaningful work.

Preserve user changes. In particular, inspect `git status` before editing and do not overwrite an existing modification to `main.cc`.

## Stable codebase facts

- `assembler` streams a coordinate-sorted BAM, forms plus/minus `bundle_base` objects, and processes bundles in batches.
- `process_gnn()` matches transcripts directly against `bundle_base`, then passes accepted transcript models into `bundle` for splice-graph construction.
- GTF exon coordinates and BAM/graph coordinates are all 0-based, half-open intervals `[left, right)`. `item::parse()` converts GTF's 1-based inclusive start by decrementing the start only.
- `bundle_base::{lpos,rpos,strand,chrm}` describe the bundle, `mmap` contains split covered intervals, and each `hit::spos` entry packs an observed donor/acceptor pair.
- The library `genome::locate_gene()` linearly scans all genes and returns only the single gene with greatest genomic overlap. It is unsuitable for retrieving all candidate transcripts.
- `transcript_index` partitions immutable transcript records by chromosome/strand and uses a split interval map with `set<int>` codomains to fetch overlapping IDs.
- `transcript_matcher` joins `bundle_base::mmap` into the bundle exon union and stores unique `int32_t` donor/acceptor positions extracted from `hit::spos`. `assembler::process_gnn()` increments deterministic per-transcript bundle counts and supplies accepted transcript models to `bundle`.
- `bundle` adds only the assigned-exon bases absent from read-derived structural coverage, retains those newly added intervals, inserts missing annotation junctions without duplicating read junctions, and partitions regions at assigned transcript boundaries. Each partial exon and graph vertex records its number of newly added bases. After revision it restores and validates every assigned transcript as a source-to-sink graph path.

## Assignment semantics

Unless the user revises the criterion, assign a transcript to a bundle only when all are true:

1. Chromosome and strand exactly match the bundle chromosome and inferred strand.
2. The union of transcript-exon bases overlapping the joined `bundle_base::mmap` intervals covers at least half of the transcript's exonic length.
3. At least half of the transcript's unique donor/acceptor positions occur among the bundle's observed splice positions.

Treat the splice-position condition as vacuously true for a single-exon transcript, so exon overlap decides it. Count union lengths and unique positions to avoid double counting. Do not consume a transcript after a match: assignment is many-to-many unless the user requests unique ownership.

Fetch candidates first with a chromosome/strand interval index over transcript genomic bounds, then apply exact exon and splice-position scoring against `bundle_base`. Keep thresholds named/configurable even if their initial defaults are both `0.5`. The splice-position option is `--min_bundle_transcript_splicing_position_overlap`; accept the former junction-named option only as a compatibility alias.

## Implementation boundaries

Keep indexing and matching separate: the index only narrows candidates by chromosome, strand, and genomic-bound overlap; a matcher computes exact criteria from bundle coverage and splice positions. Store stable flattened transcript records or stable numeric references rather than pointers into containers that may reallocate.

Preserve read-supported junction counts when an annotation junction already exists. Because matching permits partial exon overlap, include only the missing parts of assigned exon intervals when constructing partial exons; boundaries and junctions alone cannot create a missing exon vertex. Preserve the missing-interval union and its overlap length on each partial exon and graph vertex rather than inferring annotation-derived structure from edge weight. Reapply assigned paths after graph revision so pruning cannot invalidate a label.

Add focused tests for coordinate conversion, boundary-touch non-overlap, multi-exon union overlap, odd splice-position counts, single-exon handling, strand/chromosome rejection, and multiple overlapping candidates. Also run a clean or forced rebuild when the parent build tree is writable. In a restricted `src/` workspace, use the compile flags in `src/Makefile` with `g++ -fsyntax-only` for changed translation units and record that full Automake regeneration remains pending.
