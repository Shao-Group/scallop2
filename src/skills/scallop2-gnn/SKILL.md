---
name: scallop2-gnn
description: Continue design, implementation, or review of transcript-to-bundle assignment and GNN data generation on the Scallop2 gnn branch.
---

# Scallop2 GNN module

Work from `src/` on the `gnn` branch. Read [`../../GNN_MODULE_STATUS.md`](../../GNN_MODULE_STATUS.md) before changing code; update its status and validation sections after meaningful work.

Preserve user changes. In particular, inspect `git status` before editing and do not overwrite an existing modification to `main.cc`.

## Stable codebase facts

- `assembler` streams a coordinate-sorted BAM, forms plus/minus `bundle_base` objects, and processes bundles in batches.
- `bundle` prepares intervals, junctions, regions, and partial exons in its constructor. `bundle::build(1, true)` builds and revises its `splice_graph` and hyper-set.
- GTF exon coordinates and BAM/graph coordinates are all 0-based, half-open intervals `[left, right)`. `item::parse()` converts GTF's 1-based inclusive start by decrementing the start only.
- `bundle_base::{lpos,rpos,strand,chrm}` describe the bundle. After graph construction, `splice_graph::{chrm,strand}` and interior `vertex_info::{lpos,rpos}` describe the graph.
- A graph splice junction is an edge between interior vertices whose source `rpos` is less than the target `lpos`; adjacent vertices with equal boundaries are contiguous exon pieces, not junctions.
- The library `genome::locate_gene()` linearly scans all genes and returns only the single gene with greatest genomic overlap. It is unsuitable for retrieving all candidate transcripts.
- `transcript_index` partitions immutable transcript records by chromosome/strand and uses a split interval map with `set<int>` codomains to fetch overlapping IDs.
- `transcript_matcher` extracts exon unions and exact splice junctions from the final revised graph. `assembler::process_gnn()` writes accepted IDs to `bundle::assigned_transcripts`.

## Assignment semantics

Unless the user revises the criterion, assign a transcript to a bundle only when all are true:

1. Chromosome and strand exactly match the final bundle/graph chromosome and strand.
2. The union of transcript-exon bases overlapping the union of interior graph-vertex intervals covers at least half of the transcript's exonic length: `2 * shared_exonic_bases >= transcript_exonic_bases`.
3. At least half of the transcript's exact introns occur as graph splice junctions: `2 * matched_junctions >= transcript_junctions`.

Treat the junction condition as vacuously true for a single-exon transcript, so exon overlap decides it. Count union lengths to avoid double counting, and use integer cross-multiplication to avoid threshold rounding ambiguity. Do not consume a transcript after a match: assignment is many-to-many unless the user requests unique ownership.

Use the final revised graph for exact scoring so labels describe the graph passed to the GNN. Fetch candidates first with a chromosome/strand interval index over transcript genomic bounds, then apply exact exon and junction scoring. Keep thresholds named/configurable even if their initial defaults are both `0.5`.

## Implementation boundaries

Keep indexing and matching separate: the index only narrows candidates by chromosome, strand, and genomic-bound overlap; a matcher computes exact criteria from graph intervals and junctions. Store stable flattened transcript records or stable numeric references rather than pointers into containers that may reallocate.

Add focused tests for coordinate conversion, boundary-touch non-overlap, multi-exon union overlap, odd junction counts, single-exon handling, strand/chromosome rejection, and multiple overlapping candidates. Also run a clean or forced rebuild when the parent build tree is writable. In a restricted `src/` workspace, use the compile flags in `src/Makefile` with `g++ -fsyntax-only` for changed translation units and record that full Automake regeneration remains pending.
