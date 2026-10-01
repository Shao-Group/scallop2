#!/usr/bin/env bash

set -euo pipefail

src_dir=$(cd "$(dirname "$0")/.." && pwd)
test_dir=$(mktemp -d /tmp/scallop2-transcript-tests.XXXXXX)
trap 'rm -rf "$test_dir"' EXIT

cd "$src_dir"

include_flags=(
	-std=c++11
	-ffunction-sections
	-fdata-sections
	-I.
	-I../lib/gtf
	-I../lib/graph
	-I../lib/util
	-I/home/faculty/mxs2589/shared/tools/htslib/htslib-1.5-install/include
	-I/home/faculty/mxs2589/shared/tools/boost/boost_1_70_0
)

sources=(
	transcript_index.cc
	transcript_match.cc
	bundle_base.cc
	splice_graph.cc
	vertex_info.cc
	edge_info.cc
	path.cc
)

for source_file in "${sources[@]}"; do
	g++ "${include_flags[@]}" -w -c "$source_file" -o "$test_dir/${source_file%.cc}.o"
done

g++ "${include_flags[@]}" -w -c tests/test_transcript_assignment.cc -o "$test_dir/test.o"
g++ -Wl,--gc-sections "$test_dir"/*.o \
	../lib/gtf/libgtf.a ../lib/graph/libgraph.a ../lib/util/libutil.a \
	-o "$test_dir/test_transcript_assignment"

"$test_dir/test_transcript_assignment"
