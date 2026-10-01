#!/usr/bin/env bash

set -euo pipefail

src_dir=$(cd "$(dirname "$0")/.." && pwd)
test_dir=$(mktemp -d /tmp/scallop2-assigned-graph-tests.XXXXXX)
trap 'rm -rf "$test_dir"' EXIT

make -C "$src_dir" -j2 >/dev/null

g++ -std=c++11 -w \
	-I"$src_dir" \
	-I"$src_dir/../lib/gtf" \
	-I"$src_dir/../lib/graph" \
	-I"$src_dir/../lib/util" \
	-I/home/faculty/mxs2589/shared/tools/htslib/htslib-1.5-install/include \
	-I/home/faculty/mxs2589/shared/tools/boost/boost_1_70_0 \
	-c "$src_dir/tests/test_assigned_transcript_graph.cc" -o "$test_dir/test.o"

objects=()
for object in "$src_dir"/scallop2-*.o; do
	if [[ "$object" == *scallop2-main.o ]]; then continue; fi
	objects+=("$object")
done

g++ "$test_dir/test.o" "${objects[@]}" \
	-L"$src_dir/../lib/gtf" -L"$src_dir/../lib/graph" -L"$src_dir/../lib/util" \
	-L/home/faculty/mxs2589/shared/tools/htslib/htslib-1.5-install/lib \
	-lgtf -lgraph -lutil -lhts -pthread \
	-o "$test_dir/test_assigned_transcript_graph"

"$test_dir/test_assigned_transcript_graph"
