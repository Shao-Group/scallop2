bin_PROGRAMS = scallop2


GTF_INCLUDE = $(top_srcdir)/lib/gtf
UTIL_INCLUDE = $(top_srcdir)/lib/util
GRAPH_INCLUDE = $(top_srcdir)/lib/graph

GTF_LIB = $(top_builddir)/lib/gtf
UTIL_LIB = $(top_builddir)/lib/util
GRAPH_LIB = $(top_builddir)/lib/graph

scallop2_CPPFLAGS = -std=c++11 -I$(GTF_INCLUDE) -I$(GRAPH_INCLUDE) -I$(UTIL_INCLUDE)
scallop2_LDFLAGS = -pthread -L$(GTF_LIB) -L$(GRAPH_LIB) -L$(UTIL_LIB)
scallop2_LDADD = -lgtf -lgraph -lutil

scallop2_SOURCES = splice_graph.h splice_graph.cc \
				   super_graph.h super_graph.cc \
				   sgraph_compare.h sgraph_compare.cc \
				   vertex_info.h vertex_info.cc \
				   edge_info.h edge_info.cc \
				   interval_map.h interval_map.cc \
				   config.h config.cc \
				   hit.h hit.cc \
				   partial_exon.h partial_exon.cc \
				   hyper_set.h hyper_set.cc \
				   subsetsum.h subsetsum.cc \
				   router.h router.cc \
				   region.h region.cc \
				   junction.h junction.cc \
				   fragment.h fragment.cc \
				   bundle_base.h bundle_base.cc \
				   bundle_bridge.h bundle_bridge.cc \
				   bridger.h bridger.cc \
				   bundle.h bundle.cc \
				   fcluster.h fcluster.cc \
				   path.h path.cc \
				   equation.h equation.cc \
				   gtf.h gtf.cc \
				   scallop.h scallop.cc \
				   previewer.h previewer.cc \
				   assembler.h assembler.cc \
				   transcript_set.h transcript_set.cc \
				   reference.h reference.cc \
				   filter.h filter.cc \
				   main.cc 
