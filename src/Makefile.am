bin_PROGRAMS = scallop2


GTF_INCLUDE = $(top_srcdir)/lib/gtf
UTIL_INCLUDE = $(top_srcdir)/lib/util

GTF_LIB = $(top_builddir)/lib/gtf
UTIL_LIB = $(top_builddir)/lib/util

scallop2_CPPFLAGS = -std=c++11 -I$(GTF_INCLUDE) -I$(UTIL_INCLUDE)
scallop2_LDFLAGS = -pthread -pg -L$(GTF_LIB) -L$(UTIL_LIB)
scallop2_LDADD = -lgtf -lutil

scallop2_SOURCES = interval_map.h interval_map.cc \
				   config.h config.cc \
				   hit.h hit.cc \
				   partial_exon.h partial_exon.cc \
				   region.h region.cc \
				   junction.h junction.cc \
				   fragment.h fragment.cc \
				   bundle_base.h bundle_base.cc \
				   bundle_bridge.h bundle_bridge.cc \
				   bridger.h bridger.cc \
				   fcluster.h fcluster.cc \
				   path.h path.cc \
				   previewer.h previewer.cc \
				   assembler.h assembler.cc \
				   transcript_set.h transcript_set.cc \
				   reference.h reference.cc \
				   filter.h filter.cc \
				   RO_read.h RO_read.cc \
				   main.cc 
