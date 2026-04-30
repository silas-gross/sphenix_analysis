AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -isystem$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include \
  -isystem$(OPT_SPHENIX)/include  


AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS = \
  EtaShiftStudy.h \
  BuildMetaTowers.h

lib_LTLIBRARIES = \
  libEtaShiftStudy.la

libEtaShiftStudy_la_SOURCES = \
  EtaShiftStudy.cc

libEtaShiftStudy_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -lfun4all \
  -lcalo_io \
  -lg4jets \
  -lphhepmc \
  -lcalotrigger_io \
  -lcalotrigger \
  -lVandyClasses  

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libEtaShiftStudy.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
