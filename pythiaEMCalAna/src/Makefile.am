AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS = \
  pythiaEMCalAna.h

lib_LTLIBRARIES = \
  libpythiaEMCalAna.la

libpythiaEMCalAna_la_SOURCES = \
  pythiaEMCalAna.cc

libpythiaEMCalAna_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -lHepMC \
  -lCLHEP \
  -lcalo_io \
  -lg4eval \
  -lg4jets \
  -lcalotrigger_io

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libpythiaEMCalAna.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
