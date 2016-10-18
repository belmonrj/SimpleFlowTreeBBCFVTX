AUTOMAKE_OPTIONS = foreign

INCLUDES = \
	-I$(includedir) \
	-I$(OFFLINE_MAIN)/include \
	-I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libSimpleFlowTreeBBCFVTX.la

noinst_HEADERS = \
  SimpleFlowTreeBBCFVTX.h \
  SimpleFlowTreeBBCFVTXLinkDef.h

libSimpleFlowTreeBBCFVTX_la_SOURCES = \
  SimpleFlowTreeBBCFVTX.C \
  SimpleFlowTreeBBCFVTX_Dict.C

libSimpleFlowTreeBBCFVTX_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lSubsysReco \
  -lTOAD

# Rule for generating CINT dictionaries from class headers.
SimpleFlowTreeBBCFVTX_Dict.C: \
  SimpleFlowTreeBBCFVTX.h \
  SimpleFlowTreeBBCFVTXLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*
