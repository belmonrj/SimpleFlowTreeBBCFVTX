AUTOMAKE_OPTIONS = foreign

INCLUDES = \
	-I$(includedir) \
	-I$(OFFLINE_MAIN)/include \
	-I$(ROOTSYS)/include 

lib_LTLIBRARIES = \
  libVTX_event_plane.la 

noinst_HEADERS = \
  VTX_event_plane_reco.h \
  VTX_event_planeLinkDef.h

libVTX_event_plane_la_SOURCES = \
  VTX_event_plane_reco.C \
  VTX_event_plane_Dict.C 

libVTX_event_plane_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lSubsysReco \
  -lTOAD

# Rule for generating CINT dictionaries from class headers.
VTX_event_plane_Dict.C: \
  VTX_event_plane_reco.h \
  VTX_event_planeLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*
