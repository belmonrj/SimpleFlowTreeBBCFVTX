This is a readme to run the flow analysis code.

The module, VTX_event_plane_reco, is running on the 
analysis taxi. The relevant taxi macro is here:
offline/AnalysisTrain/pat/macro/Run_VTX_event_plane.C

There are a number of switches in the macro which allow
you to run different parts of the code.

This module creates per event Ttrees.

To process these trees, one needs to run the module:
wrk/post_ana_bbc_fvtx.cpp

This is a compiled macro. You must change the input and
output file location in the macro. Also, there is a calibration
text file that needs to be written out. It takes in 2 arguments:
run number and calibration iteration. This should be run in
parallel over condor with one run number per job. It should
be run like this:
post_ana_bbc_fvtx(433107,1)
gROOT->Reset();
post_ana_bbc_fvtx(433107,2)
gROOT->Reset();
post_ana_bbc_fvtx(433107,3)

The output file creates TProfiles.
