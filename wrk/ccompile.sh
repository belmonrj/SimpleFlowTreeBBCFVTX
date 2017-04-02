#!/bin/sh

#g++ -o cumulants.exe cumulants.cpp -Wall -m32 -lstdc++ `root-config --libs` -I$ROOTSYS/include
#g++ -o cumulants_pAu.exe cumulants_pAu.cpp -Wall -m32 -lstdc++ `root-config --libs` -I$ROOTSYS/include
g++ -o cumulants_3HeAu.exe cumulants_3HeAu.cpp -Wall -m32 -lstdc++ `root-config --libs` -I$ROOTSYS/include
g++ -o cumulants_3HeAu_class.exe cumulants_3HeAu_class.cpp ntp_event.cpp -Wall -m32 -lstdc++ `root-config --libs` -I$ROOTSYS/include

