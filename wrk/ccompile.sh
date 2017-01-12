#!/bin/sh

g++ -o cumulants.exe cumulants.cpp -Wall -m32 -lstdc++ `root-config --libs` -I$ROOTSYS/include
g++ -o cumulants_pAu.exe cumulants_pAu.cpp -Wall -m32 -lstdc++ `root-config --libs` -I$ROOTSYS/include

