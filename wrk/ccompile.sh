#!/bin/sh

g++ -o cumulants.exe cumulants.cpp -Wall -m32 -lstdc++ `root-config --libs` -I$ROOTSYS/include
