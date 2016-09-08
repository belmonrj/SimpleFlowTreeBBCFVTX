#!/bin/sh

g++ -o minbias.exe minbias.cpp -m32 -lstdc++ `root-config --libs` -I$ROOTSYS/include

