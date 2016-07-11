#!/bin/sh

g++ -o flattening.exe flattening.cpp -m32 -lstdc++ `root-config --libs` -I$ROOTSYS/include
