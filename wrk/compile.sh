#!/bin/sh

g++ -o flattening.exe flattening.cpp -Wall -m32 -lstdc++ `root-config --libs` -I$ROOTSYS/include
