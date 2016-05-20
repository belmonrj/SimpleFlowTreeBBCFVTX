#!/bin/sh

#g++ -std=c++11 -o flatten.exe flattening.cpp -lstdc++ `root-config --libs` -I$ROOTSYS/include
g++ -o flattening.exe flattening.cpp -m32 -lstdc++ `root-config --libs` -I$ROOTSYS/include

