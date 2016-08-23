#!/bin/sh

g++ -o weight.exe weight.cpp -m32 -lstdc++ `root-config --libs` -I$ROOTSYS/include

