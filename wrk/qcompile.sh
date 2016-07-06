#!/bin/sh

g++ -o qsimple.exe qsimple.cpp -m32 -lstdc++ `root-config --libs` -I$ROOTSYS/include

