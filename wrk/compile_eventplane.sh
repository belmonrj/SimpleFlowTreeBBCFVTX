#!/bin/sh

g++ -o eventplane.exe eventplane.cpp -Wall -m32 -lstdc++ `root-config --libs` -I$ROOTSYS/include
