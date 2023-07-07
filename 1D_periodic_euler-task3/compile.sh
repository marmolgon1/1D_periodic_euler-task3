#!/bin/bash

SOURCE_FILES="\
main.cpp \
DataStructs.cpp \
rk4.cpp \
FluxFunctions.cpp \
RHSoperator.cpp"

mpic++ -pg -g -O3 -D_DOUBLE_ $SOURCE_FILES -Iincludes -o double.p

mpic++ -pg -g -O3 $SOURCE_FILES -Iincludes -o single.p
