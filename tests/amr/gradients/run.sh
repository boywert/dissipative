#!/bin/bash

rm -rf output/*
mpirun -np 4 ./Arepo.FD param.txt |& tee run.log
