#!/bin/bash
curr_dir=$(pwd)
cd  ../../../
make CONFIG=$curr_dir/Config.sh EXEC=$curr_dir/Arepo BUILD_DIR=$curr_dir/.build
