#!/bin/bash
curr_dir=$(pwd)
cd  ../../../
make CONFIG=$curr_dir/Config_static.sh EXEC=$curr_dir/Arepo_static BUILD_DIR=$curr_dir/.build_static
make CONFIG=$curr_dir/Config_amr.sh EXEC=$curr_dir/Arepo_amr BUILD_DIR=$curr_dir/.build_amr
