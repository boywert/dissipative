#!/bin/bash
curr_dir=$(pwd)
cd  ../../../
make CONFIG=$curr_dir/Config_2d.sh EXEC=$curr_dir/Arepo_2d BUILD_DIR=$curr_dir/.build_2d
make CONFIG=$curr_dir/Config_2d_noamr.sh EXEC=$curr_dir/Arepo_2d_noamr BUILD_DIR=$curr_dir/.build_2d_noamr
make CONFIG=$curr_dir/Config_3d.sh EXEC=$curr_dir/Arepo_3d BUILD_DIR=$curr_dir/.build_3d
