#!/bin/bash
curr_dir=$(pwd)
cd  ../../../
make -j 8 CONFIG=$curr_dir/Config.sh EXEC=$curr_dir/Arepo BUILD_DIR=$curr_dir/.build
make -j 8 CONFIG=$curr_dir/Config_arepo.sh EXEC=$curr_dir/Arepo_arepo BUILD_DIR=$curr_dir/.build_arepo
make -j 8 CONFIG=$curr_dir/Config_albada.sh EXEC=$curr_dir/Arepo_albada BUILD_DIR=$curr_dir/.build_albada
make -j 8 CONFIG=$curr_dir/Config_superbee.sh EXEC=$curr_dir/Arepo_superbee BUILD_DIR=$curr_dir/.build_superbee
make -j 8 CONFIG=$curr_dir/Config_minbee.sh EXEC=$curr_dir/Arepo_minbee BUILD_DIR=$curr_dir/.build_minbee
make -j 8 CONFIG=$curr_dir/Config_vanleer.sh EXEC=$curr_dir/Arepo_vanleer BUILD_DIR=$curr_dir/.build_vanleer
