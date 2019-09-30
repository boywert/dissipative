#!/bin/bash

make CONFIG=ForceLawTests/TreePM/Config.sh BUILD_DIR=ForceLawTests/TreePM/build EXEC=ForceLawTests/TreePM/Arepo

make CONFIG=ForceLawTests/TreePMWithHighResRegion/Config.sh BUILD_DIR=ForceLawTests/TreePMWithHighResRegion/build EXEC=ForceLawTests/TreePMWithHighResRegion/Arepo

make CONFIG=ForceLawTests/TreeEwald/Config.sh BUILD_DIR=ForceLawTests/TreeEwald/build EXEC=ForceLawTests/TreeEwald/Arepo
