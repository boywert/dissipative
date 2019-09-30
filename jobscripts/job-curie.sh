#!/bin/bash
#MSUB -r MyJob                                 # Request name
#MSUB -n 32                                    # Number of cores to use
#MSUB -T 1800                                  # Elapsed time limit in seconds
#MSUB -o example_%I.o                          # Standard output. %I is the job id
#MSUB -e example_%I.e                          # Error output. %I is the job id
#MSUB -@ volker.springel@h-its.org:begin,end   # mail address
#MSUB -q standard                              # thin curie nodes
#MSUB -A ra0844                                # Project ID

set -x                                                                                                                       
cd ${BRIDGE_MSUB_PWD}                                                                                                        
echo ${BRIDGE_MSUB_PWD}
                  
export OMPI_MCA_btl_openib_use_eager_rdma=1
export OMPI_MCA_mpi_leave_pinned=1


ccc_mprun ./Arepo param.txt

