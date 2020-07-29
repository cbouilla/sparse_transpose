# MKL
source /usr/intel/mkl/bin/mklvars.sh intel64
# source /usr/intel/mkl/bin/mklvars.sh intel64 ilp64  <--- si on veut des entiers 64 bits

# TBB
source /usr/intel/tbb/bin/tbbvars.sh intel64
# Print TBB version
export TBB_VERSION=0

# ICC / ICPC
source /usr/intel/compilers_and_libraries/linux/bin/compilervars.sh intel64
