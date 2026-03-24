EXES = openacc_c_main

CXX   = nvc++
CUDAC = nvcc

GPU_ARCH ?= 90

#RNGONGPU_ROOT = ../third_party/RNGonGPU
#RNGONGPU_INC  = $(RNGONGPU_ROOT)/src/include

###
RNGONGPU_PREFIX = $(HOME)/local/rngongpu-aarch64
RNGONGPU_INC    = $(RNGONGPU_PREFIX)/include/RNGonGPU-1.0
GPUNTT_INC      = $(RNGONGPU_PREFIX)/include/GPUNTT-1.0
RNGONGPU_LIB    = $(RNGONGPU_PREFIX)/lib
###

NVHPC_ROOT = /home/tliberati/spack/opt/spack/linux-neoverse_v2/nvhpc-25.7-uhp3oqexp2v2hccjjtn5zxhhvwtdgxqs/Linux_aarch64/25.7
CUDA_HOME  = $(NVHPC_ROOT)/cuda/12.9
CUDA_LIB   = $(CUDA_HOME)/targets/sbsa-linux/lib
CUDA_INC   = $(CUDA_HOME)/targets/sbsa-linux/include
#
CURAND_LIB = $(NVHPC_ROOT)/math_libs/12.9/targets/sbsa-linux/lib
#
###
OPENSSL_LIB = /usr/lib/aarch64-linux-gnu
OPENSSL_INC = /usr/include
###

CXXFLAGS  = -O2 -acc -gpu=cc$(GPU_ARCH) -Minfo=accel
#CUDAFLAGS = -O2 -arch=sm_$(GPU_ARCH) -I$(CUDA_INC)
CUDAFLAGS = -O2 -arch=sm_$(GPU_ARCH) -I$(CUDA_INC) -I$(RNGONGPU_INC) -I$(GPUNTT_INC) -I$(OPENSSL_INC)
#LDFLAGS   = -acc -gpu=cc$(GPU_ARCH) -L$(CUDA_LIB) -Wl,-rpath,$(CUDA_LIB) -lcudart
LDFLAGS = -acc -gpu=cc$(GPU_ARCH) \
          -L$(CUDA_LIB) -Wl,-rpath,$(CUDA_LIB) \
          -L$(CURAND_LIB) -Wl,-rpath,$(CURAND_LIB) \
          -L$(RNGONGPU_LIB) -Wl,-rpath,$(RNGONGPU_LIB) \
          -L$(OPENSSL_LIB) -Wl,-rpath,$(OPENSSL_LIB)

#LIBS = -lrngongpu-1.0 -lntt-1.0 -lcudart
LIBS = -lrngongpu-1.0 -lntt-1.0 -lcudart -lcurand -lssl -lcrypto
all: $(EXES)

openacc_c_main: openacc_c_main.o saxpy_cuda.o
#openacc_c_main: openacc_c_main.o saxpy_cuda_rnd.o
#	$(CXX) -o $@ $^ $(LDFLAGS)
	$(CXX) -o $@ $^ $(LDFLAGS) $(LIBS)

openacc_c_main.o: openacc_c_main.cpp
	$(CXX) $(CXXFLAGS) -c $<

saxpy_cuda.o: saxpy_cuda.cu
#saxpy_cuda_rnd.o: saxpy_cuda_rnd.cu
	$(CUDAC) $(CUDAFLAGS) -c $<

clean:
	rm -f *.o $(EXES)