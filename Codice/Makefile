#CXX      = nvc++
#CXXFLAGS = -std=c++20 -O2 -acc -gpu=cc90,managed -Minfo=accel
#LDFLAGS  = -lcrypto

#SRCS = Codice_modificato.cpp pke.cpp utils.cpp hash_openssl.cpp noise.cpp kem.cpp

#all: lwe_kem

#clean:
#	rm -f lwe_kem

#lwe_kem: $(SRCS) xorshift.h
#	$(CXX) $(CXXFLAGS) $(SRCS) $(LDFLAGS) -o lwe_kem

# ---------- benchmark ----------
#N_LIST = 10
#DIM_LIST = 512 1024 2048 4096

#bench: lwe_kem
#	@echo "n;avg_keygen_us;avg_encaps_us;avg_decaps_us;total_s;mismatches"
#	@for d in $(DIM_LIST); do \
	  echo ">>> n=$$d"; \
	  ./lwe_kem $(N_LIST) $$d; \
	done


CXX   = nvc++
CUDAC = nvcc

GPU_ARCH ?= 90

RNGONGPU_PREFIX = $(HOME)/local/rngongpu-aarch64
RNGONGPU_INC    = $(RNGONGPU_PREFIX)/include/RNGonGPU-1.0
GPUNTT_INC      = $(RNGONGPU_PREFIX)/include/GPUNTT-1.0
RNGONGPU_LIB    = $(RNGONGPU_PREFIX)/lib

NVHPC_ROOT = /home/tliberati/spack/opt/spack/linux-neoverse_v2/nvhpc-25.7-uhp3oqexp2v2hccjjtn5zxhhvwtdgxqs/Linux_aarch64/25.7
CUDA_HOME  = $(NVHPC_ROOT)/cuda/12.9
CUDA_LIB   = $(CUDA_HOME)/targets/sbsa-linux/lib
CUDA_INC   = $(CUDA_HOME)/targets/sbsa-linux/include

CURAND_LIB = $(NVHPC_ROOT)/math_libs/12.9/targets/sbsa-linux/lib

OPENSSL_LIB = /usr/lib/aarch64-linux-gnu
OPENSSL_INC = /usr/include

CXXFLAGS = -std=c++20 -O2 -acc -gpu=cc$(GPU_ARCH),managed -Minfo=accel
CUDAFLAGS = -O2 -arch=sm_$(GPU_ARCH) -I$(CUDA_INC) -I$(RNGONGPU_INC) -I$(GPUNTT_INC) -I$(OPENSSL_INC)

LDFLAGS = -acc -gpu=cc$(GPU_ARCH),managed \
          -L$(CUDA_LIB) -Wl,-rpath,$(CUDA_LIB) \
          -L$(CURAND_LIB) -Wl,-rpath,$(CURAND_LIB) \
          -L$(RNGONGPU_LIB) -Wl,-rpath,$(RNGONGPU_LIB) \
          -L$(OPENSSL_LIB) -Wl,-rpath,$(OPENSSL_LIB)

LIBS = -lrngongpu-1.0 -lntt-1.0 -lcudart -lcurand -lssl -lcrypto

CPP_SRCS = Codice_modificato.cpp pke.cpp utils.cpp hash_openssl.cpp noise.cpp kem.cpp
CU_SRCS  = rngongpu_adapter.cu

CPP_OBJS = $(CPP_SRCS:.cpp=.o)
CU_OBJS  = $(CU_SRCS:.cu=.o)

all: lwe_kem

clean:
	rm -f *.o lwe_kem

lwe_kem: $(CPP_OBJS) $(CU_OBJS)
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) $(LIBS) -o $@

%.o: %.cpp xorshift.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o: %.cu
	$(CUDAC) $(CUDAFLAGS) -c $< -o $@

# ---------- benchmark ----------
N_LIST = 10
DIM_LIST = 512 1024 2048 4096

bench: lwe_kem
	@echo "n;avg_keygen_us;avg_encaps_us;avg_decaps_us;total_s;mismatches"
	@for d in $(DIM_LIST); do \
	  echo ">>> n=$$d"; \
	  ./lwe_kem $(N_LIST) $$d; \
	done