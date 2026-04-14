# ---------------- Compiler Settings ----------------
CXX   = nvc++ -g
CUDAC = nvcc

GPU_ARCH ?= 80

# ---------------- NVTX Toggle ----------------
# Usage:
#   make                -> no NVTX
#   make nvtx=1         -> enable NVTX
#   make profile=off    -> disable NVTX explicitly

ifeq ($(nvtx),1)
    USE_NVTX = 1
endif

ifeq ($(profile),off)
    USE_NVTX =
endif

# ---------------- Include Paths ----------------
RNGONGPU_PREFIX = ../RNGonGPU
RNGONGPU_INC    = $(RNGONGPU_PREFIX)/src/include
GPUNTT_INC      = $(RNGONGPU_PREFIX)/thirdparty/GPU-NTT/src/include
OPENSSL_INC     = /usr/include

# ---------------- Flags ----------------
BASE_CXXFLAGS = -std=c++20 -O2 -acc -gpu=cc$(GPU_ARCH) -Minfo=accel \
                -I$(RNGONGPU_INC) -I$(GPUNTT_INC) -I$(OPENSSL_INC)

CUDAFLAGS = -O2 -arch=sm_$(GPU_ARCH) \
            -I$(RNGONGPU_INC) -I$(GPUNTT_INC) -I$(OPENSSL_INC)

# NVTX flags
ifeq ($(USE_NVTX),1)
    CXXFLAGS = $(BASE_CXXFLAGS) -DUSE_NVTX
    NVTX_LIB = -lnvToolsExt
else
    CXXFLAGS = $(BASE_CXXFLAGS)
    NVTX_LIB =
endif

# ---------------- Library Paths ----------------
RNGONGPU_LIB = $(RNGONGPU_PREFIX)/build/src
GPUNTT_LIB   = $(RNGONGPU_PREFIX)/build/thirdparty/GPU-NTT/src
OPENSSL_LIB  = /usr/lib/aarch64-linux-gnu

LDFLAGS = -acc -gpu=cc$(GPU_ARCH),managed -cuda \
          -L$(RNGONGPU_LIB) -L$(GPUNTT_LIB) -L$(OPENSSL_LIB)

LIBS = -lrngongpu-1.0 -lntt-1.0 -lssl -lcrypto $(NVTX_LIB)

# ---------------- Source Files ----------------
CPP_SRCS = Codice_modificato.cpp pke.cpp utils.cpp hash_openssl.cpp noise.cpp kem.cpp
CU_SRCS  = rngongpu_adapter.cu

CPP_OBJS = $(CPP_SRCS:.cpp=.o)
CU_OBJS  = $(CU_SRCS:.cu=.o)

# ---------------- Targets ----------------
all: lwe_kem

clean:
	rm -f *.o lwe_kem

lwe_kem: $(CPP_OBJS) $(CU_OBJS)
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) $(LIBS) -o $@

# Compile C++ sources
%.o: %.cpp xorshift.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile CUDA sources
%.o: %.cu
	$(CUDAC) $(CUDAFLAGS) -c $< -o $@

# ---------------- Benchmark ----------------
N_LIST = 10
DIM_LIST = 512 1024 2048 4096

bench: lwe_kem
	@echo "n;avg_keygen_us;avg_encaps_us;avg_decaps_us;total_s;mismatches"
	@for d in $(DIM_LIST); do \
	  echo ">>> n=$$d"; \
	  ./lwe_kem $(N_LIST) $$d; \
	done
