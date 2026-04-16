# ---------------- Compiler Settings ----------------
CXX   = nvc++
CUDAC = nvcc

# ---------------- GPU Architectures ----------------
# Targets: V100 (70), A100 (80), H100/GH200 (90)
# Include compute_90 PTX as a virtual arch for future-proofing
GENCODES = -gencode arch=compute_70,code=sm_70 \
           -gencode arch=compute_80,code=sm_80 \
           -gencode arch=compute_90,code=sm_90 \
           -gencode arch=compute_90,code=compute_90

NVC_GPU_FLAGS = -acc -gpu=cc70,cc80,cc90,managed

# ---------------- Make Options ----------------
# Usage:
#   make                -> no NVTX
#   make nvtx=1         -> enable NVTX
#   make profile=off    -> disable NVTX explicitly
#   make info=1         -> show accelerator info

ifeq ($(nvtx),1)
    USE_NVTX = 1
endif

ifeq ($(profile),off)
    USE_NVTX =
endif

ifeq ($(info),1)
    NVC_GPU_FLAGS += -Minfo=accel
endif

# ---------------- OpenSSL Configuration ----------------
# Use OPENSSL_DIR if provided, or fall back to system defaults.
ifdef OPENSSL_DIR
    OPENSSL_INC     = -I$(OPENSSL_DIR)/include
    OPENSSL_LDFLAGS = -L$(OPENSSL_DIR)/lib -L$(OPENSSL_DIR)/lib64
endif

# ---------------- Include Paths ----------------
RNGONGPU_PREFIX = ../RNGonGPU
RNGONGPU_INC    = -I$(RNGONGPU_PREFIX)/src/include
GPUNTT_INC      = -I$(RNGONGPU_PREFIX)/thirdparty/GPU-NTT/src/include
BASE_INC = $(RNGONGPU_INC) $(GPUNTT_INC) $(OPENSSL_INC)
# ---------------- Flags ----------------
CXXFLAGS  = -std=c++20 -O2 $(NVC_GPU_FLAGS) $(BASE_INC)
CUDAFLAGS = -O2 $(GENCODES) $(BASE_INC)

ifeq ($(USE_NVTX),1)
    CXXFLAGS += -DUSE_NVTX
    NVTX_LIB  = -lnvToolsExt
endif

# ---------------- Library Paths ----------------
RNGONGPU_LIB = -L$(RNGONGPU_PREFIX)/build/src
GPUNTT_LIB   = -L$(RNGONGPU_PREFIX)/build/thirdparty/GPU-NTT/src

LDFLAGS = $(NVC_GPU_FLAGS) -cuda $(RNGONGPU_LIB) $(GPUNTT_LIB) $(OPENSSL_LDFLAGS)
LIBS    = -lrngongpu-1.0 -lntt-1.0 -lssl -lcrypto $(NVTX_LIB)

# ---------------- Files ----------------
CPP_SRCS = Codice_modificato.cpp pke.cpp utils.cpp hash_openssl.cpp noise.cpp kem.cpp
CU_SRCS  = rngongpu_adapter.cu

CPP_OBJS = $(CPP_SRCS:.cpp=.o)
CU_OBJS  = $(CU_SRCS:.cu=.o)

# ---------------- Targets ----------------
.PHONY: all clean

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
