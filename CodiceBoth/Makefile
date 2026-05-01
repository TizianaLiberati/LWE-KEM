# ============================================================================
#  Makefile  –  LWE-KEM GPU benchmark
#  Builds two executables:
#    - lwe_kem_acc  (OpenACC backend)
#    - lwe_kem_omp  (OpenMP Target backend)
# ============================================================================

CXX   = nvc++
CUDAC = nvcc

# sm_70 = V100 | sm_80 = A100 | sm_90 = H100/GH200
GENCODES = -gencode arch=compute_80,code=compute_80

# GPU target (shared)
GPU_TARGET = -gpu=cc80,managed

# Optional NVTX (same for both backends)
ifeq ($(nvtx),1)
    NVTX_DEFINE = -DUSE_NVTX
    NVTX_LIB    = -lnvToolsExt
endif

# OpenSSL
ifdef OPENSSL_DIR
    OPENSSL_INC     = -I$(OPENSSL_DIR)/include
    OPENSSL_LDFLAGS = -L$(OPENSSL_DIR)/lib -L$(OPENSSL_DIR)/lib64
endif

# RNGonGPU
RNGONGPU_PREFIX = ../RNGonGPU
RNGONGPU_INC    = -I$(RNGONGPU_PREFIX)/src/include
GPUNTT_INC      = -I$(RNGONGPU_PREFIX)/thirdparty/GPU-NTT/src/include

BASE_INC = $(RNGONGPU_INC) $(GPUNTT_INC) $(OPENSSL_INC)

# ---------------------------------------------------------------------------
#  Per-backend compiler / linker flag sets
# ---------------------------------------------------------------------------

# OpenACC backend
ACC_DEFINE = -DUSE_OPENACC
ACC_FLAG   = -acc
ACC_INFO   = -Minfo=accel

ACC_CXXFLAGS = \
    -std=c++20 -O2 \
    $(ACC_FLAG) $(GPU_TARGET) \
    $(ACC_INFO) \
    $(ACC_DEFINE) \
    $(NVTX_DEFINE) \
    $(BASE_INC)

ACC_LDFLAGS = \
    $(ACC_FLAG) $(GPU_TARGET) \
    -cuda \
    -L$(RNGONGPU_PREFIX)/build/src \
    -L$(RNGONGPU_PREFIX)/build/thirdparty/GPU-NTT/src \
    $(OPENSSL_LDFLAGS)

# OpenMP Target backend
OMP_DEFINE = -DUSE_OMP_TARGET
OMP_FLAG   = -mp=gpu
OMP_INFO   = -Minfo=mp

OMP_CXXFLAGS = \
    -std=c++20 -O2 \
    $(OMP_FLAG) $(GPU_TARGET) \
    $(OMP_INFO) \
    $(OMP_DEFINE) \
    $(NVTX_DEFINE) \
    $(BASE_INC)

OMP_LDFLAGS = \
    $(OMP_FLAG) $(GPU_TARGET) \
    -cuda \
    -L$(RNGONGPU_PREFIX)/build/src \
    -L$(RNGONGPU_PREFIX)/build/thirdparty/GPU-NTT/src \
    $(OPENSSL_LDFLAGS)

# CPU-only flags (for utils.cpp)
CXXFLAGS_CPU = \
    -std=c++20 -O2 \
    $(BASE_INC)

CUDAFLAGS = -O2 $(GENCODES) $(BASE_INC)

LIBS = -lrngongpu-1.0 -lntt-1.0 -lssl -lcrypto $(NVTX_LIB)

# ---------------------------------------------------------------------------
#  Sources and objects (shared)
# ---------------------------------------------------------------------------

CPP_SRCS_GPU = Codice_modificato.cpp pke.cpp hash_openssl.cpp noise.cpp kem.cpp
CPP_SRCS_CPU = utils.cpp
CPP_SRCS     = $(CPP_SRCS_GPU) $(CPP_SRCS_CPU)

CU_SRCS  = rngongpu_adapter.cu

# Objects are per backend to avoid mixing flags
ACC_OBJS_GPU = $(CPP_SRCS_GPU:.cpp=_acc.o)
OMP_OBJS_GPU = $(CPP_SRCS_GPU:.cpp=_omp.o)

# CPU object is shared (no backend-specific flags)
CPP_OBJS_CPU = $(CPP_SRCS_CPU:.cpp=.o)

CU_OBJS      = $(CU_SRCS:.cu=.o)

# ---------------------------------------------------------------------------
#  Phony targets
# ---------------------------------------------------------------------------

.PHONY: all clean info bench

all: lwe_kem_acc lwe_kem_omp
	@echo ""
	@echo "  Built: lwe_kem_acc  (OpenACC)"
	@echo "  Built: lwe_kem_omp  (OpenMP Target)"
	@echo "  GPU target : $(GPU_TARGET)"
	@echo "  NVTX       : $(if $(nvtx),enabled,disabled)"
	@echo ""

# ---------------------------------------------------------------------------
#  Link targets
# ---------------------------------------------------------------------------

lwe_kem_acc: $(ACC_OBJS_GPU) $(CPP_OBJS_CPU) $(CU_OBJS)
	$(CXX) $(ACC_CXXFLAGS) $^ $(ACC_LDFLAGS) $(LIBS) -o $@

lwe_kem_omp: $(OMP_OBJS_GPU) $(CPP_OBJS_CPU) $(CU_OBJS)
	$(CXX) $(OMP_CXXFLAGS) $^ $(OMP_LDFLAGS) $(LIBS) -o $@

# ---------------------------------------------------------------------------
#  Compile rules
# ---------------------------------------------------------------------------

# OpenACC GPU-enabled C++ sources (depend on gpu_backend.h)
$(ACC_OBJS_GPU): %_acc.o: %.cpp gpu_backend.h
	$(CXX) $(ACC_CXXFLAGS) -c $< -o $@

# OpenMP Target GPU-enabled C++ sources (depend on gpu_backend.h)
$(OMP_OBJS_GPU): %_omp.o: %.cpp gpu_backend.h
	$(CXX) $(OMP_CXXFLAGS) -c $< -o $@

# CPU-only C++ sources (no GPU/offload flags)
$(CPP_OBJS_CPU): %.o: %.cpp
	$(CXX) $(CXXFLAGS_CPU) -c $< -o $@

# CUDA sources (shared)
%.o: %.cu
	$(CUDAC) $(CUDAFLAGS) -c $< -o $@

# ---------------------------------------------------------------------------
#  Utility targets
# ---------------------------------------------------------------------------

clean:
	rm -f *.o lwe_kem_acc lwe_kem_omp

info:
	@echo "ACC_CXXFLAGS : $(ACC_CXXFLAGS)"
	@echo "ACC_LDFLAGS  : $(ACC_LDFLAGS)"
	@echo "OMP_CXXFLAGS : $(OMP_CXXFLAGS)"
	@echo "OMP_LDFLAGS  : $(OMP_LDFLAGS)"
	@echo "CXXFLAGS_CPU : $(CXXFLAGS_CPU)"

# Simple benchmark uses the OpenACC binary by default
N_LIST   = 10
DIM_LIST = 512 1024 2048 4096

bench: lwe_kem_acc
	@echo "n;avg_keygen_us;avg_encaps_us;avg_decaps_us;total_s;mismatches"
	@for d in $(DIM_LIST); do \
	    echo ">>> n=$$d  backend=OpenACC"; \
	    ./lwe_kem_acc $(N_LIST) $$d; \
	done
