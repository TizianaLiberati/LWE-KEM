# Builder stage
FROM nvcr.io/nvidia/nvhpc:24.3-devel-cuda12.3-ubuntu22.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive

ARG openmp_branch=OpenMP
ARG openacc_branch=OpenACC

RUN apt-get update && \
    apt-get install -y libssl-dev wget git && \
    apt-get clean

# Install CMake
RUN wget -qO- https://github.com/Kitware/CMake/releases/download/v4.3.1/cmake-4.3.1-linux-$(uname -m).tar.gz \
    | tar --strip-components=1 -xz -C /usr/local

WORKDIR /workspace

# OpenMP version
RUN git clone -b $openmp_branch --depth 1 \
    https://github.com/TizianaLiberati/LWE-KEM OpenMP

WORKDIR /workspace/OpenMP/Codice

RUN make clean && \
    make

# OpenACC version
WORKDIR /workspace
RUN git clone -b $openacc_branch --depth 1 --recurse-submodules \
    https://github.com/TizianaLiberati/LWE-KEM OpenACC

WORKDIR /workspace/OpenACC/RNGonGPU
RUN . /etc/profile.d/lmod.sh && \
    module load nvhpc && \
    export NVHPC_CUDA_HOME=${NVHPC_ROOT}/cuda/12.3 && \
    cmake -S . -D CMAKE_CUDA_ARCHITECTURES="70;80;90;90-virtual" -B build && \
    cmake --build ./build/ -j

WORKDIR /workspace/OpenACC/Codice
RUN . /etc/profile.d/lmod.sh && \
    module load nvhpc && \
    export NVHPC_CUDA_HOME=${NVHPC_ROOT}/cuda/12.3 && \
    make clean && \
    make


# Runtime stage
FROM nvcr.io/nvidia/nvhpc:24.3-runtime-cuda12.3-ubuntu22.04

WORKDIR /workspace

COPY --from=builder /workspace/OpenMP/Codice ./OpenMP
COPY --from=builder /workspace/OpenACC/Codice ./OpenACC

CMD ["/bin/bash"] 
