cmake_minimum_required(VERSION 3.26)
project(acc_rng LANGUAGES CXX CUDA)

find_package(CUDAToolkit REQUIRED)
find_package(RNGonGPU REQUIRED)

add_executable(acc_rng main.cpp copy_kernel.cu)

target_link_libraries(acc_rng
    PRIVATE
    RNGonGPU::rngongpu
    CUDA::cudart
)

set_target_properties(acc_rng PROPERTIES CUDA_SEPARABLE_COMPILATION ON)