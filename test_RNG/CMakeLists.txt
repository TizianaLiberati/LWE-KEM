cmake_minimum_required(VERSION 3.26)
project(acc_rng_probe LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)

message(STATUS "CMAKE_SYSTEM_PROCESSOR = ${CMAKE_SYSTEM_PROCESSOR}")

if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64|arm64)")
    set(RNGONGPU_PREFIX "$ENV{HOME}/local/rngongpu-aarch64")
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64)")
    set(RNGONGPU_PREFIX "$ENV{HOME}/local/rngongpu")
else()
    message(FATAL_ERROR "Architettura non gestita: ${CMAKE_SYSTEM_PROCESSOR}")
endif()

message(STATUS "Uso RNGonGPU da: ${RNGONGPU_PREFIX}")

# prova prima lib64, poi lib
if(EXISTS "${RNGONGPU_PREFIX}/lib64/cmake/RNGonGPU-1.0/RNGonGPUConfig.cmake")
    set(RNGonGPU_DIR "${RNGONGPU_PREFIX}/lib64/cmake/RNGonGPU-1.0")
    set(GPUNTT_DIR   "${RNGONGPU_PREFIX}/lib64/cmake/GPUNTT-1.0")
elseif(EXISTS "${RNGONGPU_PREFIX}/lib/cmake/RNGonGPU-1.0/RNGonGPUConfig.cmake")
    set(RNGonGPU_DIR "${RNGONGPU_PREFIX}/lib/cmake/RNGonGPU-1.0")
    set(GPUNTT_DIR   "${RNGONGPU_PREFIX}/lib/cmake/GPUNTT-1.0")
else()
    message(FATAL_ERROR "RNGonGPUConfig.cmake non trovato sotto ${RNGONGPU_PREFIX}")
endif()

message(STATUS "RNGonGPU_DIR = ${RNGonGPU_DIR}")
message(STATUS "GPUNTT_DIR   = ${GPUNTT_DIR}")

find_package(CUDAToolkit REQUIRED)
find_package(RNGonGPU REQUIRED CONFIG)

#add_library(rng_wrap OBJECT wrapper_rngongpu.cu)
#target_link_libraries(rng_wrap PRIVATE RNGonGPU::rngongpu CUDA::cudart)
#set_target_properties(rng_wrap PROPERTIES CUDA_SEPARABLE_COMPILATION ON CUDA_RESOLVE_DEVICE_SYMBOLS ON)

add_executable(acc_probe1
    acc_probe1.cpp
    #$<TARGET_OBJECTS:rng_wrap>
)

target_link_libraries(acc_probe1 PRIVATE RNGonGPU::rngongpu CUDA::cudart)
target_compile_options(acc_probe1 PRIVATE -acc -gpu=cc90 -static-nvidia)
target_link_options(acc_probe1 PRIVATE -acc -gpu=cc90 -static-nvidia)
