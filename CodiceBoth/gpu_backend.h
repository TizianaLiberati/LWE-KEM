/*  gpu_backend.h  –  Compile-time switch between OpenACC and OpenMP Target
 *
 *  Usage:
 *    -DUSE_OPENACC    → all pragmas expand to OpenACC
 *    -DUSE_OMP_TARGET → all pragmas expand to OpenMP Target offload
 *
 *  If neither is defined, the build fails with a clear error.
 *
 *  MAPPING TABLE
 *  ─────────────────────────────────────────────────────────────────────
 *  Concept                  OpenACC                  OpenMP Target
 *  ─────────────────────────────────────────────────────────────────────
 *  Persistent alloc         enter data create        target enter data alloc
 *  Persistent free          exit  data delete        target exit  data release
 *  Update host←device       update self              target update from
 *  Update device←host       update device            target update to
 *  Data region (scoped)     data present/create      target data map(...)
 *  Parallel loop (1-D)      parallel loop            target teams distribute parallel for
 *  Parallel loop (gang/vec) parallel loop gang vector  target teams distribute parallel for
 *  Inner vector loop        loop vector              simd / parallel for simd
 *  Inner reduction loop     loop vector reduction    parallel for reduction
 *  Async queue N            async(N)                 nowait + depend (simplified: nowait)
 *  Wait queues              wait(N,M)                taskwait
 *  Host-side device ptr     host_data use_device     use_device_ptr
 *  present() clause         present(ptr[0:len])      is_device_ptr(ptr) or map(present,...)
 *  ─────────────────────────────────────────────────────────────────────
 *
 *  NOTES ON ASYNC / NOWAIT
 *  OpenACC async(N) assigns work to numbered queues; OpenMP's equivalent
 *  is `nowait` on individual target regions combined with `taskwait`.
 *  Because NVHPC maps both onto CUDA streams, the generated PTX is
 *  functionally equivalent when the same compiler is used.
 *
 *  NOTES ON present() vs is_device_ptr()
 *  OpenACC `present()` asserts the buffer is already on-device (error if
 *  not). The closest OpenMP equivalent when using managed/unified memory
 *  is `map(present,alloc:ptr[0:N])`, which is compiler-hint only.
 *  For NVHPC managed memory (`-gpu=managed`) both models see the same
 *  physical pages, so the difference is only in diagnostic strictness.
 */

#pragma once

#if !defined(USE_OPENACC) && !defined(USE_OMP_TARGET)
#  error "Define exactly one of USE_OPENACC or USE_OMP_TARGET at compile time."
#endif

/* ======================================================================
 *  OpenACC backend
 * ====================================================================== */
#ifdef USE_OPENACC

/* --- Persistent device allocation / deallocation -------------------- */
#define GPU_ENTER_DATA_CREATE(clauses)   _Pragma("acc enter data " #clauses)
#define GPU_EXIT_DATA_DELETE(clauses)    _Pragma("acc exit data delete" #clauses)

/* --- Explicit data transfers ---------------------------------------- */
#define GPU_UPDATE_FROM_DEVICE(clauses)  _Pragma("acc update self" #clauses)
#define GPU_UPDATE_TO_DEVICE(clauses)    _Pragma("acc update device" #clauses)

/* --- Scoped data regions -------------------------------------------- */
/* Use as: GPU_DATA_REGION_BEGIN(present(p[0:n])) { ... } GPU_DATA_REGION_END */
#define GPU_DATA_REGION_BEGIN(clauses)   _Pragma("acc data " #clauses)  {
#define GPU_DATA_REGION_END              }

/* --- Host-side native device pointer -------------------------------- */
#define GPU_HOST_DATA_BEGIN(ptrs)        _Pragma("acc host_data use_device" #ptrs)  {
#define GPU_HOST_DATA_END                }

/* --- Parallel compute kernels --------------------------------------- */
/* 1-D loop with present clause */
#define GPU_PARALLEL_LOOP(clauses) \
    _Pragma("acc parallel loop " #clauses)

/* Gang+vector 1-D loop (outer) */
#define GPU_PARALLEL_LOOP_GANG_VEC(clauses) \
    _Pragma("acc parallel loop gang vector " #clauses)

/* Inner vector loop (no reduction) */
#define GPU_INNER_LOOP_VEC \
    _Pragma("acc loop vector")

/* Inner vector loop with reduction */
#define GPU_INNER_LOOP_VEC_RED(red) \
    _Pragma("acc loop vector reduction(" #red ")")

/* --- Async / synchronisation ---------------------------------------- */
#define GPU_ASYNC(queue)   async(queue)          /* appended to pragma clause */
#define GPU_WAIT(queues)   _Pragma("acc wait" #queues)

/* --- present() helper ----------------------------------------------- */
/* Use directly inside pragma clause strings, e.g.:
 *   GPU_PARALLEL_LOOP(GPU_PRESENT(p[0:n]) firstprivate(n))  */
#define GPU_PRESENT(ptrs)  present(ptrs)

#endif /* USE_OPENACC */


/* ======================================================================
 *  OpenMP Target backend
 * ====================================================================== */
#ifdef USE_OMP_TARGET

/* --- Persistent device allocation / deallocation -------------------- */
#define GPU_ENTER_DATA_CREATE(clauses)   _Pragma("omp target enter data " #clauses)
#define GPU_EXIT_DATA_DELETE(clauses)    _Pragma("omp target exit data " #clauses)

/* --- Explicit data transfers ---------------------------------------- */
#define GPU_UPDATE_FROM_DEVICE(clauses)  _Pragma("omp target update from" #clauses)
#define GPU_UPDATE_TO_DEVICE(clauses)    _Pragma("omp target update to" #clauses)

/* --- Scoped data regions -------------------------------------------- */
#define GPU_DATA_REGION_BEGIN(clauses)   _Pragma("omp target data " #clauses)  {
#define GPU_DATA_REGION_END              }

/* --- Host-side native device pointer -------------------------------- */
/* OpenMP: use_device_ptr goes on the target data directive            */
#define GPU_HOST_DATA_BEGIN(ptrs)        _Pragma("omp target data use_device_ptr" #ptrs)  {
#define GPU_HOST_DATA_END                }

/* --- Parallel compute kernels --------------------------------------- */
#define GPU_PARALLEL_LOOP(clauses) \
    _Pragma("omp target teams distribute parallel for " #clauses)

#define GPU_PARALLEL_LOOP_GANG_VEC(clauses) \
    _Pragma("omp target teams distribute parallel for " #clauses)

/* Inner simd loop (inside a teams distribute parallel for) */
#define GPU_INNER_LOOP_VEC \
    _Pragma("omp simd")

#define GPU_INNER_LOOP_VEC_RED(red) \
    _Pragma("omp parallel for reduction(" #red ")")

/* --- Async / synchronisation ---------------------------------------- */
/* OpenMP nowait requires a matching taskwait; we omit nowait for      */
/* correctness unless the caller knows to insert taskwait.             */
#define GPU_ASYNC(queue)   /* nowait — omitted for safety */
#define GPU_WAIT(queues)   _Pragma("omp taskwait")

/* --- map(present,...) helper ---------------------------------------- */
/* Closest OMP equivalent: hint that memory is already mapped.         */
#define GPU_PRESENT(ptrs)  map(present, alloc: ptrs)

#endif /* USE_OMP_TARGET */
