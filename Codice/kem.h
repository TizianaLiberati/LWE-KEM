#pragma once
#include <vector>
#include <cstdint>

/* Forward declaration */
struct DevicePool;

/* ══════════════════════════════════════════════════════════════════════════
 *  GPU KEM API — v4 (aligned with OPT-5 / OPT-6 / OPT-8)
 * ══════════════════════════════════════════════════════════════════════════ */

/* ── Encapsulation ──────────────────────────────────────────────────────── */
/* OPT-5: A_flat is uint16_t*
 * OPT-8: r_buf_T / e1_buf_T are l-major (transposed)
 * OPT-6: uses pool.AT internally
 */
void Encaps_GPU_Aflat(uint64_t key_seed,
                      uint32_t n,
                      uint32_t q,
                      uint16_t* A_flat,        // ✅ FIXED TYPE
                      int32_t*  t_ptr,
                      int32_t*  c_out,
                      std::vector<int32_t>& Hash_K,
                      int32_t*  m_buf,
                      int32_t*  r_buf_T,
                      int32_t*  e1_buf_T,
                      int32_t*  e2_buf,
                      int32_t*  ptxt_buf,
                      DevicePool& pool);

/* ── Decapsulation ──────────────────────────────────────────────────────── */
void Decaps_GPU_Aflat(uint64_t key_seed,
                      uint32_t n,
                      uint32_t q,
                      uint16_t* A_flat,        // ✅ FIXED TYPE
                      int32_t*  t_ptr,
                      int32_t*  s_ptr,
                      int32_t*  z_ptr,
                      int32_t*  c_in,
                      std::vector<int32_t>& Hash_K,
                      int32_t*  mp_buf,
                      int32_t*  dec_buf,
                      int32_t*  r_buf_T,
                      int32_t*  e1_buf_T,
                      int32_t*  e2_buf,
                      int32_t*  ptxt_buf,
                      int32_t*  c_chk,
                      DevicePool& pool);

/* ── CPU fallback stubs (unchanged) ─────────────────────────────────────── */
void Encaps(uint32_t, uint32_t,
            std::vector<int32_t>&,
            std::vector<int32_t>&,
            const std::vector<std::vector<int32_t>>&,
            const std::vector<std::vector<int32_t>>&,
            std::vector<int32_t>&);

void Decaps(uint32_t, uint32_t,
            const std::vector<int32_t>&,
            const std::vector<int32_t>&,
            const std::vector<int32_t>&,
            std::vector<int32_t>&,
            const std::vector<std::vector<int32_t>>&,
            const std::vector<std::vector<int32_t>>&);
