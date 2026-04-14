/*  Codice_modificato.cpp  –  LWE-KEM GPU benchmark driver (NVTX optional)
 *
 *  Compile (no NVTX):
 *    nvc++ ... Codice_modificato.cpp ...
 *
 *  Compile (with NVTX):
 *    nvc++ -DUSE_NVTX ... Codice_modificato.cpp ... -lnvToolsExt
 */

#include <iostream>
#include <vector>
#include <cstdint>
#include <cstdlib>
#include <chrono>

#ifdef USE_NVTX
  #include <nvtx3/nvToolsExt.h>

  static constexpr uint32_t nvtxColor(uint8_t r, uint8_t g, uint8_t b)
  {
      return (0xFFu << 24) | (r << 16) | (g << 8) | b;
  }

  struct NvtxRange {
      NvtxRange(const char* name, uint32_t color)
      {
          nvtxEventAttributes_t attr = {};
          attr.version = NVTX_VERSION;
          attr.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
          attr.colorType = NVTX_COLOR_ARGB;
          attr.color = color;
          attr.messageType = NVTX_MESSAGE_TYPE_ASCII;
          attr.message.ascii = name;
          nvtxRangePushEx(&attr);
      }
      ~NvtxRange() { nvtxRangePop(); }
  };

#else
  struct NvtxRange {
      NvtxRange(const char*, uint32_t) {}
  };

  static constexpr uint32_t nvtxColor(uint8_t, uint8_t, uint8_t)
  {
      return 0;
  }
#endif

#include "pke.h"
#include "kem.h"
#include "utils.h"
#include <openssl/rand.h>

/* Colors */
constexpr uint32_t C_TOTAL    = nvtxColor(180,180,180);
constexpr uint32_t C_KEYGEN   = nvtxColor(80,160,255);
constexpr uint32_t C_ENCAPS   = nvtxColor(80,200,120);
constexpr uint32_t C_DECAPS   = nvtxColor(255,140,80);
constexpr uint32_t C_A_MGMT   = nvtxColor(200,120,255);
constexpr uint32_t C_PREFETCH = nvtxColor(150,150,255);

/* ================= Double Buffer ================= */
struct DoubleBufA {
    std::vector<uint16_t> ping_buf, pong_buf;
    uint16_t* ping=nullptr;
    uint16_t* pong=nullptr;

    uint64_t next_rho=0;
    uint32_t n_=0,q_=0;
    bool prefetch_valid=false;

    void init(uint32_t n,uint32_t q,uint64_t first_rho){
        NvtxRange r("DoubleBufA::init",C_A_MGMT);

        n_=n; q_=q;
        size_t An=(size_t)n*n;

        ping_buf.assign(An,0);
        pong_buf.assign(An,0);
        ping=ping_buf.data();
        pong=pong_buf.data();

        #pragma acc enter data create(ping[0:An],pong[0:An])

        {
            NvtxRange r2("GenerateA_initial",C_A_MGMT);
            GenerateA_GPU_async(first_rho,n,q,ping,STREAM_A);
            #pragma acc wait(STREAM_A)
        }

        next_rho=first_rho;
        prefetch_valid=false;
    }

    void prefetch_next(uint64_t rho_next){
        NvtxRange r("A_prefetch_async",C_PREFETCH);
        next_rho=rho_next;
        GenerateA_GPU_async(rho_next,n_,q_,pong,STREAM_PREFETCH);
        prefetch_valid=true;
    }

    uint16_t* acquire_next(){
        NvtxRange r("A_acquire",C_A_MGMT);

        if(prefetch_valid){
            #pragma acc wait(STREAM_PREFETCH)
            std::swap(ping,pong);
            prefetch_valid=false;
        }else{
            GenerateA_GPU_async(next_rho,n_,q_,ping,STREAM_A);
            #pragma acc wait(STREAM_A)
        }
        return ping;
    }

    void release(){
        if(!n_) return;

        NvtxRange r("DoubleBufA::release",C_A_MGMT);

        size_t An=(size_t)n_*n_;
        #pragma acc wait(STREAM_PREFETCH)
        #pragma acc exit data delete(ping[0:An],pong[0:An])
        n_=0;
    }

    ~DoubleBufA(){ release(); }
};

/* ================= MAIN ================= */
int main(int argc,char** argv)
{
    NvtxRange total_range("TOTAL_EXECUTION",C_TOTAL);

    if(argc<3){
        std::cerr<<"Usage: "<<argv[0]<<" <N> <n>\n";
        return 1;
    }

    int N=std::atoi(argv[1]);
    uint32_t n=std::atoi(argv[2]);
    uint32_t q=3329;
    uint32_t MSG=256;

    size_t total_c=(size_t)MSG*(n+1);

    std::vector<int32_t> s_buf(n),t_buf(n),z_buf(256),
                         m_buf(MSG),mp_buf(MSG),dec_buf(MSG),
                         ptxt_buf(MSG),e2_buf(MSG),
                         c_buf(total_c),c_chk(total_c);

    int32_t* sp=s_buf.data();
    int32_t* tp=t_buf.data();
    int32_t* zp=z_buf.data();
    int32_t* mp=m_buf.data();
    int32_t* mpp=mp_buf.data();
    int32_t* decp=dec_buf.data();
    int32_t* ptxtp=ptxt_buf.data();
    int32_t* e2p=e2_buf.data();
    int32_t* cp=c_buf.data();
    int32_t* cchkp=c_chk.data();

    #pragma acc enter data create(sp[0:n],tp[0:n],zp[0:256],mp[0:MSG],mpp[0:MSG],\
                                 decp[0:MSG],ptxtp[0:MSG],e2p[0:MSG],\
                                 cp[0:total_c],cchkp[0:total_c])

    DevicePool pool;
    pool.init(n,MSG);

    int32_t* rp_T=pool.r_buf_T.data();
    int32_t* e1p_T=pool.e1_buf_T.data();

    uint64_t key_seed=0;
    RAND_bytes((unsigned char*)&key_seed,sizeof(key_seed));
    uint64_t rho_seed=key_seed^0xCAFEBABE12345678ULL;

    DoubleBufA dba;
    dba.init(n,q,rho_seed);

    long long sum_keygen_us=0,sum_encaps_us=0,sum_decaps_us=0;
    int mismatches=0;

    auto startTot=std::chrono::steady_clock::now();

    for(int k=0;k<N;++k){
        NvtxRange iter("Iteration",C_TOTAL);

        RAND_bytes((unsigned char*)&key_seed,sizeof(key_seed));
        rho_seed=key_seed^0xCAFEBABE12345678ULL;

        uint16_t* Ap=dba.acquire_next();

        /* KeyGen */
        {
            NvtxRange r("KeyGen",C_KEYGEN);
            auto t0=std::chrono::steady_clock::now();

            KeyGen_GPU_rngongpu_Aflat(key_seed,n,q,Ap,sp,tp,pool);

            for(int i=0;i<256;++i) zp[i]=getRandomInt(0,1);
            #pragma acc update device(zp[0:256])

            auto t1=std::chrono::steady_clock::now();
            sum_keygen_us+=std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count();
        }

        /* Prefetch */
        if(k+1<N){
            NvtxRange r("Prefetch_overlap",C_PREFETCH);
            uint64_t rho_next=0;
            RAND_bytes((unsigned char*)&rho_next,sizeof(rho_next));
            dba.prefetch_next(rho_next);
        }

        std::vector<int32_t> K_enc;
        std::vector<int32_t> K_dec;

        /* Encaps */
        {
            NvtxRange r("Encaps",C_ENCAPS);

            auto t0=std::chrono::steady_clock::now();
            Encaps_GPU_Aflat(key_seed,n,q,Ap,tp,cp,K_enc,
                             mp,rp_T,e1p_T,e2p,ptxtp,pool);
            auto t1=std::chrono::steady_clock::now();

            sum_encaps_us+=std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count();
        }

        /* Decaps */
        {
            NvtxRange r("Decaps",C_DECAPS);

            auto t0=std::chrono::steady_clock::now();
            Decaps_GPU_Aflat(key_seed,n,q,Ap,tp,sp,zp,cp,K_dec,
                             mpp,decp,rp_T,e1p_T,e2p,ptxtp,cchkp,pool);
            auto t1=std::chrono::steady_clock::now();

            sum_decaps_us+=std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count();

            bool same=(K_enc.size()==K_dec.size());
            for(size_t i=0;i<K_enc.size() && same;i++)
                if(K_enc[i]!=K_dec[i]) same=false;

            if(!same) mismatches++;
        }
    }

    pool.release();
    dba.release();

    #pragma acc exit data delete(sp[0:n],tp[0:n],zp[0:256],mp[0:MSG],mpp[0:MSG],\
                                decp[0:MSG],ptxtp[0:MSG],e2p[0:MSG],\
                                cp[0:total_c],cchkp[0:total_c])

    auto endTot=std::chrono::steady_clock::now();
    double total_s=std::chrono::duration_cast<std::chrono::microseconds>(endTot-startTot).count()/1e6;

    std::cout<<"mismatches="<<mismatches<<"\n";
    std::cout<<n<<";"<<(double)sum_keygen_us/N<<";"
             <<(double)sum_encaps_us/N<<";"
             <<(double)sum_decaps_us/N<<";"
             <<total_s<<"\n";

    return 0;
}
