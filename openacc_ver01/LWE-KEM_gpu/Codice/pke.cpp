#include <vector>
#include <cstdint>
#include <iostream>

#include <omp.h> // openMP

#include "pke.h"
#include "utils.h"

/////////////////////////////////////   KeyGen  /////////////////////////////////////

void KeyGen(uint32_t n, uint32_t q, std::vector<std::vector<int32_t>> &A, std::vector<int32_t> &s_k, std::vector<int32_t> &t)
{
    A = GenerateRandomMatrixInt32(n, q - 1);
    std::vector<int32_t> Aflat = flatten_matrix(A);
    std::vector<int32_t> s = sample_vector_binomial(n);
    std::vector<int32_t> e = GenerateGaussianVector(n);
    std::vector<int32_t> prod(n, q);

    if(Aflat.size()!= n*n || s.size()!= n || e.size() != n)
    {
        std::cerr << "KeyGen: Vector size mismatch! Aflat=" << Aflat.size() << " s=" << s.size() << " e=" << e.size() << std::endl;
        return; 
    }

    std::vector<int32_t> z(256); //256 bits

    //#pragma acc kernels
    //#pragma omp parallel for //OKKKS
    for (int i = 0; i < 256; ++i)
        z[i] = getRandomInt(0, 1);

    // secret key s_k data dalla concat di s e z
    s_k = concat(s, z);
    
    // auto start_prod = std::chrono::steady_clock::now();
        
    //#pragma acc data copyin(Aflat[0:n*n], s[0:n]) copyout(prod[0:n]) ///////////
    //{
    #pragma acc enter data copyin(Aflat.data()[0:n*n], s.data()[0:n], e.data()[0:n])
    #pragma acc enter data create(prod.data()[0:n])
        #pragma acc kernels
        for (uint32_t i = 0; i < n; ++i)
        {
            long long acc = 0;

            //#pragma acc loop reduction(+:acc)
            for (uint32_t j = 0; j < n; ++j)
            {
                acc += (long long)Aflat.data()[i*n + j] * (long long)s.data()[j];
            }

            long long r = acc % (long long)q;
            if (r < 0) r += (long long)q;
            prod.data()[i] = (int32_t)r;
            
        }
    //}
    // auto end_prod = std::chrono::steady_clock::now();
    // auto elapsed_prod = std::chrono::duration_cast<std::chrono::microseconds>(end_prod - start_prod);
    // std::cout << "Prod time: " << elapsed_prod.count() << " mus\n";

    std::vector<int32_t> t1(n, 0);
    //#pragma acc data copyin(prod[0:n], e[0:n]) copyout(t1[0:n]) ///////////
    //{
        #pragma acc parallel loop
        for (uint32_t i = 0; i < n; ++i)
        {
            long long v = (long long)prod.data()[i] + (long long)e.data()[i];
            long long r = v % (long long)q;
            if (r < 0) r += (long long)q;
            t1[i] = (int32_t)r;
        }
        #pragma acc exit data copyout(prod.data()[0:n])
        #pragma acc exit data delete(Aflat.data()[0:n*n], s.data()[0:n], e.data()[0:n])
    //}
    t = t1;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////   Encrypt /////////////////////////////////////
void Encrypt(uint32_t n, uint32_t q, std::vector<int32_t> &t, std::vector<int32_t> &u, int32_t &v_i, uint32_t plaintext_i, std::vector<int32_t> &r, std::vector<int32_t> &e1, int32_t &e2, const std::vector<std::vector<int32_t>> &AT)
{
    std::vector<int32_t> ATflat = flatten_matrix(AT);
    std::vector<int32_t> prod(n, q);

    u.resize(n);
    //#pragma acc data copyin(ATflat[0:n*n], r[0:n], e1[0:n]) copyout(u[0:n]) /////////
    //{
    #pragma acc enter data copyin(ATflat.data()[0:n*n], r.data()[0:n], e1.data()[0:n])
    #pragma acc enter data create(u.data()[0:n])
        #pragma acc kernels
        for (uint32_t i = 0; i < n; ++i)
        {
            long long acc = 0;

            //#pragma acc loop reduction(+:acc)
            for (uint32_t j = 0; j < n; ++j)
            {
                acc += (long long)ATflat.data()[i*n + j] * (long long)r.data()[j];
            }

            long long val = acc + (long long)e1.data()[i];
            long long m = val % (long long)q;
            if (m < 0) m += (long long)q;
            u.data()[i] = (int32_t)m;
            #pragma acc exit data copyout(u.data()[0:n])
            #pragma acc exit data delete(ATflat.data()[0:n*n], r.data()[0:n], e1.data()[0:n])
        }
    //}

    long long dot = 0;
    #pragma acc enter data copyin(t.data()[0:n], r.data()[0:n]) ///////////////
        //#pragma acc parallel loop reduction(+:dot)
        #pragma acc kernels loop reduction(+:dot)
        for (uint32_t i = 0; i < n; ++i)
            dot += (long long)t.data()[i] * (long long)r.data()[i];
    #pragma acc exit data delete(t.data()[0:n], r.data()[0:n])

    //#pragma acc update self(dot)

    long long res = dot % (long long)q;
    if (res < 0) res += (long long)q;

    long long vv = res + (long long)e2 + (long long)plaintext_i;
    vv %= (long long)q;
    if (vv < 0) vv += (long long)q;
    v_i = (int32_t)vv;
}
/////////////////////////////////////   Decrypt /////////////////////////////////////

void Decrypt(int32_t v_i, const std::vector<int32_t> &u, const std::vector<int32_t> &s_k, uint32_t q, int32_t &decrypt_i)
{
    long long dot = 0;

    // Considero solo la prima parte di s_k, quella relativa ad s
    const size_t n = u.size(); // la dimensione di s (senza z) è n = dimensione di u
    std::vector<int32_t> s(s_k.begin(), s_k.begin() + n); // definisco s uguale al primo elemento di s_k fino all'elemento n-esimo di s_k

    //#pragma omp parallel for reduction(+:dot)  //OKKKS
    
    #pragma acc enter data copyin(s.data()[0:n], u.data()[0:n])
        #pragma acc loop seq 
        for (size_t i = 0; i < s.size(); ++i)
            dot += (long long)s.data()[i] * (long long)u.data()[i];
    //#pragma acc exit data delete(s.data()[0:n], u.data()[0:n])
    
    long long r = ((long long)v_i - dot) % (long long)q;
    if (r < 0)
        r += (long long)q;
    int32_t mu = (int32_t)r;

    const int32_t bound = (int32_t)q / 4;
    decrypt_i = (mu <= bound || mu >= (int32_t)q - bound) ? 0 : (int32_t)q / 2;
}