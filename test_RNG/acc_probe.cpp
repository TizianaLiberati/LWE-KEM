#include <cstdio>
#include <cstdlib>

extern "C" void rngongpu_fill_normal(double* d_out, int n, double sigma);

int main()
{
    const int n = 16;
    double* h = (double*)std::malloc(n * sizeof(double));
    //if (!h) return 1;

    //#pragma acc data create(h[0:n])
    //{
        // RNGonGPU riempie il buffer sulla GPU
       // #pragma acc host_data use_device(h)
       // {
         //   rngongpu_fill_normal(h, n, 1.0);
        //}

        // Copio su CPU solo per stampare i valori generati
       // #pragma acc update self(h[0:n])

        std::printf("Valori generati da RNGonGPU:\n");
       // for (int i = 0; i < 8; ++i) {
        //    std::printf("h[%d] = %.8f\n", i, h[i]);
       // }

        // Ora OpenACC usa quegli stessi valori sulla GPU
        //#pragma acc parallel loop present(h[0:n])
        //for (int i = 0; i < n; ++i) {
          //  h[i] = 2.0 * h[i];
       // }

        //#pragma acc update self(h[0:n])

        //std::printf("\nDopo il loop OpenACC (dovrebbero essere raddoppiati):\n");
        //for (int i = 0; i < 8; ++i) {
          //  std::printf("h[%d] = %.8f\n", i, h[i]);
       // }
   // }

    std::free(h);
    return 0;
}
