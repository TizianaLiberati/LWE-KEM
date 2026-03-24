#include <stdio.h>
#include <stdlib.h>

//extern "C" void saxpy(int n, float a, float *x, float *y);
extern "C" void rngongpu_fill_normal(double* d_out, int n, double sigma);

//OpenACC chiama CUDA (senza libreria per la generazione dei numeri crittograficamente sicuri)
/* int main()
{
    int n = 1 << 20;
    //Creo vettori x e y su CPU
    float *x = (float*) malloc(n * sizeof(float));
    float *y = (float*) malloc(n * sizeof(float));

    //Creo copia su GPU di x e y con create e poi copia con copyout su OpenACC
    #pragma acc data create(x[0:n]) copyout(y[0:n])
    {
        //OpenACC esegue ciclo su x e y su GPU
        #pragma acc parallel loop
        for (int i = 0; i < n; i++) {
            x[i] = 1.0f;
            y[i] = 0.0f;
        }

        //OpenACC chiama CUDA 
        #pragma acc host_data use_device(x,y)
        {
            //saxpy(n, a, *x, *y) = y[i] += a * x[i] 
            //-> y[i] = y[i] + a* x[i] = 0 + 2 * 1 = 2
            saxpy(n, 2.0f, x, y); //y[i] = 2
        }
    }

    printf("y[0]   = %f\n", y[0]);
    printf("y[5]   = %f\n", y[5]);
    printf("y[n-1] = %f\n", y[n-1]);

    free(x);
    free(y);
    return 0;
} */

int main()
{
    // riempie a "coppie", quindi metti dimensione array sempre pari (al massimo scarti un elemento se non ti serve)
    int n = 4;
    double* h = (double*) malloc(n * sizeof(double));

    //tolgo cpyout perchè già fa le copie esplicit con self
    #pragma acc data create(h[0:n]) copyout(h[0:n])
    {
        //chiama funzione libreria
        #pragma acc host_data use_device(h)
        {
            rngongpu_fill_normal(h, n, 1.0);
        }
        /*
        // Copia i valori generati dalla GPU alla CPU per poterli stampare
        //self copia GPU -> CPU
        #pragma acc update self(h[0:n])

        printf("Valori dopo RNGonGPU:\n");
        for (int i = 0; i < n; ++i) {
            printf("h[%d] = %.8f\n", i, h[i]);
        }
        // Modifica dell'array sulla GPU con OpenACC
        #pragma acc parallel loop
        for (int i = 0; i < n; ++i) {
            h[i] += 2.0;
        }
        // Riporta i valori modificati dalla GPU alla CPU
        #pragma acc update self(h[0:n])

        printf("\nValori dopo modifica OpenACC (+2):\n");
        for (int i = 0; i < n; ++i) {
            printf("h[%d] = %.8f\n", i, h[i]);
        }
        */
    }

    
    for (int i = 0; i < n; ++i) {
        printf("h[%d] = %.8f\n", i, h[i]);
    }
    

    free(h);
    return 0;
}