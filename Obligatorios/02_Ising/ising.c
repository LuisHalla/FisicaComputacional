#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "gsl_rng.h"

gsl_rng *tau;

#define seed = 7466024


int main (void)
{
    FILE *fespines;
    fespines=fopen("ising_desordenada_7466024_T=1.dat","w");

    extern gsl_rng *tau;
    tau=gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(tau, seed);

    bool ordenada;
    double x, T, E, C, p;
    int i, j, x, N, pMC, s[N][N], n, m, na, nd, ma, md;


    //Parametros
    T=1.0;

    N=50;

    pMC=400;


    //Configuracion inicial
    if(ordenada)
    {
        for(i=0; i<N; i++)
        {
            for(j=0; j<N; j++)
            {
                s[i][j]=1;    
            }
        }
    }
    else
    {
        for(i=0; i<N; i++)
        {
            for(j=0; j<N; j++)
            {
                x = gsl_rng_uniform_int(tau, 2);

                if(x == 0)
                {
                    s[i][j] = x-1;
                }
                else s[i][j] = x;         
            }
        }
    }

   
    for(i=0; i<N; i++) //Escritura en fichero
    {
        for(j=0; j<N-1; j++)
        {
            fprintf(fespines, "%i,", s[i][j]);
        }
        fprintf(fespines, "%i\n", s[i][N-1]);
    }
    fprintf(fespines, "\n");


    //Evolucion
    for(i=0; i<pMC; i++)
    {
        for(j=0; j<=N*N; j++)
        {
            n = gsl_rng_uniform_int(tau, N);
            m = gsl_rng_uniform_int(tau, N);

            //Condiciones periódicas

            //Filas
            if(n == 0)
            {
                na = N-1;
                nd = 1;
            }
            else if(n == N-1)
            {
                na = N-2;
                nd = 0;
            }
            else
            {
                na = n-1;
                nd = n+1;
            }

            //Columnas
            if(m == 0)
            {
                ma = N-1;
                md = 1;
            }
            else if(m == N-1)
            {
                ma = N-2;
                md = 0;
            }
            else
            {
                ma = m-1;
                md = m+1;
            }

            //Condicion de cambio
            E = 2 * s[n][m] * ( s[nd][m] + s[na][m] + s[n][ma] + s[n][md] );
            C = exp(-E/T);

            if(1<=C)
            {
                p=1;
            }
            else
            {
                p=C;
            }
            
            x = gsl_rng_uniform(tau);
            if(x < p)
            {
                s[n][m] = -s[n][m];
            }
        }

        for(i=0; i<N; i++) //Escritura en fichero
        {
            for(j=0; j<N-1; j++)
            {
                fprintf(fespines, "%i,", s[i][j]);
            }
            fprintf(fespines, "%i\n", s[i][N-1]);
        }
        fprintf(fespines, "\n");
    }

    fclose(fespines);

    return 0;
}
