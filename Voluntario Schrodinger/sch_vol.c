#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "complex.h"
#include "gsl_rng.h"

gsl_rng *tau;

#define PI 3.1415926535
#define seed 176895

#define n 1000
#define nD 200


#define N 2000
#define lambda 1



int main(void)
{
    FILE *fout;
    fout = fopen("aux.dat", "a");

    extern gsl_rng *tau;
    tau = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(tau, seed);

    fcomplex aux1, aux2, aux3, aux4, phi[N+1], xhi[N+1], alpha[N], betta[N], gamma[N], b[N+1];

    bool bucle;
    int i, j, k, mT, nciclos;
    double ko, s, modulo, fase, v[N+1], modulophi[N+1], norma, norma_ini, Pder, Pizq, x, T;


    //Parametros
    nciclos = N/10;
    ko = (2.0*PI*nciclos)/N;
    s = 1/(ko*ko*4.0);

    //Potencial V
    for(j=0; j<N+1; j++)
    {
        if((j<2.0*N/5)||(j>3.0*N/5))
        {
            v[j] = 0;
        }
        else
        {
            v[j] = lambda*ko*ko;
        }
    }

    //Contorno
    phi[0] = Complex(0.0, 0.0);
    phi[N] = Complex(0.0, 0.0);
    xhi[0] = Complex(0.0, 0.0);
    xhi[N] = Complex(0.0, 0.0);

    alpha[N-1]=Complex(0.0, 0.0);
    betta[N-1]=Complex(0.0, 0.0);

    //Phi inicial
    for (j=1; j<N; j++)
    {
        modulo = exp( (-8.0*pow( (4*j-N),2) ) / (N*N) );
        fase = j*ko;

        phi[j] = Cgauss(fase, modulo);
    }

    //Norma
    norma = 0;
    for(j=0; j<=N; j++)
    {
        modulophi[j] = pow(Cabs(phi[j]), 2);

        norma += modulophi[j];
    }
    norma_ini = norma;

    //Normalizacion
    for(j=0; j<=N; j++)
    {
        phi[j] = RCmul(1.0/sqrt(norma_ini), phi[j]);
        modulophi[j] = pow(Cabs(phi[j]), 2);
    }

    //Gammas y Alphas
    for(j=N-2; j>=0; j--)
    {
            aux1 = Complex(-2-v[j+1], 2/s);
            aux2 = Cadd(aux1, alpha[j+1]);
            aux3 = Complex(1.0, 0.0);

            gamma[j] = Cdiv(aux3, aux2);
            alpha[j] = RCmul(-1.0, gamma[j]);
    }


    //Evolucion
    aux4 = Complex(0.0, 4/s);
    mT = 0;
    for(k=0; k<n; k++)
    {
        bucle = true;
        while(bucle)
        {
            for(i=1; i<=nD; i++)
            {

                for(j=0; j<=N; j++) //Ec. 19
                {
                    b[j] = Cmul(aux4, phi[j]); 
                }
                for(j=N-2; j>=0; j--) //Ec. 22
                {
                    aux1 = Csub(b[j+1], betta[j+1]);
                    betta[j] = Cmul(gamma[j], aux1);
                }
                for(j=0; j<N; j++) //Ec. 20
                {
                    aux2 = Cmul(alpha[j], xhi[j]);
                    xhi[j+1] = Cadd(aux2, betta[j]); 
                }
                for(j=0;j<=N;j++) //Phis
                {
                    phi[j] = Csub(xhi[j], phi[j]);
                }
            }

            Pder = 0;
            Pizq = 0;

            for(j=4*N/5; j<N; j++) //Inicio deteccion a la derecha
            {
                modulophi[j] = pow(Cabs(phi[j]), 2);
                Pder += modulophi[j];
            }
            x = gsl_rng_uniform(tau);

            if(x < Pder)
            {
                mT++;

                //Phi inicial normalizada
                phi[0] = Complex(0.0,0.0);
                phi[N] = Complex(0.0,0.0);

                for(j=1; j<N; j++)
                {
                    modulo=exp((-8.0*pow((4*j-N),2))/(N*N));
                    fase=j*ko;
                    
                    phi[j] = Cgauss(fase, modulo);
                    phi[j] = RCmul(1.0/sqrt(norma_ini), phi[j]);
                }

                bucle = false;
                break; //Ha habido deteccion a la derecha
            }


            for(j=4*N/5; j<=N; j++) // No puede estar en la derecha
            {
                phi[j] = Complex(0.0, 0.0);
            }
            norma = 0;
            for(j=0; j<=N; j++) //Nueva norma
            {
                norma += Cabs(phi[j])*Cabs(phi[j]);
            }
            for(j=0; j<=N; j++) //Nueva normalizacion
            {
                phi[j] = RCmul(1.0/sqrt(norma), phi[j]);
            }


            for(j=0; j<=N/5 ;j++) //Inicio deteccion a la izquierda
            {
                modulophi[j] = pow(Cabs(phi[j]), 2);
                Pizq += modulophi[j];
            }
            x = gsl_rng_uniform(tau);

            if(x < Pizq)
            {
                //Phi inicial normalizada
                phi[0] = Complex(0.0,0.0);
                phi[N] = Complex(0.0,0.0);

                for(j=1; j<N; j++)
                {
                    modulo=exp((-8.0*pow((4*j-N),2))/(N*N));
                    fase=j*ko;

                    phi[j] = Cgauss(fase, modulo);
                    phi[j] = RCmul(1.0/sqrt(norma_ini), phi[j]);
                }

                bucle = false;
                break; //Ha habido deteccion a la izquierda
            }

            for(j=0; j<=N/5; j++) //No puede estar a la izquierda
            {
                phi[j] = Complex(0.0, 0.0);
            }
            norma = 0;
            for(j=0; j<=N; j++) //Nueva norma
            {
                norma += Cabs(phi[j])*Cabs(phi[j]);
            }
            for(j=0; j<=N; j++) //Nueva normalizacion
            {
                phi[j] = RCmul(1.0/sqrt(norma), phi[j]);
            }
        }
    }
    
    T = mT*1.0/n;
    fprintf(fout, "%lf\n", T);

    fclose(fout);

    return 0;
}
