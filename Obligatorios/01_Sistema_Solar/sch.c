#include <stdio.h>
#include <math.h>
#include "complex.h"

#define PI 3.1415926535
#define N 1000
#define nciclos 50
#define lambda 5
#define n 500


int main(void)
{
    FILE *fout, *fnorma;
    fout = fopen("sch.dat", "w");
    fnorma = fopen("norm.dat", "w");


    fcomplex phi[N+1], xhi[N+1], alpha[N], betta[N], gamma[N], b[N+1];
    fcomplex den1, den, num, cte, factor, mult;

    int i, j, k;
    double ko, s, modulo, fase, v[N+1], modulophi[N+1], norma;

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
        
        phi[j] = RCmul(1/sqrt(110.778365), phi[j]); //Escala  
    }


    for(j=0; j<=N; j++) //Modulo y escritura
    {
        modulophi[j] = pow(Cabs(phi[j]), 2);

        fprintf(fout, "%i,%lf,%lf\n", j, modulophi[j] ,v[j]);
    }
    fprintf(fout, "\n");


    //Gammas y Alphas
    for(j=N-2; j>=0; j--)
    {
            den1 = Complex(-2-v[j+1], 2/s);
            den = Cadd(den1, alpha[j+1]);
            num = Complex(1.0, 0.0);

            gamma[j] = Cdiv(num, den);
            alpha[j] = RCmul(-1.0, gamma[j]);
    }


    // Norma
    norma=0;

    for(j=0;j<=N;j++)
    {
        norma += modulophi[j];
    }
    fprintf(fnorma, "%lf\n", norma);


    //Evolucion

    cte = Complex(0.0, 4/s);

    for(k=0; k<=n; k++)
    {

        for(j=0; j<=N; j++) //Ec. 19
        {
            b[j] = Cmul(cte, phi[j]);
            
        }
        for(j=N-2; j>=0; j--) //Ec. 22
        {
            factor = Csub(b[j+1], betta[j+1]);
            betta[j] = Cmul(gamma[j], factor);
        }
        for(j=0; j<N; j++) //Ec. 20
        {

            mult = Cmul(alpha[j], xhi[j]);
            xhi[j+1] = Cadd(mult, betta[j]); 
        }


        for(j=0;j<=N;j++) //Phis
        {
            phi[j] = Csub(xhi[j], phi[j]);
            modulophi[j] = pow(Cabs(phi[j]), 2);

            fprintf(fout, "%i,%lf,%lf\n", j, modulophi[j], v[j]);
        }
        fprintf(fout, "\n");

        //Norma
        norma=0;

        for(j=0;j<=N;j++)
        {
            norma=norma+modulophi[j];
        }
        fprintf(fnorma, "%lf\n", norma);
    }


    fclose(fout);
    fclose(fnorma);

    return 0;
}
