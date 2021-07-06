#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define G 6.67E-11
#define Mt 5.9736E24
#define Ml 0.07349E24
#define dtl 3.844E8
#define w 2.6617E-6
#define Rt 6.378160E6
#define Rl 1.7374E6
#define m 4.5932E4
#define pi 3.141592

#define delta (1.0*G*Mt/(pow(dtl,3)))
#define mu (1.0*Ml/Mt)


void Kutta(double kutta[], double vector[], double tiempo, double h);



int main (void)
{
    FILE *fposiciones, *fconservacion;

    fposiciones = fopen("pos.txt","w");
    fconservacion = fopen("Hprima.txt","w");


    int i, j, l, N, n;
    double t, h, iter, r, phi, pr, pphi, theta, vescape, H, factor;
    double x[n], y[n], f[N], k1[N], k2[N], k3[N], k4[N], aux[N];

    t = 0.0;
    h = 60;
    iter = 100; 
    N = 4;
    n = 3;


    //Condiciones iniciales

    theta = M_PI/4.0; //45 grados
    vescape = 11116; //Vesc de la Tierra


    x[0] = 0;
    y[0] = 0;

    x[2] = 1;
    y[2] = 0;


    f[0] = 1.0*Rt/dtl;
    f[1] = M_PI/4.0;
    f[2] = vescape*cos(theta-f[1])/dtl;
    f[3] = f[0]*vescape*sin(theta-f[1])/dtl;

    printf("%lf, %lf, %lf, %lf\n", f[0], f[1], f[2], f[3]);
 
    x[1]=f[0]*cos(f[1]); //Coordenadas polares
    y[1]=f[0]*sin(f[1]); 

    for(l=0; l<n; l++)
        {
            fprintf(fposiciones, "%lf, %lf\n", x[l], y[l]);
        }
    fprintf(fposiciones, "\n");


 
    //Hprima
    r = (f[0])*dtl;
    phi = f[1];
    pr = (f[2])*m*dtl;
    pphi = (f[3])*m*dtl*dtl;

    factor = pow( r*r + dtl*dtl - 2*r*dtl*cos( -w*t+phi ), 1/2 );

    H = 0.5*( (pr*pr)/m + (pphi*pphi)/(m*r*r) ) - (G*Mt*m/r) - G*Ml*m/(factor) - w*pphi;

    fprintf(fconservacion, "%e\n",H);


    //Evolucion
    for(j=0; j<iter; j++)
    {
        for(i=0 ;i<100; i++)
        {
            Kutta(k1, f, t, h); //k1
        

            for(l=0; l<N ;l++) //k2
            {
               aux[l] = f[l] + (k1[l]/2.0);
            }   
            Kutta(k2, aux, t+h/2.0, h);


            for(l=0; l<N; l++) //k3
            {
               aux[l] = f[l] + k2[l]/2.0;
            }   
            Kutta(k3, aux, t+h/2.0, h);


            for(l=0; l<N; l++) //k4
            {
               aux[l] = f[l] + k3[l]/2.0;
            }   
            Kutta(k4, aux, t+h, h);



            for(l=0;l<N;l++) //Nuevas coordenadas
            {
                f[l] += ( k1[l] + 2*(k2[l] + k3[l]) + k4[l] ) / 6.0;
            }

            t += h;
        }
    
        //Nuevas coordenadas
        x[2] = cos(w*t);
        y[2] = sin(w*t);

        x[1] = f[0]*cos(f[1]);
        y[1] = f[0]*sin(f[1]);
        
        for(l=0; l<n; l++)
        {
            fprintf(fposiciones, "%lf, %lf\n", x[l], y[l]);
        }
        fprintf(fposiciones, "\n");


        //Hprima
        r = (f[0])*dtl;
        phi = f[1];
        pr = (f[2])*m*dtl;
        pphi = (f[3])*m*dtl*dtl;

        factor = pow( r*r + dtl*dtl - 2*r*dtl*cos(-w*t+phi), 1/2 );

        H = 0.5*( (pr*pr)/m + (pphi*pphi)/(m*r*r) ) - G*Mt*m/r - G*Ml*m/(factor) - w*pphi;

        fprintf(fconservacion, "%e\n",H);
    }
    fclose(fposiciones);
    fclose(fconservacion);

    return 0;
}



void Kutta (double kutta[], double vector[], double tiempo, double h)
{
    double rprima, aux2;

    rprima = sqrt( 1 + vector[0]*vector[0] - 2*vector[0]*cos(vector[1] - w*tiempo) );
    aux2 = pow(rprima, 3);

    kutta[0] = h * vector[2];

    kutta[1] = h * vector[3] / (vector[0]*vector[0]);

    kutta[2] = h * ( vector[3]*vector[3] / pow(vector[0], 3 ) - delta*( (1/(vector[0]*vector[0])) + (mu/aux2)*( vector[0] - cos(vector[1] - w*tiempo) ) ) );

    kutta[3] = h * ( ( -delta*mu*vector[0]/aux2 )*sin(vector[1] - w*tiempo) );

        return;
}
