#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


void Evolucion(double x[], double y[], double vx[], double vy[], double masas[], double ax[], double ay[], double aux_x[], double aux_y[], int cuerpos, double h, double h2, int saltos);
void Aceleracion(double x[], double y[], double ax[], double ay[], double masas[], int cuerpos);


int main (void)
{
    FILE *finiciales, *fposiciones, *fposicionestierra, *fperiodos;


    char planetas[8][10] = {"Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"};
    double Msol, medida, t, h, h2;
    int i, k, l, cuerpos, saltos, vueltas[cuerpos];
         
    double masas[cuerpos], x[cuerpos], y[cuerpos], vx[cuerpos], vy[cuerpos], ax[cuerpos], ay[cuerpos], aux_x[cuerpos], aux_y[cuerpos];
    double posicionx[cuerpos], posiciony[cuerpos], anteriorx[cuerpos], anteriory[cuerpos], xt[cuerpos], yt[cuerpos];

    t = 0;
    h = 0.001;
    h2 = h/2;
    medida = 0.1;
    saltos = round(medida/h);
    cuerpos = 9;
    Msol = 1988500;
   

    //Condiciones iniciales    
    finiciales = fopen("ci2.dat","r");
    fperiodos = fopen("periodos2.txt", "w");

    for(i=0; i<cuerpos; i++)
    {
        fscanf(finiciales,"%lf, %lf, %lf, %lf, %lf", &x[i], &y[i], &vx[i], &vy[i], &masas[i]);
       
        //Posiciones auxiliares
        posicionx[i] = x[i];
        posiciony[i] = y[i];
        anteriorx[i] = x[i];
        anteriory[i] = y[i];

        //Reescalamiento  
        vx[i] = vx[i]*58.1;
        vy[i] = vy[i]*58.1;
        masas[i] /= Msol; 
    }  
    fclose(finiciales);

    for(l=0; l<cuerpos; l++)
    {
        vueltas[l] = 0;
    }
    


    fposiciones = fopen("pos2.txt","w");
    fposicionestierra = fopen("posTierra2.txt", "w");

    //Evolucion
    for(k=0; k<2500; k++)
    {
        Evolucion(x, y, vx, vy, masas, ax, ay, aux_x, aux_y, cuerpos, h, h2, saltos);    

        for(l=0; l<cuerpos; l++) //Nuevas posiciones
        {
            fprintf(fposiciones, "%lf, %lf\n", x[l], y[l]);
        }
        fprintf(fposiciones, "\n");

        t += saltos*h;  


        for (l=1; l<cuerpos; l++) //Periodos
        {
            if(anteriory[l]<posiciony[l] && y[l]>posiciony[l] && posicionx[l]>0 && vueltas[l]<1)
            {
                fprintf(fperiodos, "%s: %lf years\n", planetas[l-1], t*58.1/365.256);
                vueltas[l]++;
            }
            if(anteriory[l]>posiciony[l] && y[l]<posiciony[l] && posicionx[l]<0 && vueltas[l]<1)
            {
                fprintf(fperiodos, "%s: %lf years\n", planetas[l-1], t*58.1/365.256);
                vueltas[l]++;
            }
        }
        
        for(l=1; l<cuerpos; l++)
        {
            anteriory[l] = y[l];
        }
        

        for(l=0; l<cuerpos; l++) //Coordenadas Tierra
        {
            xt[l] = x[l]-x[3];
            yt[l] = y[l]-y[3];

            fprintf(fposicionestierra, "%lf, %lf\n", xt[l], yt[l]);
        }
        fprintf(fposicionestierra, "\n");  
    }
   fclose(fposiciones);
   fclose(fposicionestierra);
   fclose(fperiodos);

    return 0;
}



void Evolucion(double x[], double y[], double vx[], double vy[], double masas[], double ax[], double ay[], double aux_x[], double aux_y[], int cuerpos, double h, double h2, int saltos)
    {
        int i, j;
        Aceleracion(x, y, ax, ay, masas, cuerpos);

        for(i=0; i<saltos; i++)
        {
            for(j=0; j<cuerpos; j++)
            {
                aux_x[j] = vx[j] + h2*ax[j];
                aux_y[j] = vy[j] + h2*ay[j];

                x[j] = x[j] + h*aux_x[j];
                y[j] = y[j] + h*aux_y[j];
            }

            Aceleracion(x, y, ax, ay, masas, cuerpos);

            for(j=0;j<cuerpos; j++)
            {
                vx[j] = aux_x[j] + h2*ax[j];
                vy[j] = aux_y[j] + h2*ay[j];
            }
        }
        return;
    }



void Aceleracion(double x[], double y[], double ax[], double ay[], double masas[], int cuerpos)
{
    int i, j;
    double fuerzax, fuerzay, fuerzatotal;
    double dx, dy, d, modulo;

    for(i=0; i<cuerpos; i++)
    {
        ax[i] = 0;
        ay[i] = 0;
    }

    for(i=0; i<cuerpos; i++)
    {
        for(j=i+1; j<cuerpos; j++)
        {
            dx = x[i]-x[j];
            dy = y[i]-y[j];
            d = dx*dx + dy*dy;
            modulo = pow(d, 3./2.);

            fuerzatotal = masas[i]*masas[j]/modulo;
            fuerzax = (-1)*fuerzatotal*dx;
            fuerzay = (-1)*fuerzatotal*dy;

            ax[i] += fuerzax/masas[i];
            ay[i] += fuerzay/masas[i];
            ax[j] -= fuerzax/masas[j];
            ay[j] -= fuerzay/masas[j];
        }
    }
    return;
}
