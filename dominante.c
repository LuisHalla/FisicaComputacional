#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Metodo de las potencias CON ESCALADO: 
dada una matriz 2x2 simetrica "M" y un vector semilla "X0", se itera la operacion Xn = MXn-1, diviendo en cada    paso el vector Xn por su maxima componente (factor de escala). Dicho factor converge al valor propio dominante de la matriz M.*/

int main(void)
{
    /*Declaracion de variables*/
    int i, n;
    double m[3], x0[2], xn[2], vp, maximo, dif;
    FILE *f;
    
    f = fopen("vp_datos.txt", "w"); /*Abrir fichero*/

    /*Presentacion del programa*/
    printf("\nMatriz 2x2 simetrica\t(a  b)\n");
    printf("                      \t(b  c)\n");

    /*Peticion de valalores por teclado*/
    printf("\nInserte los coeficientes\n");

    printf("a: ");
    scanf("%lf", &m[0]);
    printf("b: ");
    scanf("%lf", &m[1]);
    printf("c: ");
    scanf("%lf", &m[2]);

    printf("\nInserte las componentes del vector semilla\n");

    printf("x0[0]: ");
    scanf("%lf", &x0[0]);
    printf("x0[1]: ");
    scanf("%lf", &x0[1]);

    printf("\nNumero de iteraciones (n>=3)\n");
    scanf("%i", &n);

    /*Calculo exacto del valor propio dominante*/
    vp = m[2] + 0.5*( (m[0]-m[2]) + sqrt( (m[2]-m[0])*(m[2]-m[0]) + 4*m[1]*m[1] ) );
    
    /*Ejecucion del metodo de las potencias y escritura en fichero de cada iteracion para posterior representacion*/
    for(i=1; i<=n; i++)
    {
        xn[0] = m[0]*x0[0] + m[1]*x0[1];
        xn[1] = m[1]*x0[0] + m[2]*x0[1];

        maximo = fmax(xn[0],xn[1]);
   
        x0[0] = xn[0] / maximo;
        x0[1] = xn[1] / maximo;

        dif = maximo - vp;
        dif = fabs(dif);
        
        printf("\n%lf\n", dif);
        
        /*Se suprime la primera iteracion desde aqui porque no se como hacerlo desde matplotlib. De no hacerlo, las graficas se visualizan muy mal*/
        if(i!=1)
        {
            fprintf(f,"%i %lf %lf\n", i, maximo, dif);
        }
    }

    /*Resultados por pantalla*/
    printf("\nResultados:\n");
    printf(" vp exacto = %lf\n", vp);
    printf(" vp power method = %lf\n", maximo);
    printf("\n Diferencia = %lf\n\n", dif);

    fclose(f); /*Cerrar fichero*/

    return 0;
}
