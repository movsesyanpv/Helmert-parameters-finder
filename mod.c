#include <stdio.h>
#include <stdlib.h>
#include <lapacke.h>
#include <math.h>

int main()
{
    FILE *in;
    int i=0, j=0, l=0;
    double sum=0;
    double cov[7][7];
    //double **a;     // матрица системы
    double a[7][9];
    //double m[7];    // столбец неизвестных
/*    double *x;      // столбец свободных членов*/
    double w[7];
    double x[9];
/*    double *xs;     // X_s, Y_s, Z_s*/
    double xs[9];
    lapack_int n=3;        // количество точек с известными координатами
    lapack_int k;          // количество уравнений в системе 3*n
    
    lapack_int info=5;
    lapack_int rank=0;
    
    k = 3 * n;
    
/*    a = (double**) malloc(7*sizeof(double*));*/
/*    for(int i=0; i<7; i++)*/
/*    {*/
/*        a[i] = (double*) malloc(k * sizeof(double));  */
/*    }*/
    
/*    x = (double*) malloc(k * sizeof(double));*/
/*    xs = (double*) malloc(k * sizeof(double));*/
    
    in = fopen("xxs.dat", "r");
    for(i = 0; i<k; i++)
    {
        fscanf(in,"%15lf", &xs[i]);
/*        xs[i] = i+1;*/
    }
    
    for(i = 0; i<k; i++)
    {
        fscanf(in,"%15lf", &x[i]);
/*        x[i] = i+1;*/
    }
    fclose(in);
    
    // заполнение матрицы системы
    for(i=0; i<n; i++)
    {
        a[0][3*i] = xs[3*i];
        a[1][3*i] = -xs[3*i+1];
        a[2][3*i] = xs[3*i+2];
        a[3][3*i] = 1;
        a[4][3*i] = a[5][3*i] = a[6][3*i] = 0;
        
        a[0][3*i+1] = xs[3*i+1];
        a[1][3*i+1] = xs[3*i];
        a[4][3*i+1] = -xs[3*i+2];
        a[5][3*i+1] = 1;
        a[2][3*i+1] = a[3][3*i+1] = a[6][3*i] = 0;
        
        a[0][3*i+2] = xs[3*i+2];
        a[2][3*i+2] = -xs[3*i];
        a[4][3*i+2] = xs[3*i+1];
        a[6][3*i+2] = 1;
        a[1][3*i+2] = a[3][3*i+2] = a[5][3*i+2] = 0;
    }
    
        for(i = 0; i<k; i++)
        {   
            for(j = 0; j<7; j++)
                { 
                    printf("%16lf", a[j][i]);
                }
            printf(" xs=%15lf x=%15lf\n",xs[i], x[i]);
        }
        printf("\n\n");
    
/*      info = LAPACKE_dgels(LAPACK_COL_MAJOR, 'N', 9, 7, 1, *a, 9, x, 9);*/
    info = LAPACKE_dgelss(LAPACK_COL_MAJOR, 9, 7, 1, *a, 9, x, 9, w, -1, &rank);
      
      for(i = 0; i<k; i++)
        {   
            for(j = 0; j<7; j++)
                { 
                    printf("%16lf", a[j][i]);
                }
/*            printf(" xs=%15lf x=%15lf\n",xs[i], x[i]);*/
                printf("\n");
        }
        printf("\n");
    
    printf("m =%.10e\nRx=%lf\nRy=%lf\nRz=%lf\ndX=%lf\ndY=%lf\ndZ=%lf\n", x[0]-1, x[4]/x[0], x[2]/x[0], x[1]/x[0], x[3], x[5], x[6]);
    
        for(j = 0; j<7; j++)
        {
            for(l = 0; l<7; l++)
            {
                sum = 0;   
                for(i = 0; i<7; i++)
                {
                    if(w[i]>1.0E-8)
                    {
                        sum += a[i][j]*a[i][k]/w[i]/w[i];
                    }
                }
                cov[j][l] = sum;
            }
        }
        
        for(i = 0; i<7; i++)
        {   
            for(j = 0; j<7; j++)
                { 
                    printf("%13.3e", cov[j][i]);
                }
/*            printf(" xs=%15lf x=%15lf\n",xs[i], x[i]);*/
                printf(" w=%lf\n",w[i]);
        }
        printf("\n");
    
/*    printf("ready to free\n");*/
    
/*    for(int i=0; i<7; i++)*/
/*    {*/
/*        free(a[i]);*/
/*        printf("a[i]");*/
/*    }*/
/*    free(a);*/
/*    printf("a");*/
    
/*    free(x);*/
/*    printf("x");*/
/*    free(xs);*/
}
