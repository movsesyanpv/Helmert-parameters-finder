#include <stdio.h>
#include <stdlib.h>
//#include <lapacke.h>
#include <math.h>
#include <mkl.h>

int main()
{
	FILE *in;
    int i=0, j=0, l=0;
    double sum=0;
    double cov[7][7];
	double s = 0;
    //double **a;     // матрица системы
	double a[7][9], ba[7][9];
    //double m[7];    // столбец неизвестных
/*    double *x;      // столбец свободных членов*/
    double w[7];
    double x[9], y[9];
/*    double *xs;     // X_s, Y_s, Z_s*/
    double xs[9];
    lapack_int n=3;        // количество точек с известными координатами
    lapack_int k;          // количество уравнений в системе 3*n
    
    lapack_int info=5;
    lapack_int rank=0;
    
    k = 3 * n;
    
    info = fopen_s(&in,"E:\s.dat", "r");
    for(i = 0; i<k; i++)
    {
        fscanf_s(in,"%15lf", &xs[i]);
        //xs[i] = i+1;
    }
	//xs[6] = 0;
	//xs[2] = -10;
    
    for(i = 0; i<k; i++)
    {
        fscanf_s(in,"%15lf", &x[i]);
        //x[i] = i+1;
		y[i] = x[i];
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
        a[2][3*i+1] = a[3][3*i+1] = a[6][3*i+1] = 0;
        
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
					ba[j][i] = a[j][i];
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
                    printf("%20.10e", a[j][i]);
                }
/*            printf(" xs=%15lf x=%15lf\n",xs[i], x[i]);*/
                printf("\n");
        }
        printf("\n");
    
        for(j = 0; j<7; j++)
        {
            for(l = 0; l<7; l++)
            {
                sum = 0;   
                for(i = 0; i<7; i++)
                {
                    if(abs(w[i])>1.0E-16)
                    {
                        sum += a[i][j]*a[i][l]/w[i]/w[i];
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

		for (i = 0; i < 9; i++)
		{
			sum = 0;
			for (j = 0; j < 7; j++)
			{
				sum += ba[j][i] * x[j];
			}
			s += pow((y[i] - sum), 2);
		}
		s = sqrt(s / (9.0 - 7.0));

		for (i = 0; i < 7; i++)
		{
			cov[i][i] = sqrt(cov[i][i])*s;
		}

		printf("m =%18.10e +- %.10e\nRx=%18lf +- %.10e\nRy=%18lf +- %.10e\nRz=%18lf +- %.10e\ndX=%18lf +- %.10e\ndY=%18lf +- %.10e\ndZ=%18lf +- %.10e\n", x[0] - 1, cov[0][0], x[4] / x[0], cov[1][4], x[2] / x[0], cov[2][2], x[1] / x[0],cov[1][1], x[3],cov[3][3], x[5],cov[5][5], x[6],cov[6][6]);
    
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
