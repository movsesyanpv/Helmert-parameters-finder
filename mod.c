#include <stdio.h>
#include <stdlib.h>
#include <lapacke.h>
#include <math.h>
#include "proj_api.h"
//#include <mkl.h>

#define numpoints 5

void printCov(double cov[7][7], double w[7])
{
	for(int i = 0; i<7; i++)
	{
			for(int j = 0; j<7; j++)
					{
							printf("%13.3e", cov[j][i]);
					}
/*            printf(" xs=%15lf x=%15lf\n",xs[i], x[i]);*/
					printf(" w=%lf\n",w[i]);
	}
	printf("\n");
}

int main()
{
	FILE *inb, *ina;
    int i=0, j=0, l=0;
    double sum=0;
    double cov[7][7];
	double s = 0;
    //double **a;     // матрица системы
	double a[7][numpoints*3], ba[7][numpoints*3];
    //double m[7];    // столбец неизвестных
/*    double *x;      // столбец свободных членов*/
    double w[7];
    double x[numpoints*3], y[numpoints*3];
/*    double *xs;     // X_s, Y_s, Z_s*/
    double xs[numpoints*3];
    lapack_int n=numpoints;        // количество точек с известными координатами
    lapack_int k;          // количество уравнений в системе 3*n

    lapack_int info=5;
    lapack_int rank=0;

    k = 3 * n;

		inb = fopen("xb.dat", "r");
    for(i = 0; i<n; i++)
    {
        fscanf(inb,"%15lf %15lf %15lf", &xs[3*i], &xs[3*i+1], &xs[3*i+2]);
        //xs[i] = i+1;
    }
		fclose(inb);
	//xs[6] = 0;
	//xs[2] = -10;
		ina = fopen("xa.dat", "r");
    for(i = 0; i<n; i++)
    {
        fscanf(ina,"%15lf %15lf %15lf", &x[3*i], &x[3*i+1], &x[3*i+2]);
        //x[i] = i+1;
				y[3*i] = x[3*i];
				y[3*i+1] = x[3*i+1];
				y[3*i+2] = x[3*i+2];
    }
    fclose(ina);

    // заполнение матрицы системы
    for(i=0; i<n; i++)
    {
        a[0][3*i] = xs[3*i];
        a[1][3*i] = 0;
        a[2][3*i] = -xs[3*i+2];
        a[3][3*i] = xs[3*i+1];
        a[4][3*i] = 1;
				a[5][3*i] = 0;
				a[6][3*i] = 0;

				a[0][3*i+1] = xs[3*i+1];
        a[1][3*i+1] = xs[3*i+2];
        a[2][3*i+1] = 0;
        a[3][3*i+1] = -xs[3*i];
        a[4][3*i+1] = 0;
				a[5][3*i+1] = 1;
				a[6][3*i+1] = 0;

				a[0][3*i+2] = xs[3*i+2];
        a[1][3*i+2] = -xs[3*i+1];
        a[2][3*i+2] = xs[3*i];
        a[3][3*i+2] = 0;
        a[4][3*i+2] = 0;
				a[5][3*i+2] = 0;
				a[6][3*i+2] = 1;
    }

        for(i = 0; i<k; i++)
        {
            for(j = 0; j<7; j++)
                {
                    //printf("%16lf", a[j][i]);
										ba[j][i] = a[j][i];
                }
            //printf(" xs=%15lf x=%15lf\n",xs[i], x[i]);
        }
        //printf("\n\n");

/*      info = LAPACKE_dgels(LAPACK_COL_MAJOR, 'N', 9, 7, 1, *a, 9, x, 9);*/
		if((info=LAPACKE_dgelss(LAPACK_COL_MAJOR, k, 7, 1, *a, k, x, k, w, -1, &rank))==0)
		{

		// 	for(i = 0; i<k; i++)
		// 		{
		// 				for(j = 0; j<7; j++)
		// 						{
		// 								printf("%16lf", a[j][i]);
		// 						}
		// /*            printf(" xs=%15lf x=%15lf\n",xs[i], x[i]);*/
		// 						printf("\n");
		// 		}
		// 		printf("\n");

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

		//printCov(cov,w);

		//printf("%10lf\n\n", s);

		for (i = 0; i < k; i++)
		{
			sum = 0;
			for (j = 0; j < 7; j++)
			{
				sum += ba[j][i] * x[j];
			}
			s += pow((y[i] - sum), 2);
			//printf("%15lf, %10lf, %10lf\n", s, y[i], sum);

		}
		//printf("\n%10lf\n\n", s);

		s = sqrt(s / (1.0*k - 7.0));

		//printf("%10lf\n\n", s);

		for (i = 0; i < 7; i++)
		{
			cov[i][i] = sqrt(cov[i][i])*s;
		}

		//printCov(cov,w);

		printf("m =%18.10e +- %.10e\nRx=%18lf +- %.10e\nRy=%18lf +- %.10e\nRz=%18lf +- %.10e\ndX=%18lf +- %.10e\ndY=%18lf +- %.10e\ndZ=%18lf +- %.10e\n",
		x[0] - 1, cov[0][0], -x[1] / x[0]*206265.0, cov[1][1]*206265.0, -x[2] / x[0]*206265.0, cov[2][2]*206265.0, -x[3] / x[0]*206265.0,cov[3][3]*206265.0, x[4],cov[4][4], x[5],cov[5][5], x[6],cov[6][6]);

		}
		else
			printf("%3i\n",info);
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
