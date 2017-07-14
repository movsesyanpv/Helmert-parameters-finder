#include <stdio.h>
#include <math.h>

int main()
{
    double dX = 23.57;
    double dY = -140.95;
    double dZ = -79.8;
    double m = -0.00000022;
    double Rx = 0;
    double Ry = -7.0/72000.0;
    double Rz = -79.0/360000.0;
    
    double lat[6];
    double lng[6];
    double x[6];
    double y[6];
    double z[6];
    double a = 6378137;
    double b = 6356752.3142;
    double rho;
    
    FILE *o;
    
    lat[0] = 60;
    lat[1] = 50.5;
    lat[2] = 55.8;
    lat[3] = 60.0000317786;
    lat[4] = 50.500060278;
    lat[5] = 558.79996037846;
    lng[0] = 30;
    lng[1] = 36.5;
    lng[2] = 37.4;
    lng[3] = 30.002262925;
    lng[4] = 36.50164478;
    lng[5] = 37.4018805332;
    
    for(int i=0; i<3; i++)
    {
        rho = pow(a,2)/(sqrt(pow(a*cos(lat[i]),2)+pow(b*sin(lat[i]),2)));
        x[i] = rho * cos(lat[i]) * cos(lng[i]);
        y[i] = rho * cos(lat[i]) * sin(lng[i]);
        z[i] = rho * pow(a/b,2) * sin(lat[i]);
    }
    
    a = 6378425;
    b = 6356863.019;
    for(int i=3; i<6; i++)
    {
/*        rho = pow(a,2)/(sqrt(pow(a*cos(lat[i]),2)+pow(b*sin(lat[i]),2)));*/
/*        x[i] = rho * cos(lat[i]) * cos(lng[i]);*/
/*        y[i] = rho * cos(lat[i]) * sin(lng[i]);*/
/*        z[i] = rho * pow(a/b,2) * sin(lat[i]);*/
        x[i] = (1+m)*(    x[i-3]-Rz*y[i-3]+Ry*z[i-3])+dX;
        y[i] = (1+m)*( Rz*x[i-3]+   y[i-3]-Rx*z[i-3])+dY;
        z[i] = (1+m)*(-Ry*x[i-3]+Rx*y[i-3]+   z[i-3])+dZ;
    }
    
    o = fopen("xxs.dat","w");
    for(int i = 0; i<6; i++)
    {
        fprintf(o,"%15lf\n%15lf\n%15lf\n", x[i], y[i], z[i]);
    }
    fclose(o);
    
    o = fopen("params.dat","w");
    fprintf(o,"m =%.10e\nRx=%lf\nRy=%lf\nRz=%lf\ndX=%lf\ndY=%lf\ndZ=%lf\n", m,Rx,Ry,Rz,dX,dY,dZ);
    fclose(o);
}
