#include <thor_scsi/math/interpolation.h>
#include <exception>
#include <iostream>
//namespace ts = thor_scsi;
//namespace tsc = thor_scsi::core;
//namespace tsm = thor_scsi::math;


template<typename T>
void spline(const double x[], const T y[], int const n,
	    double const yp1, const double ypn, T y2[])
{
  int    i, k;
  double sig;
  // Variable length arrays for user defined data types is not supported by
  // ANSI C++.
  // T      p, u[n], qn, un;
  T      p, qn, un;
  T      *u = new T[n];

  if (yp1 > 0.99e30)
    y2[1] = u[1] = 0.0;
  else {
    y2[1] = -0.5; u[1] = (3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
  }
  for (i = 2; i <= n-1; i++) {
    sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]); p = sig*y2[i-1]+2.0;
    y2[i] = (sig-1.0)/p;
    u[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i] = (6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  if (ypn > 0.99e30)
    qn = un = 0.0;
  else {
    qn = 0.5; un = (3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
  }
  y2[n] = (un-qn*u[n-1])/(qn*y2[n-1]+1.0);
  for (k = n-1; k >= 1; k--)
    y2[k] = y2[k]*y2[k+1]+u[k];

  delete [] u;
}


template<typename T, typename U>
void splint(const double xa[], const U ya[], const U y2a[],
	    const int n, const T &x, T &y)
{
  int     klo, khi, k;
  double  h;
  T       a, b;

  klo = 1; khi = n;
  while (khi-klo > 1) {
    k = (khi+klo) >> 1;
    if (xa[k] > x)
      khi = k;
    else
      klo = k;
  }
  h = xa[khi]-xa[klo];
  if (h == 0.0){
    std::cerr << "Bad xa input to routine splint (step difference 0)" << std::endl;
    throw std::invalid_argument("Bad xa input to routine splint");
  }
  a = (xa[khi]-x)/h; b = (x-xa[klo])/h;
  y = a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}


template<typename T>
void splin2(const double x1a[], const double x2a[],
	    double **ya, double **y2a, const int m, const int n,
	    const T &x1, const T &x2, T &y)
{
  int j;
  // Variable length arrays for user defined data types is not supported by
  // ANSI C++.
  // T   ytmp[m+1], yytmp[m+1];
  T   *ytmp = new T[m+1], *yytmp = new T[m+1];

  for (j = 1; j <= m; j++)
    splint(x2a, ya[j], y2a[j], n, x2, yytmp[j]);
  spline(x1a, yytmp, m, 1.0e30, 1.0e30, ytmp);
  splint(x1a, yytmp, ytmp, m, x1, y);

  delete [] ytmp; delete [] yytmp;
}


void splie2(double x1a[], double x2a[], double **ya,
	    int m, int n, double **y2a)

{
  int  j;

  for (j = 1; j <= m; j++)
    spline(x2a, ya[j], n, 1.0e30, 1.0e30, y2a[j]);
}
