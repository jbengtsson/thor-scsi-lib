#include <thor_scsi/core/elements.h>
#include <cstdio>
const static bool traceID = false;

namespace tse = thor_scsi::elements;

template<typename T>
void LinearInterpolation2(T &X, T &Z, T &TX, T &TZ, T &B2,
			  tse::ElemType *Cell, bool &out, int order)
{
  int            i, ix = 0, iz = 0;
  T              T1, U, THX = 0.0, THZ = 0.0;
  double         xstep = 0.0;
  double         zstep = 0.0;
  int            nx = 0, nz = 0;
  tse::InsertionType  *WITH;

  WITH = dynamic_cast<tse::InsertionType*>(Cell);
  nx = WITH->nx; nz = WITH->nz;

  xstep = WITH->tabx[1]-WITH->tabx[0]; /* increasing values */
  zstep = WITH->tabz[0]-WITH->tabz[1]; /* decreasing values */

  if (traceID) printf("xstep = % f zstep = % f\n", xstep, zstep);

  /* test wether X and Z within the transverse map area */
  if (X < WITH->tabx[0] || X > WITH->tabx[nx-1]) {
    printf("LinearInterpolation2: X out of borders \n");
    printf("X = % lf but tabx[0] = % lf and tabx[nx-1] = % lf\n",
	   is_double<T>::cst(X), WITH->tabx[0], WITH->tabx[nx-1]);
    out = true;
    return;
  }

  if (Z > WITH->tabz[0] || Z < WITH->tabz[nz-1]) {
    printf("LinearInterpolation2: Z out of borders \n");
    printf("Z = % lf but tabz[0] = % lf and tabz[nz-1] = % lf\n",
	   is_double<T>::cst(Z),  WITH->tabz[0], WITH->tabz[nz-1]);
    out = true;
    return;
  }

  out = false;

  /* looking for the index for X */
  i = 0;
  while (X >= WITH->tabx[i]  && i <= nx-1) {
    i++;
    if (traceID)
      printf("%2d % lf % lf % lf \n",
	     i, is_double<T>::cst(X), WITH->tabx[i], WITH->tabx[i+1]);
  }
  ix = i - 1;

  /* looking for the index for Z */
  i = 0;
  while (Z <= WITH->tabz[i] && i <= nz-1) {
    i++;
    if (traceID)
      printf("%2d % lf % lf % lf \n",
	     i, is_double<T>::cst(Z), WITH->tabz[i], WITH->tabz[i+1]);
  }
  iz = i - 1;

  if (traceID) printf("Indices are ix=%d and iz=%d\n", ix, iz);

  /** Bilinear Interpolation **/
  U = (X - WITH->tabx[ix])/xstep; T1 = -(Z - WITH->tabz[iz])/zstep;

  if (order == 1) { // first order kick map interpolation
    if (traceID) printf("first order kick map interpolation\n");
    if (ix >= 0 && iz >= 0) {
      THX = (1.0-U)*(1.0-T1)*WITH->thetax1[iz][ix]
	    + U*(1.0-T1)*WITH->thetax1[iz][ix+1]
	    + (1.0-U)*T1*WITH->thetax1[iz+1][ix]
	    + U*T1*WITH->thetax1[iz+1][ix+1];

      THZ = (1.0-U)*(1.0-T1)*WITH->thetaz1[iz][ix]
	    + U*(1.0-T1)*WITH->thetaz1[iz][ix+1]
	    + (1.0-U)*T1*WITH->thetaz1[iz+1][ix]
	    + U*T1*WITH->thetaz1[iz+1][ix+1];
    }

    if (traceID) {
      printf("X=% f interpolation : U= % lf T =% lf\n",
	     is_double<T>::cst(X), is_double<T>::cst(U),
	     is_double<T>::cst(T1));
      printf("THX = % lf 11= % lf 12= %lf 21 = %lf 22 =%lf \n",
	     is_double<T>::cst(THX),
	     WITH->thetax1[iz][ix], WITH->thetax1[iz][ix+1],
	     WITH->thetax1[iz+1][ix], WITH->thetax1[iz+1][ix+1]);
      printf("Z=% f interpolation : U= % lf T =% lf\n",
	     is_double<T>::cst(Z), is_double<T>::cst(U),
	     is_double<T>::cst(T1));
      printf("THZ = % lf 11= % lf 12= %lf 21 = %lf 22 =%lf \n",
	     is_double<T>::cst(THZ),
	     WITH->thetaz1[iz][ix], WITH->thetaz1[iz][ix+1],
	     WITH->thetaz1[iz+1][ix],WITH->thetaz1[iz+1][ix+1]);
    }
  }

  if (order == 2) { // second order kick map interpolation
    if (traceID) printf("second order kick map interpolation\n");
    if (ix >= 0 && iz >= 0) {
      THX =
	(1.0-U)*(1.0-T1)*WITH->thetax[iz][ix]
	+ U*(1.0-T1)*WITH->thetax[iz][ix+1]
	+ (1.0-U)*T1*WITH->thetax[iz+1][ix]
	+ U*T1*WITH->thetax[iz+1][ix+1];

      THZ =
	(1.0-U)*(1.0-T1)*WITH->thetaz[iz][ix]
	+ U*(1.0-T1)*WITH->thetaz[iz][ix+1]
	+ (1.0-U)*T1*WITH->thetaz[iz+1][ix]
	+ U*T1*WITH->thetaz[iz+1][ix+1];

      if (WITH->long_comp)
	B2 =
	  (1.0-U)*(1.0-T1)*WITH->B2[iz][ix]
	  + U*(1.0-T1)*WITH->B2[iz][ix+1]
	  + (1.0-U)*T1*WITH->B2[iz+1][ix]
	  + U*T1*WITH->B2[iz+1][ix+1];
    }

    if (traceID) {
      printf("X=% f interpolation : U= % lf T =% lf\n",
	     is_double<T>::cst(X), is_double<T>::cst(U),
	     is_double<T>::cst(T1));
      printf("THX = % lf 11= % lf 12= %lf 21 = %lf 22 =%lf \n",
	     is_double<T>::cst(THX),
	     WITH->thetax[iz][ix], WITH->thetax[iz][ix+1],
	     WITH->thetax[iz+1][ix], WITH->thetax[iz+1][ix+1]);
      printf("Z=% f interpolation : U= % lf T =% lf\n",
	     is_double<T>::cst(Z), is_double<T>::cst(U),
	     is_double<T>::cst(T1));
      printf("THZ = % lf 11= % lf 12= %lf 21 = %lf 22 =%lf \n",
	     is_double<T>::cst(THZ),
	     WITH->thetaz[iz][ix], WITH->thetaz[iz][ix+1],
	     WITH->thetaz[iz+1][ix], WITH->thetaz[iz+1][ix+1]);
      printf("B2 = % lf 11= % lf 12= %lf 21 = %lf 22 =%lf \n",
	     is_double<T>::cst(THZ),
	     WITH->B2[iz][ix], WITH->B2[iz][ix+1],
	     WITH->B2[iz+1][ix],WITH->B2[iz+1][ix+1]);
    }
  }
  TX = THX; TZ = THZ;
}


/****************************************************************************/
/* void SplineInterpolation2(double X, double Z, double &TX, double &TZ,
                             ElemType *Cell, bool &out)

   Purpose:
        Computes thx and thz in X and Z values using a bilinear interpolation
        interpolation of the array thetax(x, z) and thetaz(x, z)

   Input:
       X, Z location of the interpolation
       Cell elment containing ID device

   Output:
       TX, TZ thetax and thetaz interpolated at X and Z
       out true if interpolation out of table

   Return:
       none

   Global variables:
       none

   Specific functions:

   Comments:
       none

****************************************************************************/
template<typename T>
void SplineInterpolation2(T &X, T &Z, T &thetax, T &thetaz,
			  tse::ElemType *Cell, bool &out)
{
    int            nx, nz;
    tse::InsertionType  *WITH;
//    int kx, kz;

    WITH = dynamic_cast<tse::InsertionType*>(Cell);
    nx = WITH->nx; nz = WITH->nz;

    /* test wether X and Z within the transverse map area */
    if (X < WITH->tabx[0] || X > WITH->tabx[nx-1] ||
	Z > WITH->tabz[0] || Z < WITH->tabz[nz-1]) {
        printf("SplineInterpDeriv2: out of borders in element s= %4.2f %*s\n",
	       Cell->S, 5, Cell->Name.c_str());
        printf("X = % lf but tabx[0] = % lf and tabx[nx-1] = % lf\n",
	       is_double<T>::cst(X), WITH->tabx[0], WITH->tabx[nx-1]);
        printf("Z = % lf but tabz[0] = % lf and tabz[nz-1] = % lf\n",
	       is_double<T>::cst(Z), WITH->tabz[0], WITH->tabz[nz-1]);
        out = true;
        return;
    }

    out = false;
    splin2(WITH->tab2-1, WITH->tab1-1, WITH->tx, WITH->f2x, nz, nx,
	   Z, X, thetax);
/*    if (fabs(temp) > ZERO_RADIA)
      *thetax = (double) temp;
    else
      *thetax = 0.0;*/
    splin2(WITH->tab2-1, WITH->tab1-1, WITH->tz, WITH->f2z, nz, nx,
	   Z, X, thetaz);
/*    if (fabs(temp) > ZERO_RADIA)
      *thetaz = (double) temp;
    else
      *thetaz = 0.0;*/

/*    FILE * fic0;
    char *fic="fit.out";
    fic0 = fopen(fic, "w");
    for (kz = 1; kz <= nz; kz++) {
      for (kx = 1; kx <= nx; kx++)
	fprintf(fic0, "% 12.3e", tz[kz][kx]);
      fprintf(fic0, "\n");
    }
    fclose(fic0);*/
}
