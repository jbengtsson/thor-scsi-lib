#include <armadillo>


#define PS_DIM 6

arma::vec arrtovec(const double a[])
{
  int       k;
  arma::vec a_vec(PS_DIM);

  for (k = 0; k < PS_DIM; k++)
    a_vec(k) = a[k];
  return a_vec;
}


void vectoarr(const arma::vec &a_vec, double a[])
{
  int k;

  for (k = 0; k < PS_DIM; k++)
    a[k] = a_vec(k);
}


arma::mat arrtomat(const double a[])
{
  int       j, k;
  arma::mat A(PS_DIM, PS_DIM);

  for (j = 0; j < PS_DIM; j++)
    for (k = 0; k < PS_DIM; k++)
      A(j, k) = a[j*PS_DIM+k];
  return A;
}


void mattoarr(const arma::mat &A, double a[])
{
  int j, k;

  for (j = 0; j < PS_DIM; j++)
    for (k = 0; k < PS_DIM; k++)
      a[j*PS_DIM+k] = A(j, k);
}


void IdentityPass(double ps_n[],
		  const double t1[], const double t2[],
		  const double r1[], const double r2[],
		  const double limits[], const double *axesptr,
		  const int num_particles)
{
  int       j, k;
  arma::vec ps(PS_DIM);

  const arma::vec
    t1_vec = arrtovec(t1),
    t2_vec = arrtovec(t2);
  const arma::mat
    R1 = arrtomat(r1),
    R2 = arrtomat(r2);

  for (j = 0; j < num_particles; j++) {	// Loop over particles.
    for (k = 0; k < PS_DIM; k++)
      ps(k) = ps_n[j*PS_DIM+k];
    if (!isnan(ps(0))) {
      //  misalignment at entrance.
      if (t1) ps = ps + t1_vec;
      if (r1) ps = R1*ps;
      // Check physical apertures.
      // if (limits) checkiflostRectangularAp(r6,limits);
      // if (axesptr) checkiflostEllipticalAp(r6,axesptr);
      // Misalignment at exit.
      if (r2) ps = R2*ps;
      if (t2) ps = ps + t2_vec;
    }
    for (k = 0; k < PS_DIM; k++)
      ps_n[j*PS_DIM+k] = ps(k);
  }
}
