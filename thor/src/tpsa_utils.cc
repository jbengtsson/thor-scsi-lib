#include <tps/utils.h>
#include <exception>
#include <iostream>

double eps_tps  = 1e-25;

long int fact(long int n)
{
  if (n > 0)
    return n*fact(n-1);
  else if (n == 0)
    return 1;
  else {
    std::cout << "fact: neg. argument: " << n << "\n";
    throw std::domain_error("fact: neg. argument!");

    // avoid compiler warning
    return -1;
  }
}

long int nok(long int n, long int k)
{
  long int j;
  double   u;

  u = 1.0;
  for (j = 0; j < k; j++)
    u *= (double)(n-j)/(double)(k-j);
  return (long int)(u+0.5);
}
