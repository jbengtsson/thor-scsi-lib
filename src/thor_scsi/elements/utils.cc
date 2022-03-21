#include <thor_scsi/elements/utils.h>
#include <cmath>

namespace tse = thor_scsi::elements;

double tse::thirdroot(const double a)
{
	const double exponent = 1.0 / 3.0;
	return pow(a, exponent);
}

/* a bit historical
double tse::thirdroot(const double a)
{
	int    i;
	double x;

	x = 1e0; i = 0;

	do {
		i++;
		x = (x+a)/(x*x+1e0);
	} while (i != 250);
	return x;
}
*/
