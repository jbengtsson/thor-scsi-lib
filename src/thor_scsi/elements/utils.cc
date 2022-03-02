#include <thor_scsi/elements/utils.h>

namespace tse = thor_scsi::elements;

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
