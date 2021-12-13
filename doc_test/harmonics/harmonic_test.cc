#include "harmonics.h"
#include <iostream>
#include "math.h"

int main(int argc, char * argv[])
{

	namespace ts = thor_scsi;
	ts::Harmonics h;

	h.setHarmonic(10, 4);
	auto h2 = h;
	h.setHarmonic(1, 3);
	h.setHarmonic(.1, 1);
	h.setHarmonic(.5, 2);

	for (auto coeff: *h.getCoeffs()){
		std::cout <<  coeff << ", ";
	}
	std::cout << std::endl;

	double angle = M_PI/2./3.;
	h.applyRollAngle(angle);

	std::cout << "Converted by " << angle * 180 / M_PI << " degree ";
	for (auto coeff: *h.getCoeffs()){
		std::cout <<  coeff << ", ";
	}
	std::cout << std::endl;

	h2.applyTranslation(.1);
	std::cout << "Translated by .1 ";
	for (auto coeff: *h2.getCoeffs()){
		std::cout <<  coeff << ", ";
	}
	std::cout << std::endl;

	auto h3 = h2 + h;
}
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
