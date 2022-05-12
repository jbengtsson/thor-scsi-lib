#include <thor_scsi/elements/mpole.h>
#include <ostream>

namespace tse = thor_scsi::elements;

int main(int argc, char *argv[])
{
	Config C;
	C.set<std::string>("name", "test");

	{
		C.set<double>("L", 0.0);
		tse::FieldKick kick(C);
	}
	return 0;
}
