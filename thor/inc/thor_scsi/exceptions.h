#ifndef _THOR_SCSI_EXCEPTIONS_H_
#define _THOR_SCSI_EXCEPTIONS_H_
#include <exception>

namespace thor_scsi {

	class NotImplemented: public std::exception
	{
		virtual const char* what() const throw()
			{
				return "not implemented";
			}
	};
}
#endif
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
