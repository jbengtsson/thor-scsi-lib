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

	class InvalidPosition: public std::exception
	{
		virtual const char* what() const throw()
			{
				return "Invalid position";
			}
	};

/**  Last position eg. 0 or 1
 *   Most probably none reached
 */
	class InvalidLastPosition: public InvalidPosition
	{
		virtual const char* what() const throw()
			{
				return "Invalid last position";
			}
	}invalid_last_position;

	class LatticeParseError: public std::exception
	{
		virtual const char* what() const throw()
			{
				return "Lattice parse error";
			}

	}lattice_parse_error;

}
#endif
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
