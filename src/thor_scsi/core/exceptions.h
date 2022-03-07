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
	};

	class LatticeParseError: public std::exception
	{
		virtual const char* what() const throw()
			{
				return "Lattice parse error";
			}

	};

	class SanityCheckError:  public std::exception
	{
		virtual const char* what() const throw()
			{
				return "Sanity check error";
			}
	};

        class LoadException: public std::exception {
	private:
		std::string m_message;
	protected:
		std::string m_prefix;
	public:
		explicit LoadException(const std::string& message) : m_message(message) {m_prefix = "";}
		const char* what() const noexcept override {
			return (m_prefix + " " + m_message).c_str();
		}
	};

	class PhysicsViolation : LoadException
	{
	public:
		explicit PhysicsViolation(const std::string& message) : LoadException(message) {
			m_prefix = "Physics domain violated";
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
