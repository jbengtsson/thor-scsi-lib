#ifndef _THOR_SCSI_EXCEPTIONS_H_
#define _THOR_SCSI_EXCEPTIONS_H_
#include <exception>

namespace thor_scsi {

#if 0
	class NotImplemented: public std::exception
	{
		virtual const char* what() const throw()
			{
				return "not implemented";
			}
	};
#endif
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

	class LoadException: public std::exception {
	private:
		std::string m_message;
	public:
		explicit LoadException(const std::string& message, const std::string& prefix = "") :
			m_message(prefix + " " + message)
			{}
		const char* what() const noexcept override {
			return this->m_message.c_str();
		}
	};

	class SanityCheckError:  public LoadException
	{
	public:
		explicit SanityCheckError(const std::string& message) :
			LoadException(message, "Sanity check error") { }
/*
		virtual const char* what() const throw()
			{
				return "Sanity check error";
			}
*/
	};

	class NotImplemented:  public LoadException
	{
	public:
		explicit NotImplemented(const std::string& message) :
			LoadException(message, "Sanity check error") { }
		/*
		virtual const char* what() const throw()
			{
				return "Sanity check error";
			}
		*/
	};


	class PhysicsViolation : public LoadException
	{
	public:
		explicit PhysicsViolation(const std::string& message) :
			LoadException(message, "Physics domain violated") { }
	};



}
#endif
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
