#ifndef _TPS_EXCEPTIONS_H_
#define _TPS_EXCEPTIONS_H_
#include <exception>
class NotImplemented: public std::exception
{
  virtual const char* what() const throw()
  {
    return "not implemented";
  }
};

#endif /* _TPS_EXCEPTIONS_H_ */
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
