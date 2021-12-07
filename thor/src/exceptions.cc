#include <exception>

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
