 /* Author:	 Johan Bengtsson

   Definitions:  Polymorphic number class.              */

#include <tps/ss_vect.h>
#include <tps/tps_type.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <math.h>


std::string tps::repr(void)
{
	std::stringstream strm;
	this->show(strm, 16);
	return strm.str();
}

std::string tps::pstr(void)
{
	std::stringstream strm;
	this->show(strm, 6);
	return strm.str();
}

#if NO_TPSA == 1

void tps::show(std::ostream& stream, int precision) const
{
  stream << "cst:" << std::endl
	 << std::setw(14) << ltps[0] << std::endl << "linear:" << std::endl;
  for (int i = 1; i <= nv_tps; i++)
    stream << std::scientific << std::setprecision(precision)
	      << std::setw(14) << ltps[i];
  stream << std::endl;
}

void tps::print(const std::string &str)
{
  std::cout << std::scientific << std::setprecision(6) << str << "";
  tps::show(std::cout);
#if 0
  std::cout << std::scientific << std::setprecision(6) << str << "cst:\n"
	    << std::setw(14) << ltps[0] << "\nlinear:\n";
  for (int i = 1; i <= nv_tps; i++)
    std::cout << std::scientific << std::setprecision(6)
	      << std::setw(14) << ltps[i];
  std::cout << "\n";
#endif
}

template<>
void ss_vect<tps>::show(std::ostream& stream, int precision) const
{

  for (int i = 1; i <= nv_tps; i++) {
    for (int j = 1; j <= nv_tps; j++)
      stream << std::scientific << std::setprecision(precision)
		<< std::setw(14) << get_m_ij(*this, i, j);
    stream << std::endl;
  }
}

template<>
void ss_vect<tps>::print(const std::string &str)
{
  std::cout << str;
  ss_vect<tps>::show(std::cout, 0);
}


#if 0

template<>
void ss_vect<tps>::print(const std::string &str)
{


  std::cout << std::scientific << std::setprecision(6) << str << "cst:\n"
	    << std::setw(14) << ltps[0] << "\nlinear:\n";
  for (int i = 1; i <= nv_tps; i++)
    std::cout << std::scientific << std::setprecision(6)
	      << std::setw(14) << ltps[i];
  std::cout << "\n";

}
#endif


#else /* NO_TPSA */
#error "Code completly unchecked here ... needs to be checked"

void tps::print(const std::string &str) { std::cout << str << *this; }

template<>
void ss_vect<tps>::print(const std::string &str) { std::cout << str << *this; }

#endif

template<typename T>
std::string ss_vect<T>::repr(void)
{
	std::stringstream strm;
	this->show(strm, 16);
	return strm.str();
}

template<typename T>
std::string ss_vect<T>::pstr(void)
{
	std::stringstream strm;
	this->show(strm, 6);
	return strm.str();
}


template<>
void ss_vect<double>::show(std::ostream& stream, int precision) const
{
  std::cout << std::scientific << std::setprecision(precision) << std::setw(14)
	    << *this << std::endl;
}

template<>
void ss_vect<double>::print(const std::string &str)
{
  std::cout << std::scientific << std::setprecision(6) << str;
  ss_vect<double>::show(std::cout);
}

// template std::string ss_vect<double>::repr(void);
// template std::string ss_vect<tps>::repr(void);

template<typename T>
ss_vect<T>& ss_vect<T>::operator=(const ss_vect<T> &x)
{
  for (int i = 0; i < ps_dim; i++)
    (*this)[i] = x[i];
  return *this;
}

template<typename T>
ss_vect<T>& ss_vect<T>::operator+=(const ss_vect<T> &a)
{
  for (int i = 0; i < ps_dim; i++)
    ss[i] += a[i];
  return *this;
}

template<typename T>
ss_vect<T>& ss_vect<T>::operator-=(const ss_vect<T> &a)
{
  for (int i = 0; i < ps_dim; i++)
    ss[i] -= a[i];
  return *this;
}

template<typename T>
ss_vect<T>& ss_vect<T>::operator*=(const double a)
{
  for (int i = 0; i < ps_dim; i++)
    ss[i] *= a;
  return *this;
}

template<>
ss_vect<tps>& ss_vect<tps>::operator*=(const tps &a)
{
  for (int i = 0; i < ps_dim; i++)
    ss[i] *= a;
  return *this;
}


template<typename T>
ss_vect<T> operator+(const ss_vect<T> &x) { return ss_vect<T>(x); }

template<typename T>
ss_vect<T> operator-(const ss_vect<T> &x) { return ss_vect<T>(x) *= -1; }

// instantiate
template class ss_vect<double>;
template class ss_vect<tps>;

template ss_vect<double> operator-(const ss_vect<double> &);
template ss_vect<tps> operator-(const ss_vect<tps> &);

ss_vect<double> operator+(const ss_vect<double> &a, const ss_vect<double> &b)
{ return ss_vect<double>(a) += b; }

ss_vect<tps> operator+(const ss_vect<tps> &a, const ss_vect<double> &b)
{ return ss_vect<tps>(a) += b; }

ss_vect<tps> operator+(const ss_vect<double> &a, const ss_vect<tps> &b)
{ return ss_vect<tps>(a) += b; }

ss_vect<tps> operator+(const ss_vect<tps> &a, const ss_vect<tps> &b)
{ return ss_vect<tps>(a) += b; }

ss_vect<double> operator-(const ss_vect<double> &a, const ss_vect<double> &b)
{ return ss_vect<double>(a) -= b; }

ss_vect<tps> operator-(const ss_vect<tps> &a, const ss_vect<double> &b)
{ return ss_vect<tps>(a) -= b; }

ss_vect<tps> operator-(const ss_vect<double> &a, const ss_vect<tps> &b)
{ return ss_vect<tps>(a) -= b; }

ss_vect<tps> operator-(const ss_vect<tps> &a, const ss_vect<tps> &b)
{ return ss_vect<tps>(a) -= b; }

ss_vect<double> operator*(const ss_vect<double> &a, const double b)
{ return ss_vect<double>(a) *= b; }

ss_vect<double> operator*(const double a, const ss_vect<double> &b)
{ return ss_vect<double>(b) *= a; }

ss_vect<tps> operator*(const ss_vect<tps> &a, const double b)
{ return ss_vect<tps>(a) *= b; }

ss_vect<tps> operator*(const double a, const ss_vect<tps> &b)
{ return ss_vect<tps>(b) *= a; }


template<typename T>
ss_vect<T> ss_vect<T>::zero(void)
{
  for (int i = 0; i < ps_dim; i++)
    ss[i] = 0e0;
  return *this;
}

template ss_vect<tps> ss_vect<tps>::zero(void);
template ss_vect<double> ss_vect<double>::zero(void);

template<>
ss_vect<double> ss_vect<double>::identity(void)
{
  printf("\nidentity: not implemented for ss_vect<double>\n");
  exit(1);
}

template<>
ss_vect<tps> ss_vect<tps>::identity(void)
{
  for (int i = 0; i < ps_dim; i++)
   ss[i] = tps(0e0, i+1);
  return *this;
}


template<typename CharT, class Traits>
std::basic_istream<CharT, Traits>&
operator>>(std::basic_istream<CharT, Traits> &is, ss_vect<tps> &a)
{
  tps b;

  for (int i = 0; i < ps_dim; i++)
    is >> a[i];
  return is;
}

// instantiate
template std::basic_istream<char, std::char_traits<char> >&
operator>>(std::basic_istream<char, std::char_traits<char> > &, ss_vect<tps> &);

template<typename CharT, class Traits>
std::basic_ostream<CharT, Traits>&
operator<<(std::basic_ostream<CharT, Traits> &os, const ss_vect<double> &a)
{
  std::basic_ostringstream<CharT, Traits>  s;

  s.flags(os.flags()); s.imbue(os.getloc());
  for (int i = 0; i < ps_dim; i++)
    s << std::setprecision(os.precision()) << std::setw(os.width()) << a[i];
//   s << endl;
  return os << s.str();
}

// instantiate
template std::basic_ostream<char, std::char_traits<char> >&
operator<<(std::basic_ostream<char, std::char_traits<char> > &,
	   const ss_vect<double> &);

template<typename CharT, class Traits>
std::basic_ostream<CharT, Traits>&
operator<<(std::basic_ostream<CharT, Traits> &os, const ss_vect<tps> &a)
{
  std::basic_ostringstream<CharT, Traits>  s;

  s.flags(os.flags()); s.imbue(os.getloc());
  for (int i = 0; i < ps_dim; i++)
    s << std::setprecision(os.precision()) << std::setw(os.width()) << a[i];
  return os << s.str();
}

// instantiate
template std::basic_ostream<char, std::char_traits<char> >&
operator<<(std::basic_ostream<char, std::char_traits<char> > &,
	   const ss_vect<tps> &);

double xabs(long n, ss_vect<double> &x)
{
  long    i;
  double  sum;

  sum = 0.0;
  for (i = 0; i < n; i++)
    sum += sqr(x[i]);

  return sqrt(sum);
}
