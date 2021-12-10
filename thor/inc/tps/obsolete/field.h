 /* Author:	 Johan Bengtsson

   Definitions:  Polymorphic number class.              */

#ifndef _TPS_FIELD_H_
#define _TPS_FIELD_H_ 1
#include <tps/ss_vect.h>


// Class for single particle phase space dynamics

// pre-declare template friend functions: f<>()
template<typename T>
ss_vect<T> operator+(const ss_vect<T> &);
template<typename T>
ss_vect<T> operator-(const ss_vect<T> &);



typedef struct MNF_struct
{
  tps
    K,       // New effective Hamiltonian.
    g;       // Generator for nonlinear transformation to Floquet space.
  ss_vect<tps>
    A0,      // Transformation to fixed point.
    A1,      // Transformation (linear) to Floquet space.
    map_res; // Residual map.
} MNF_struct;


template<>
inline ss_vect<double> ss_vect<tps>::cst(void) const
{
  int             i;
  ss_vect<double> x;

  for (i = 0; i < ps_dim; i++)
    x[i] = (*this)[i].cst();
  return x;
}


// partial template-class specialization
// primary version: is_double<>
template<typename T>
class is_double { };

// partial specialization
template<>
class is_double<double> {
 public:
  static inline double cst(const double x) { return x; }
};

// partial specialization
template<>
class is_double<tps> {
 public:
  static inline double cst(const tps &x) { return x.cst(); }
};

// partial specialization
template<>
class is_double< ss_vect<double> > {
 public:
  static inline ss_vect<double> cst(const ss_vect<double> &x) { return x; }
  static inline ss_vect<double> ps(const ss_vect<tps> &x) { return x.cst(); }
};

// partial specialization
template<>
class is_double< ss_vect<tps> > {
 public:
  static inline ss_vect<double> cst(const ss_vect<tps> &x) { return x.cst(); }
  static inline ss_vect<tps> ps(const ss_vect<tps> &x) { return x; }
};


template<typename T>
inline T sqr(const T &a) { return a*a; }


inline tps operator+(const tps &a, const tps &b) { return tps(a) += b; }

inline tps operator+(const tps &a, const double b) { return tps(a) += b; }

inline tps operator+(const double a, const tps &b) { return tps(a) += b; }

inline tps operator-(const tps &a, const tps &b) { return tps(a) -= b; }

inline tps operator-(const tps &a, const double b) { return tps(a) -= b; }

inline tps operator-(const double a, const tps &b) { return tps(a) -= b; }

inline tps operator*(const tps &a, const tps &b) { return tps(a) *= b; }

inline tps operator*(const tps &a, const double b) { return tps(a) *= b; }

inline tps operator*(const double a, const tps &b) { return tps(a) *= b; }

inline tps operator/(const tps &a, const tps &b) { return tps(a) /= b; }

inline tps operator/(const tps &a, const double b) { return tps(a) /= b; }

inline tps operator/(const double a, const tps &b) { return tps(a) /= b; }


inline tps operator+(const tps &x) { return tps(x); }

inline tps operator-(const tps &x) { return tps(x) *= -1.0; }


inline bool operator>(const tps &a, const tps &b) { return a.cst() > b.cst(); }

inline bool operator>(const tps &a, const double b) { return a.cst() > b; }

inline bool operator>(const double a, const tps &b) { return a > b.cst(); }


inline bool operator<(const tps &a, const tps &b) { return a.cst() < b.cst(); }

inline bool operator<(const tps &a, const double b) { return a.cst() < b; }

inline bool operator<(const double a, const tps &b) { return a < b.cst(); }


inline bool operator>=(const tps &a, const tps &b)
{ return a.cst() >= b.cst(); }

inline bool operator>=(const tps &a, const double b) { return a.cst() >= b; }

inline bool operator>=(const double a, const tps &b) { return a >= b.cst(); }


inline bool operator<=(const tps &a, const tps &b)
{ return a.cst() <= b.cst(); }

inline bool operator<=(const tps &a, const double b) { return a.cst() <= b; }

inline bool operator<=(const double a, const tps &b) { return a <= b.cst(); }


inline bool operator==(const tps &a, const tps &b)
{ return a.cst() == b.cst(); }

inline bool operator==(const tps &a, const double b) { return a.cst() == b; }

inline bool operator==(const double a, const tps &b) { return a == b.cst(); }


inline bool operator!=(const tps &a, const tps &b)
{ return a.cst() != b.cst(); }

inline bool operator!=(const tps &a, const double b) { return a.cst() != b; }

inline bool operator!=(const double a, const tps &b) { return a != b.cst(); }


template<typename T>
inline ss_vect<T>::ss_vect(void) { ss.resize(ps_dim); }

template<typename T>
inline ss_vect<T>::ss_vect(const double x, const double px, const double y,
			   const double py, const double delta, const double ct)
{ ss = {x, px, y, py, delta, ct}; }

template<typename T>
template<typename U>
inline ss_vect<T>::ss_vect(const U &a) : ss(a) { }

template<typename T>
template<typename U>
inline ss_vect<T>::ss_vect(const ss_vect<U> &a)
{ for (int i = 0; i < ps_dim; i++) ss.push_back(a[i]); }


template<typename T>
ss_vect<T> ss_vect<T>::zero(void);
/*
{
  for (int i = 0; i < ps_dim; i++)
    ss[i] = 0e0;
  return *this;
}
*/
#endif /* _TPS_FIELD_H_ */
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
