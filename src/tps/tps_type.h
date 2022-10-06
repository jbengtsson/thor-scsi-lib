#ifndef _TPS_TPS_TYPE_H_
#define _TPS_TPS_TYPE_H_ 1

/**
 *
 * Todo:
 *    add a little explanation what tps is all about
 */
#include <tps/config.h>
#include <string>
#include <vector>
#include <array>
#include <ostream>
#include <tps/enums.h>
#include <tps/forward_decl.h>

typedef std::array<long int, ss_dim> tpsa_index;


class tps {
 public:
  tps(void);
  tps(const double);
  tps(const double, const int);
  // replace const long int with std::array<const long int, 7>
  // check length (nv_tps)
  tps(const double, const long int []);
  tps(const tps &);

  tps clone(void) const;
  ~tps(void);

  void print(const std::string &);
    /**
   * @brief show / print to ostream
   * @ærg precision: number of digits to use
   * @ærg with_endl: add an std::endl at the end of the stream (defaults to true)
   * @todo: consider if adding an unused  level flag for consistency to flame
   */
  void show(std::ostream&, int precision=6, bool with_endl=false) const;

  // For python interface
  std::string repr(void);
  // For python interface
  std::string pstr(void);

  // initialize TPSA library
   friend void TPSAEps(const double);
  // trace level for TPSALib and LieLib
  friend void idprset(const int);

  // Constant ness of cst return ignored by callers
  // const double cst(void) const;
  double cst(void) const;
  double operator[](const int) const;


  double operator[](const long int []) const;
  void pook(const long int [], const double);

  double operator[](const tpsa_index) const;
  void pook(const tpsa_index, const double);

  /**
   * Todo:
   *     Rename to export
   */
  void exprt(double [], long int [], long int [], char []) const;
  /**
   * Todo:
   *     Rename to import
   */
  void imprt(const int, double [], const long int [], const long int []);

  tps& operator=(const double);
  tps& operator+=(const double);
  tps& operator-=(const double);
  tps& operator*=(const double);
  tps& operator/=(const double);

  tps& operator=(const tps &);
  tps& operator+=(const tps &);
  tps& operator-=(const tps &);
  tps& operator*=(const tps &);
  tps& operator/=(const tps &);

#if NO_TPSA == 1
  // consider to make these functions inline
  friend double get_m_ij(const ss_vect<tps> &map, const int i, const int j);
  friend void put_m_ij(ss_vect<tps> &map, const int i, const int j,
		       const double r);
  // Same functionallity as above but using .at internally.
  friend  double get_m_ij_save(const ss_vect<tps> &map, const int i, const int j);
  friend void put_m_ij_save(ss_vect<tps> &map, const int i, const int j,
			    const double r);

  friend void dacct_(const ss_vect<tps> &x, const int i,
                     const ss_vect<tps> &y, const int j,
                     ss_vect<tps> &z, const int k);
#endif

  friend std::istream& operator>>(std::istream &, tps &);
  friend std::ostream& operator<<(std::ostream &, const tps &);

  friend double abs(const tps &);
  friend double abs2(const tps &);
  friend tps sqrt(const tps &);
  // refrain from defining sqr as friend class: it is defined as inline below
  // otherwise linking could fail
  friend tps pow(const tps &, const int);
  friend tps exp(const tps &);
  friend tps log(const tps &);
  friend tps sin(const tps &);
  friend tps cos(const tps &);
  friend tps tan(const tps &);
  friend tps asin(const tps &);
  friend tps acos(const tps &);
  friend tps atan(const tps &);
  friend tps sinh(const tps &);
  friend tps cosh(const tps &);

  friend tps Der(const tps &, const int);
  friend tps LieExp(const tps &, const tps &);
  friend tps LieFlo(const ss_vect<tps> &, const tps &);
  friend tps PB(const tps &, const tps &);
  friend tps Take(const tps &, const int);

  // R(nd2, nv) = P(nd2, nd2)*Q(nd2, nv)
  friend ss_vect<tps> operator*(const ss_vect<tps> &, const ss_vect<tps> &);
  // R(nv, nv) = P(nv, nv)*Q(nv, nv)
  friend void CCT(const tps [], const int, const tps [], const int,
		  tps [], const int);
  friend ss_vect<tps> MTREE(const ss_vect<tps> &);
  friend ss_vect<double> PPUSH(const ss_vect<tps> &, ss_vect<double> &);
  friend tps operator*(const tps &, const ss_vect<tps> &);

  friend ss_vect<tps> FExpo(const tps &, const ss_vect<tps> &,
			    const int, const int, const int);
  friend ss_vect<tps> LieExp(const tps &, const ss_vect<tps> &);
  friend ss_vect<tps> LieFlo(const ss_vect<tps> &, const ss_vect<tps> &);
  // Q(nv, nv) = P(nd2, nd2)^-1
  friend ss_vect<tps> Inv(const ss_vect<tps> &);
  // Q(nv, nv) = P(nv, nv)^-1
  friend ss_vect<tps> Inv_Ext(const ss_vect<tps> &);
  friend ss_vect<tps> PInv(const ss_vect<tps> &, const long int[]);
  friend void GoFix(const ss_vect<tps> &, ss_vect<tps> &,
		    ss_vect<tps> &, const int);
  friend MNF_struct MapNorm(const ss_vect<tps> &, const int);
  friend ss_vect<tps> MapNormF(const ss_vect<tps> &, ss_vect<tps> &,
			       ss_vect<tps> &, ss_vect<tps> &,
			       ss_vect<tps> &, const int, const int);
  friend ss_vect<tps> dHdJ(const tps &);
  friend void CtoR(const tps &, tps &, tps &);
  friend tps RtoC(const tps &, const tps &);
  friend tps LieFact_DF(const ss_vect<tps> &, ss_vect<tps> &);
  friend ss_vect<tps> FlowFact(const ss_vect<tps> &);
  friend tps Intd(const ss_vect<tps> &, const double);
  friend ss_vect<tps> Difd(const tps &, const double);
  friend ss_vect<tps> Taked(const ss_vect<tps> &, const int);

#if NO_TPSA != 1
  // accessing elements
  // functions should be inlined?
  friend double getmat(const ss_vect<tps> &map, const int i, const int j);
  friend void putmat(ss_vect<tps> &map, const int i, const int j,
		     const double r);
#endif

private:
#if NO_TPSA == 1
  // Linear TPSA: cst. & linear terms.
  std::vector<double> ltps{0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0};
#else
  long int intptr; // Index used by Fortran implementation.
#endif
  double r;        // Floating-point calc. if intptr = 0.
}; /* class tps */


inline
std::ostream& operator<<(std::ostream& strm, const tps& s)
{
    s.show(strm);
    return strm;
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



#endif /* _TPS_TPS_TYPE_H_ */
