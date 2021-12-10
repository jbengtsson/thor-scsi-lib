#ifndef _THOR_SCSI_LEGACY_H_
#define _THOR_SCSI_LEGACY_H_ 1

#define Cell_nLocMax    20000 // maximum number of LEGO blocks (Cell_nLoc).

#ifndef LONG_MAX
# define LONG_MAX       ((long)(((unsigned long) -1) >> 1))
# define LONG_MIN       (~LONG_MAX)
#endif

#define maxincl         5
#define maxfil          10

// Dynamic aperture (chk_if_lost).
#define px_0            0e0
#define py_0            0e0


// Macros.

#define degtorad(x) ((x)*M_PI/180.0)
#define radtodeg(x) ((x)*180.0/M_PI)

#define sqr(x)      ((x)*(x))
#define cube(x)     ((x)*(x)*(x))

#define fract(x)    ((x)-(int)(x))
#define nint(x)     ((x) < 0 ? ((long)(x-0.5)) : ((long)(x+0.5)))

#define sgn(n)      ((n > 0) ? 1 : ((n < 0) ? -1 : 0))


namespace thor_scsi {
  namespace legacy {

  }
}
#endif /* _THOR_SCSI_LEGACY_H_ */
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
