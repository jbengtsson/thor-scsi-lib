
#ifndef TRACY_H
#define TRACY_H
#include <thor_scsi/version.h>




#define HOMmax   21     // [a_n, b_n] <=> [-HOMmax..HOMmax].

#define IDXMAX  200
#define IDZMAX  100

#define fitvectmax 200
typedef long   fitvect[fitvectmax];

extern double Fdrift1, Fkick1, Fdrift2, Fkick2, crad, cfluc;


extern void t2init(void);

extern void prt_gcmat(int bpm, int corr, int plane);

extern void gcmat(int bpm, int corr, int plane);

extern void lsoc(int niter, int bpm, int corr, int plane);






// Beam line class.






void t2init(void);

void exit_(int exit_code);


#endif
