/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -

*/

#include <thor_scsi/core/constants.h>
#include <thor_scsi/exceptions.h>
#include <thor_scsi/importers/radia.h>
#include <thor_scsi/math/interpolation.h>
#include <iostream>
#include <math.h>
#include <cstring>

const static bool traceID = false;

namespace ts = thor_scsi;
namespace tse = thor_scsi::elements;
namespace tsc = thor_scsi::core;
namespace tsm = thor_scsi::math;


#define ZERO_RADIA 1e-7

void Read_IDfile(char *fic_radia, const thor_scsi::core::ConfigType &conf, thor_scsi::elements::InsertionType *ID)
{
  char dummy[5000];
  int  i, j;
  FILE *fi;

  const double Brho = conf.Energy*1e9/tsc::c0;

  /* open radia text file */
  if ((fi = fopen(fic_radia,"r")) == NULL) {
    printf("Read_IDfile: Error while opening file %s \n", fic_radia);
    throw ts::LatticeParseError();
    // exit_(1);
  }

  printf("\n");
  printf("Reading ID filename %s \n", fic_radia);
  printf("E      = %6.3f GeV\n", conf.Energy);
  printf("(Brho) = %6.3f\n", Brho);

  /* first line */
  fscanf(fi, "%[^\n]\n", dummy); /* Read a full line */
  printf("%s\n", dummy);
  /* second line */
  fscanf(fi, "%[^\n]\n", dummy);
  printf("%s\n", dummy);
  /* third line */
  fscanf(fi, "%[^\n]\n", dummy);
  printf("%s\n", dummy);
  /* fourth line : Undulator length */
  fscanf(fi, "%lf\n", &ID->PL);
  printf("Insertion de longueur L = %lf m\n", ID->PL);
  /* fifth line */
  fscanf(fi, "%[^\n]\n", dummy);
  printf("%s\n", dummy);
  /* sisxth line : Number of Horizontal points */
  fscanf(fi, "%d\n", &ID->nx);
  printf("nx = %d\n", ID->nx);
  /* seventh line */
  fscanf(fi, "%[^\n]\n", dummy);
  printf("%s\n", dummy);
  /* Number of Vertical points */
  fscanf(fi, "%d\n", &ID->nz);
  printf("nz = %d\n", ID->nz);

  /* Check dimensions */
  if (ID->nx > IDXMAX || ID->nz > IDZMAX) {
    printf("Read_IDfile:  Increase the size of insertion tables \n");
    printf("nx = % d (IDXmax = %d) and nz = %d (IDZMAX = % d) \n",
	   ID->nx, IDXMAX, ID->nz, IDZMAX);
    throw std::length_error("Increase the size of insertion table");
    //exit_(1);
  }

  /* ninth line */
  fscanf(fi, "%[^\n]\n", dummy);
  printf("%s\n", dummy);
  /* tenth line */
  fscanf(fi, "%[^\n]\n", dummy);
//   printf("%s\n", dummy);

  for (j = 0; j < ID->nx; j++)
    fscanf(fi, "%lf", &ID->tabx[j]);
  fscanf(fi, "%[^\n]\n", dummy);

  /* Array creation for thetax */
  for (i = 0; i < ID->nz; i++) {
    fscanf(fi, "%lf", &ID->tabz[i]);

    for (j = 0; j < ID->nx; j++) {
      fscanf(fi, "%lf", &ID->thetax[i][j]);
      if (fabs(ID->thetax[i][j]) < ZERO_RADIA)
	ID->thetax[i][j] = 0.0;
      if (traceID) printf("%+12.8lf ", ID->thetax[i][j]);
    }
    fscanf(fi, "\n");
    if (traceID) printf("\n");
  }

  fscanf(fi, "%[^\n]\n", dummy);
  printf("%s\n", dummy);
  fscanf(fi, "%[^\n]\n", dummy);
//   printf("%s\n", dummy);
  for (j = 0; j < ID->nx; j++) {
    fscanf(fi, "%lf", &ID->tabx[j]);
  }

  /* Array creation for thetaz */
  for (i = 0; i < ID->nz; i++) {
    fscanf(fi, "%*f");
    for (j = 0; j < ID->nx; j++) {
      fscanf(fi, "%lf", &ID->thetaz[i][j]);
      if (fabs(ID->thetaz[i][j]) < ZERO_RADIA)
	ID->thetaz[i][j] = 0.0;
      if (traceID)
	printf("%+12.8lf ", ID->thetaz[i][j]);
    }
    fscanf(fi, "\n");
    if (traceID) printf("\n");
  }

  /* Array creation for B2 */
  strcpy(dummy, "");
  fscanf(fi, "%[^\n]\n", dummy);
  printf("%s\n", dummy);
  fscanf(fi, "%[^\n]\n", dummy);
//  printf("%s\n", dummy);
  if (strstr(dummy, "START") == NULL) {
    ID->long_comp = false;
    printf("no longitudinal component\n");
  } else {
    ID->long_comp = true;
    printf("read longitudinal component\n");

    for (j = 0; j < ID->nx; j++)
      fscanf(fi, "%lf", &ID->tabx[j]);

    for (i = 0; i < ID->nz; i++) {
      fscanf(fi, "%*f");
      for (j = 0; j < ID->nx; j++) {
	fscanf(fi, "%lf", &ID->B2[i][j]);
	ID->B2[i][j] /= ID->PL*sqr(Brho);
	if (fabs(ID->B2[i][j]) < ZERO_RADIA)
	  ID->B2[i][j] = 0e0;
	if (traceID) printf("%+12.8lf ", ID->B2[i][j]);
      }
      fscanf(fi, "\n");
      if (traceID) printf("\n");
    }
  }

  if (traceID)
    for (j = 0; j < ID->nx; j++)
      printf("tabx[%d] =% lf\n", j, ID->tabx[j]);
  if (traceID)
    for (j = 0; j < ID->nz; j++)
      printf("tabz[%d] =% lf\n", j, ID->tabz[j]);

  fclose(fi);
}




void Matrices4Spline(tse::InsertionType *WITH)
{
  int kx, kz;

  for (kx = 0; kx < WITH->nx; kx++) {
    WITH->tab1[kx] = (float) WITH->tabx[kx];
  }

  /** reordering: it has to be in increasing order */
  for (kz = 0; kz < WITH->nz; kz++) {
    WITH->tab2[kz] = (float) WITH->tabz[WITH->nz-kz-1];
  }

  for (kx = 0; kx < WITH->nx; kx++) {
    for (kz = 0; kz <WITH-> nz; kz++) {
      WITH->tx[kz+1][kx+1] = (float) (WITH->thetax[WITH->nz-kz-1][kx]);
      WITH->tz[kz+1][kx+1] = (float) (WITH->thetaz[WITH->nz-kz-1][kx]);
    }
  }

  // computes second derivative matrices
  tsm::splie2(WITH->tab2-1,WITH->tab1-1,WITH->tx,WITH->nz,WITH->nx,WITH->f2x);
  tsm::splie2(WITH->tab2-1,WITH->tab1-1,WITH->tz,WITH->nz,WITH->nx,WITH->f2z);
}
