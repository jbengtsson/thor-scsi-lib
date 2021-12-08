#include <thor_scsi/tweak/errors.h>
#include <string>

namespace ts = thor_scsi;
namespace tsc = thor_scsi::core;
namespace tse = thor_scsi::elements;

const static bool normal = true; // Normal or rectangular distribution.

const static long int
  k_ = 19,
  c_ = 656329L,
  m_ = 100000001;

const static int
  maxiter = 100;

static long int
  rseed0,
  rseed;
static double
  normcut_;

void iniranf(const long i) { rseed0 = i; rseed = i; }

void newseed(void)
{
  rseed0 = (k_*rseed0+c_) % m_; rseed = (rseed0+54321) % m_;
}

// Random number generator with rectangular distribution.
double ranf(void) { rseed = (k_*rseed+c_) % m_; return (rseed/1e8); }

void setrancut(const double cut)
{
  printf("\nsetrancut: cut set to %3.1f\n", cut);
  normcut_ = cut;
}

double normranf(void)
{
  int    i, j;
  double f, w;

  j = 0;
  do {
    j++;
    w = 0.0;
    for (i = 1; i <= 12; i++)
      w += ranf();
    f = w - 6.0;
  }
  while (fabs(f) > fabs(normcut_) && j <= maxiter);

  if (j > maxiter) {
    printf("normranf: algorithm did not converge\n");
    exit(1);
  }
  return f;
}


// Misalignments.

void CheckAlignTol(tsc::LatticeType &lat, const char *OutputFile)
  // check aligment errors of individual magnets on giders
  // the dT and roll angle are all printed out
{
  int          i, j;
  int          n_girders;
  int          gs_Fnum, ge_Fnum;
  int          gs_nKid, ge_nKid;
  int          dip_Fnum,dip_nKid;
  int          loc, loc_gs, loc_ge;
  std::string name;
  //__gnu_debug::string  name;
  double       s;
  double       PdSsys[2], PdSrms[2], PdSrnd[2], dS[2], dT[2];
  tse::MpoleType    *M;
  std::fstream fout;

  gs_Fnum = lat.conf.gs;   gs_nKid = lat.GetnKid(gs_Fnum);
  ge_Fnum = lat.conf.ge;   ge_nKid = lat.GetnKid(ge_Fnum);
  if (gs_nKid == ge_nKid)
    n_girders= gs_nKid;
  else {
    std::cout << " The numbers of GS = "<< gs_nKid  << " and "
	      << "GE = " << ge_nKid << " not same. " << std::endl;
    throw ts::SanityCheckError();
    //exit (1);
  }

  fout.open(OutputFile,std::ios::out);
  if(!fout) {
    std::cout << "error in opening the file  " << OutputFile << std::endl;
    throw std::ios_base::failure("Error opening file");
  }

  fout << "Girders, Quads, Sexts:  " << std::endl;
  for (i = 1; i <= n_girders; i++){
    fout << i << ":" << std::endl;
    loc_gs = lat.Elem_GetPos(gs_Fnum, i); loc_ge = lat.Elem_GetPos(ge_Fnum, i);

    loc = loc_gs;
    M = dynamic_cast<tse::MpoleType*>(lat.elems[loc]);
    PdSsys[X_] = M->PdSsys[X_];
    PdSsys[Y_] = M->PdSsys[Y_];
    PdSrms[X_] = M->PdSrms[X_];
    PdSrms[Y_] = M->PdSrms[Y_];
    PdSrnd[X_] = M->PdSrnd[X_];
    PdSrnd[Y_] = M->PdSrnd[Y_];
    dS[X_] = lat.elems[loc]->dS[X_]; dS[Y_] = lat.elems[loc]->dS[Y_];
    dT[0] = lat.elems[loc]->dT[0]; dT[1] = lat.elems[loc]->dT[1];
    s = lat.elems[loc]->S; name = lat.elems[loc]->Name;
    fout << "  " << name << "  " << loc << "   " << s
	 << "  " <<  PdSsys[X_] << "  " <<  PdSsys[Y_]
	 << "   " << PdSrms[X_] << "  " <<  PdSrms[Y_]
	 << "   " << PdSrnd[X_] << "  " <<  PdSrnd[Y_]
         << "   " << M->PdTrms << "  "
	 << M->PdTrnd << "   " << dS[X_]     << "  " <<  dS[Y_]
	 << "   " << atan2( dT[1], dT[0] )  << std::endl;

    for (j = loc_gs+1; j < loc_ge; j++) {
      loc = j;
      M = dynamic_cast<tse::MpoleType*>(lat.elems[loc]);
      if ((lat.elems[j]->Pkind == tse::Mpole) &&
	  (M->n_design >= tse::Quad || M->n_design >= tse::Sext)) {
	PdSsys[X_] = M->PdSsys[X_];
	PdSsys[Y_] = M->PdSsys[Y_];
	PdSrms[X_] = M->PdSrms[X_];
	PdSrms[Y_] = M->PdSrms[Y_];
	PdSrnd[X_] = M->PdSrnd[X_];
	PdSrnd[Y_] = M->PdSrnd[Y_];
	dS[X_] = lat.elems[loc]->dS[X_]; dS[Y_] = lat.elems[loc]->dS[Y_];
	dT[0] = lat.elems[loc]->dT[0];   dT[1] = lat.elems[loc]->dT[1];
	s = lat.elems[loc]->S; name = lat.elems[loc]->Name;
	fout << "  " << name << "  " << loc << "   " << s
	     << "  " <<  PdSsys[X_] << "  " <<  PdSsys[Y_]
	     << "   " << PdSrms[X_] << "  " <<  PdSrms[Y_]
	     << "   " << PdSrnd[X_] << "  " <<  PdSrnd[Y_]
	     << "   " << M->PdTrms << "  "
	     << M->PdTrnd
	     << "   " << dS[X_] << "  " <<  dS[Y_]
	     << "   " << atan2( dT[1], dT[0] )  << std::endl;
      }
    }

    loc = loc_ge;
    M = dynamic_cast<tse::MpoleType*>(lat.elems[loc]);
    PdSsys[X_] = M->PdSsys[X_];
    PdSsys[Y_] = M->PdSsys[Y_];
    PdSrms[X_] = M->PdSrms[X_];
    PdSrms[Y_] = M->PdSrms[Y_];
    PdSrnd[X_] = M->PdSrnd[X_];
    PdSrnd[Y_] = M->PdSrnd[Y_];
    dS[X_] = lat.elems[loc]->dS[X_]; dS[Y_] = lat.elems[loc]->dS[Y_];
    dT[0] = lat.elems[loc]->dT[0]; dT[1] = lat.elems[loc]->dT[1];
    s=lat.elems[loc]->S; name = lat.elems[loc]->Name;
    fout << "  " << name << "  " << loc << "   " << s
	 << "  " <<  PdSsys[X_] << "  " <<  PdSsys[Y_]
	 << "   " << PdSrms[X_] << "  " <<  PdSrms[Y_]
	 << "   " << PdSrnd[X_] << "  " <<  PdSrnd[Y_]
         << "   " << M->PdTrms
	 << "  " << M->PdTrnd
         << "   " << dS[X_]     << "  " <<  dS[Y_]
	 << "   " << atan2( dT[1], dT[0] )  << std::endl;

  }

  fout << "  " << std::endl;
  fout << "Dipoles:  " << std::endl;
  dip_Fnum = lat.ElemIndex("B1"); dip_nKid = lat.GetnKid(dip_Fnum);
  for (i = 1; i <= dip_nKid; i++){
    loc = lat.Elem_GetPos(dip_Fnum, i);
    M = dynamic_cast<tse::MpoleType*>(lat.elems[loc]);
    PdSsys[X_] = M->PdSsys[X_];
    PdSsys[Y_] = M->PdSsys[Y_];
    PdSrms[X_] = M->PdSrms[X_];
    PdSrms[Y_] = M->PdSrms[Y_];
    PdSrnd[X_] = M->PdSrnd[X_];
    PdSrnd[Y_] = M->PdSrnd[Y_];
    dS[X_] = lat.elems[loc]->dS[X_]; dS[Y_] = lat.elems[loc]->dS[Y_];
    dT[0] = lat.elems[loc]->dT[0]; dT[1] = lat.elems[loc]->dT[1];
    s = lat.elems[loc]->S; name = lat.elems[loc]->Name;
    fout << "  " << name << "  " << loc << "   " << s
	 << "  " <<  PdSsys[X_] << "  " <<  PdSsys[Y_]
	 << "   " << PdSrms[X_] << "  " <<  PdSrms[Y_]
	 << "   " << PdSrnd[X_] << "  " <<  PdSrnd[Y_]
	 << "   " << M->PdTrms
	 << "  " << M->PdTrnd
	 << "   " << dS[X_]     << "  " <<  dS[Y_]
	 << "   " << atan2( dT[1], dT[0] )  << std::endl;
  }

  fout.close();
}


void ts::misalign(tsc::LatticeType &lat, const ts::dx_types dx_type, const int Fnum,
	      const int Knum, const double dx, const double dy, const double dr,
	      const bool new_rnd)
{
  long int  loc;
  tse::MpoleType *Mp;

  loc = lat.Elem_GetPos(Fnum, Knum);
  Mp  = dynamic_cast<tse::MpoleType*>(lat.elems[loc]);

  switch (dx_type) {
  case ts::dx_sys:
  Mp->PdSsys[X_] = dx;
  Mp->PdSsys[Y_] = dy;
  Mp->PdTsys     = dr;

  lat.SetdS(Fnum, Knum); lat.SetdT(Fnum, Knum);
  break;

  case ts::dx_rms:
    Mp->PdSrms[X_] = dx;
    Mp->PdSrms[Y_] = dy;
    Mp->PdTrms     = dr;
    if (new_rnd) {
      if (normal) {
	Mp->PdSrnd[X_] = normranf();
	Mp->PdSrnd[Y_] = normranf();
	Mp->PdTrnd     = normranf();
      } else {
	Mp->PdSrnd[X_] = ranf();
	Mp->PdSrnd[Y_] = ranf();
	Mp->PdTrnd     = ranf();
      }
    }
  break;
  default:
    printf("\nget_bn: undefined error type: %d\n", dx_type);
    exit(1);
    break;
  }

  lat.SetdS(Fnum, Knum); lat.SetdT(Fnum, Knum);
}

void ts::misalign(tsc::LatticeType &lat, const ts::dx_types dx_type, const int Fnum,
	      const double dx, const double dy, const double dr,
	      const bool new_rnd)
{
  int  i;

  for (i = 1; i <= lat.GetnKid(Fnum); i++)
    misalign(lat, dx_type, Fnum, i, dx, dy, dr, new_rnd);
}

void ts::misalign_rms_type(tsc::LatticeType &lat, const int type, const double dx,
		       const double dy, const double dr, const bool new_rnd)
{
  long int  k;
  tse::MpoleType *Mp;

  if ((type >= tse::All) && (type <= HOMmax)) {
    for (k = 1; k <= lat.conf.Cell_nLoc; k++) {
      Mp = dynamic_cast<tse::MpoleType*>(lat.elems[k]);
      if ((lat.elems[k]->Pkind == tse::Mpole) &&
	  ((type == Mp->n_design) || ((type == tse::All) &&
	   ((lat.elems[k]->Fnum != lat.conf.gs)
	    && (lat.elems[k]->Fnum != lat.conf.ge))))) {
	// if all: skip girders
	misalign(lat, ts::dx_rms, lat.elems[k]->Fnum, lat.elems[k]->Knum, dx, dy,
		 dr, new_rnd);
      }
    }
  } else {
    std::cerr << __FILE__ << "." << __FUNCTION__ << ":" << __LINE__
	      << " incorrect type " << type << std::endl;
    throw ts::SanityCheckError();
  }
}

void ts::misalign_sys_type(tsc::LatticeType &lat, const int type, const double dx,
		       const double dy, const double dr)
{
  long int  k;
  tse::MpoleType *Mp;

  if ((type >= tse::All) && (type <= HOMmax)) {
    for (k = 1; k <= lat.conf.Cell_nLoc; k++) {
      Mp = dynamic_cast<tse::MpoleType*>(lat.elems[k]);
      if ((lat.elems[k]->Pkind == tse::Mpole) &&
	  ((type == Mp->n_design) || ((type == tse::All) &&
	   ((lat.elems[k]->Fnum != lat.conf.gs)
	    && (lat.elems[k]->Fnum != lat.conf.ge))))) {
	// if all: skip girders
	misalign(lat, ts::dx_sys, lat.elems[k]->Fnum, lat.elems[k]->Knum, dx, dy,
		 dr, false);
      }
    }
  } else {
    std::cerr << __FILE__ << "." << __FUNCTION__ << ":" << __LINE__
	      << " incorrect type " << type << std::endl;
    throw ts::SanityCheckError();
  }
}

void ts::misalign_rms_girders(tsc::LatticeType &lat, const int gs, const int ge,
			  const double dx, const double dy,
			  const double dr, const bool new_rnd)
{
  int       i, k, n_girders, n_ge, n_gs;
  long int  loc_gs, loc_ge, j;
  double    s_gs, s_ge, dx_gs[2], dx_ge[2], s;
  tse::MpoleType *Mgs, *Mge, *Mj;

  n_gs = lat.GetnKid(gs); n_ge = lat.GetnKid(ge);

  if (n_gs == n_ge)
    n_girders = n_gs;
  else {
    std::cout << "set_girders: no of GS != no of GE" << std::endl;
    exit (1);
  }

  misalign(lat, ts::dx_rms, gs, dx, dy, dr, new_rnd);
  misalign(lat, ts::dx_rms, ge, dx, dy, dr, new_rnd);

  for (i = 1; i <= n_girders; i++) {
    loc_gs = lat.Elem_GetPos(gs, i);
    loc_ge = lat.Elem_GetPos(ge, i);
    s_gs = lat.elems[loc_gs]->S;
    s_ge = lat.elems[loc_ge]->S;

    // roll for a rigid boby
    // Note, girders needs to be introduced as gs->ge pairs
    Mgs = dynamic_cast<tse::MpoleType*>(lat.elems[loc_gs]);
    Mge = dynamic_cast<tse::MpoleType*>(lat.elems[loc_ge]);
    Mge->PdTrnd = Mgs->PdTrnd;
    lat.SetdT(ge, i);

    for (k = 0; k <= 1; k++) {
      dx_gs[k] = lat.elems[loc_gs]->dS[k]; dx_ge[k] = lat.elems[loc_ge]->dS[k];
    }

    // move elements onto mis-aligned girder
    for (j = loc_gs+1; j < loc_ge; j++) {
      if ((lat.elems[j]->Pkind == tse::Mpole)
	  || (lat.elems[j]->Fnum == lat.conf.bpm)) {
	Mj = dynamic_cast<tse::MpoleType*>(lat.elems[j]);
        s = lat.elems[j]->S;
	for (k = 0; k <= 1; k++)
	  Mj->PdSsys[k] = dx_gs[k] + (dx_ge[k]-dx_gs[k])*(s-s_gs)/(s_ge-s_gs);
	Mj->PdTsys = Mgs->PdTrms*Mgs->PdTrnd;
      }
    }
  }
}


void ts::misalign_sys_girders(tsc::LatticeType &lat, const int gs, const int ge,
			  const double dx, const double dy,
			  const double dr)
{
  int       i, k, n_girders, n_ge, n_gs;
  long int  loc_gs, loc_ge, j;
  double    s_gs, s_ge, dx_gs[2], dx_ge[2], s;
  tse::MpoleType *Mgs, *Mge, *Mj;

  n_gs = lat.GetnKid(gs); n_ge = lat.GetnKid(ge);

  if (n_gs == n_ge)
    n_girders = n_gs;
  else {
    std::cout << "set_girders: no of GS != no of GE" << std::endl;
    exit (1);
  }

  misalign(lat, ts::dx_sys, gs, dx, dy, dr, false);
  misalign(lat, ts::dx_sys, ge, dx, dy, dr, false);

  for (i = 1; i <= n_girders; i++) {
    loc_gs = lat.Elem_GetPos(gs, i); loc_ge = lat.Elem_GetPos(ge, i);
    s_gs = lat.elems[loc_gs]->S; s_ge = lat.elems[loc_ge]->S;

    // roll for a rigid boby
    // Note, girders needs to be introduced as gs->ge pairs
    Mgs = dynamic_cast<tse::MpoleType*>(lat.elems[loc_gs]);
    Mge = dynamic_cast<tse::MpoleType*>(lat.elems[loc_ge]);
    Mge->PdTrnd = Mgs->PdTrnd;
    lat.SetdT(ge, i);

    for (k = 0; k <= 1; k++) {
      dx_gs[k] = lat.elems[loc_gs]->dS[k]; dx_ge[k] = lat.elems[loc_ge]->dS[k];
    }

    // move elements onto mis-aligned girder
    for (j = loc_gs+1; j < loc_ge; j++) {
      if ((lat.elems[j]->Pkind == tse::Mpole)
	  || (lat.elems[j]->Fnum == lat.conf.bpm)) {
	Mj = dynamic_cast<tse::MpoleType*>(lat.elems[j]);
        s = lat.elems[j]->S;
	for (k = 0; k <= 1; k++)
	  Mj->PdSsys[k] = dx_gs[k] + (dx_ge[k]-dx_gs[k])*(s-s_gs)/(s_ge-s_gs);
	Mj->PdTsys = Mgs->PdTrms*Mgs->PdTrnd;
      }
    }
  }
}


// Apertures.

void ts::set_aper(tsc::LatticeType &lat, const int Fnum, const int Knum,
	      const double Dxmin, const double Dxmax,
	      const double Dymin, const double Dymax)
{
  const int k = lat.Elem_GetPos(Fnum, Knum);

  lat.elems[k]->maxampl[X_][0] = Dxmin; lat.elems[k]->maxampl[X_][1] = Dxmax;
  lat.elems[k]->maxampl[Y_][0] = Dymin; lat.elems[k]->maxampl[Y_][1] = Dymax;
}

void ts::set_aper(tsc::LatticeType &lat, const int Fnum,
	      const double Dxmin, const double Dxmax,
	      const double Dymin, const double Dymax)
{
  int k;

  for (k = 1; k <= lat.GetnKid(Fnum); k++)
    set_aper(lat, Fnum, k, Dxmin, Dxmax, Dymin, Dymax);
}

void ts::set_aper_type(tsc::LatticeType &lat, const int type, const double Dxmin,
		   const double Dxmax, const double Dymin, const double Dymax)
{
  long int  k;
  tse::MpoleType *M;

  if (type >= tse::All && type <= HOMmax) {
    for(k = 1; k <= lat.conf.Cell_nLoc; k++) {
      M = dynamic_cast<tse::MpoleType*>(lat.elems[k]);
      if (((lat.elems[k]->Pkind == tse::Mpole) && (M->n_design == type))
	  || (type == tse::All))
	set_aper(lat, lat.elems[k]->Fnum, lat.elems[k]->Knum, Dxmin, Dxmax,
		 Dymin, Dymax);
    }
  } else
    printf("set_aper_type: bad design type %d\n", type);
}


// Miscellaneous.

double ts::get_L(tsc::LatticeType &lat, const int Fnum, const int Knum)
{
  return lat.elems[lat.Elem_GetPos(Fnum, Knum)]->PL;
}


void ts::set_L(tsc::LatticeType &lat, const int Fnum, const int Knum, const double L)
{
  long int  loc;
  double    phi;
  tse::ElemType  *elemp;
  tse::MpoleType *M;

  loc = lat.Elem_GetPos(Fnum, Knum);
  elemp = lat.elems[loc];
  if (elemp->Pkind == tse::Mpole) {
    M = dynamic_cast<tse::MpoleType*>(elemp);
    if (M->Pirho != 0e0) {
      // Phi is constant.
      phi = elemp->PL*M->Pirho; M->Pirho = phi/L;
      // M->Pc0 = sin(phi/2e0);
    }
  }
  elemp->PL = L;
}


void ts::set_L(tsc::LatticeType &lat, const int Fnum, const double L)
{
  int k;

  for (k = 1; k <= lat.GetnKid(Fnum); k++)
    set_L(lat, Fnum, k, L);
}


void ts::set_dL(tsc::LatticeType &lat, const int Fnum, const int Knum, const double dL)
{
  lat.elems[lat.Elem_GetPos(Fnum, Knum)]->PL += dL;
}


void ts::set_dL(tsc::LatticeType &lat, const int Fnum, const double dL)
{
  int k;

  for (k = 1; k <= lat.GetnKid(Fnum); k++)
    set_dL(lat, Fnum, k, dL);
}


// Multipoles.

void ts::get_bn(tsc::LatticeType &lat, const ts::bn_types bn_type, const int Fnum,
	    const int Knum, const int n, double &bn, double &an)
{
  tse::ElemType  *elem;
  tse::MpoleType *M;

  if (n < 1) {
    std::cout << "get_bnL_design_elem: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  elem = lat.elems[lat.Elem_GetPos(Fnum, Knum)];
  M    = dynamic_cast<tse::MpoleType*>(elem);

  switch (bn_type) {
  case ts::bn_des:
    bn = M->PBpar[HOMmax+n];
    an = M->PBpar[HOMmax-n];
    break;
  case ts::bnL_des:
    bn = (elem->PL != 0e0)? M->PBpar[HOMmax+n]*elem->PL : M->PBpar[HOMmax+n];
    an = (elem->PL != 0e0)? M->PBpar[HOMmax-n]*elem->PL : M->PBpar[HOMmax-n];
    break;
  default:
    printf("\nget_bn: undefined error type: %d\n", bn_type);
    exit(1);
    break;
  }
}


void ts::set_bn(tsc::LatticeType &lat, const ts::bn_types bn_type, const int Fnum,
	    const int Knum, const int n, const double bn, const double an,
	    const bool new_rnd)
{
  int       nd;
  tse::ElemType* elem = lat.elems[lat.Elem_GetPos(Fnum, Knum)];
  tse::MpoleType *M   = dynamic_cast<tse::MpoleType*>(elem);

  if (n < 1) {
    std::cout << "set_bn: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  switch (bn_type) {
  case ts::bn_des:
    M->PBpar[HOMmax+n] = bn;
    M->PBpar[HOMmax-n] = an;
    break;
  case  ts::dbn_des:
    M->PBpar[HOMmax+n] += bn;
    M->PBpar[HOMmax-n] += an;
    break;
  case  ts::bnL_des:
    M->PBpar[HOMmax+n] = (elem->PL != 0e0)? bn/elem->PL : bn;
    M->PBpar[HOMmax-n] = (elem->PL != 0e0)? an/elem->PL : an;
    break;
  case  ts::dbnL_des:
    M->PBpar[HOMmax+n] += (elem->PL != 0e0)? bn/elem->PL : bn;
    M->PBpar[HOMmax-n] += (elem->PL != 0e0)? an/elem->PL : an;
    break;
  case  ts::bnL_sys:
    M->PBsys[HOMmax+n] += (elem->PL != 0e0)? bn/elem->PL : bn;
    M->PBsys[HOMmax-n] += (elem->PL != 0e0)? an/elem->PL : an;
    break;
  case  ts::bnL_rms:
    M->PBrms[HOMmax+n] += (elem->PL != 0e0)? bn/elem->PL : bn;
    M->PBrms[HOMmax-n] += (elem->PL != 0e0)? an/elem->PL : an;
    if (new_rnd) {
	M->PBrnd[HOMmax+n] = (normal)? normranf() : ranf();
	M->PBrnd[HOMmax-n] = (normal)? normranf() : ranf();
      }
    break;
  case  ts::bnr_sys:
    // Errors are relative to design values for: [Dip, Quad, Sext, ...].
    nd = M->n_design;
    M->PBsys[HOMmax+n] = bn*M->PBpar[HOMmax+nd];
    M->PBsys[HOMmax-n] = an*M->PBpar[HOMmax+nd];
    break;
  case  ts::bnr_rms:
    // Errors are relative to design values for: [Dip, Quad, Sext, ...].
    nd = M->n_design;
    if (nd == tse::Dip) {
      M->PBrms[HOMmax+n] = bn*M->Pirho;
      M->PBrms[HOMmax-n] = an*M->Pirho;
    } else {
      M->PBrms[HOMmax+n] = bn*M->PBpar[HOMmax+nd];
      M->PBrms[HOMmax-n] = an*M->PBpar[HOMmax+nd];
    }
    if(new_rnd){
	M->PBrnd[HOMmax+n] = (normal)? normranf() : ranf();
	M->PBrnd[HOMmax-n] = (normal)? normranf() : ranf();
    }
    break;
  default:
    printf("\nset_bn: undefined error type: %d\n", bn_type);
    exit(1);
    break;
  }

  lat.SetPB(Fnum, Knum, n); lat.SetPB(Fnum, Knum, -n);
}


void ts::set_bn(tsc::LatticeType &lat, const ts::bn_types bn_type, const int Fnum,
	    const int n, const double bn, const double an, const bool new_rnd)
{
  int k;

  if (n < 1) {
    std::cout << "set_bn: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  for (k = 1; k <= lat.GetnKid(Fnum); k++)
    set_bn(lat, bn_type, Fnum, k, n, bn, an, new_rnd);
}


//------------------------------------------------------------------------------

#if 0

void ts::set_bnL_design_type(tsc::LatticeType &lat, const int type, const int n,
			 const double bnL, const double anL)
{
  long int  k;
  tse::MpoleType *M;

  if (n < 1) {
    std::cout << "set_bnL_design_type: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  if ((type >= Dip) && (type <= HOMmax)) {
    for (k = 1; k <= lat.conf.Cell_nLoc; k++) {
      M = dynamic_cast<tse::MpoleType*>(lat.elems[k]);
      if ((lat.elems[k]->Pkind == tse::Mpole) && (M->n_design == type))
	set_bnL_design_elem(lat, lat.elems[k]->Fnum, lat.elems[k]->Knum, n, bnL,
			    anL);
    }
  } else
    printf("Bad type argument to set_bnL_design_type()\n");
}


void  ts::set_bnL_sys_type(tsc::LatticeType &lat, const int type, const int n,
		      const double bnL, const double anL)
{
  long int  k;
  tse::MpoleType *M;

  if (n < 1) {
    std::cout << "set_bnL_sys_type: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  if (type >= Dip && type <= HOMmax) {
    for(k = 1; k <= lat.conf.Cell_nLoc; k++) {
      M = dynamic_cast<tse::MpoleType*>(lat.elems[k]);
      if ((lat.elems[k]->Pkind == tse::Mpole) && (M->n_design == type))
	set_bnL_sys_elem(lat, lat.elems[k]->Fnum, lat.elems[k]->Knum, n, bnL,
			 anL);
    }
  } else
    printf("Bad type argument to set_bnL_sys_type()\n");
}


void  ts::set_bnL_rms_type(tsc::LatticeType &lat, const int type, const int n,
		      const double bnL, const double anL, const bool new_rnd)
{
  long int  k;
  tse::MpoleType *M;

  if (n < 1) {
    std::cout << "get_bnL_rms_type: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  if (type >= Dip && type <= HOMmax) {
    for(k = 1; k <= lat.conf.Cell_nLoc; k++) {
      M = dynamic_cast<tse::MpoleType*>(lat.elems[k]);
      if ((lat.elems[k]->Pkind == tse::Mpole) && (M->n_design == type))
	set_bnL_rms_elem(lat, lat.elems[k]->Fnum, lat.elems[k]->Knum, n, bnL,
			 anL, new_rnd);
    }
  } else
    printf("Bad type argument to set_bnL_rms_type()\n");
}


void  ts::set_bnr_sys_type(tsc::LatticeType &lat, const int type, const int n,
		      const double bnr, const double anr)
{
  long int  k;
  tse::MpoleType *M;

  if (n < 1) {
    std::cout << "set_bnr_sys_type: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  if (type >= Dip && type <= HOMmax) {
    for(k = 1; k <= lat.conf.Cell_nLoc; k++) {
      M = dynamic_cast<tse::MpoleType*>(lat.elems[k]);
      if ((lat.elems[k]->Pkind == tse::Mpole) && (M->n_design == type))
	set_bnr_sys_elem(lat, lat.elems[k]->Fnum, lat.elems[k]->Knum, n, bnr,
			 anr);
    }
  } else
    printf("Bad type argument to set_bnr_sys_type()\n");
}


void  ts::set_bnr_rms_type(tsc::LatticeType &lat, const int type, const int n,
		      const double bnr, const double anr, const bool new_rnd)
{
  long int  k;
  tse::MpoleType *M;

  if (n < 1) {
    std::cout << "set_bnr_rms_type: n < 1 (" << n << ")" << std::endl;
    exit(1);
  }

  if (type >= Dip && type <= HOMmax) {
    for(k = 1; k <= lat.conf.Cell_nLoc; k++) {
      M = dynamic_cast<tse::MpoleType*>(lat.elems[k]);
      if ((lat.elems[k]->Pkind == tse::Mpole) && (M->n_design == type))
	set_bnr_rms_elem(lat, lat.elems[k]->Fnum, lat.elems[k]->Knum, n, bnr,
			 anr, new_rnd);
    }
  } else
    printf("Bad type argument to set_bnr_rms_type()\n");
}

#endif
