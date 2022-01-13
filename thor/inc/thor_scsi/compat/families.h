#ifndef _THOR_SCSI_LEGACY_FAMILIES_H_
#define _THOR_SCSI_LEGACY_FAMILIES_H_ 1
#warning "Functions in these headers require to be implemented"

void set_bn_design_fam(const int Fnum,
		       const int n, const double bn, const double an);

void set_dbnL_design_fam(const int Fnum,
			 const int n, const double dbnL, const double danL);

void get_bn_design_elem(const int Fnum, const int Knum,
			const int n, double &bn, double &an);

void set_dbnL_design_elem(const int Fnum, const int Knum,
			  const int n, const double dbnL, const double danL);

#endif /* _THOR_SCSI_LEGACY_FAMILIES_H_ */
/*
 * Local Variables:
 * mode: c++
 * c-file-style: "python"
 * End:
 */
