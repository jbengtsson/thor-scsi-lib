#include <gtpsa/ss_vect.h>
#include <gtpsa/lielib.hpp>
#include <thor_scsi/elements/sextupole.h>
#include <iostream>
#include <assert.h>

namespace tsc = thor_scsi::core;
namespace tse = thor_scsi::elements;

#if 0

void test1(void)
{
  const int
    nv = 6, no = 3, np = 1, po = 2;

  auto desc = std::make_shared<gtpsa::desc>(nv, no, np, po);

  auto M1   = gtpsa::ss_vect<gtpsa::tpsa>(desc, no);

  M1.set_zero();
  for (auto k = 0; k < nv; k++)
    M1[k].set(k+1, 0e0, k+1e0);
  // (start location, cst, ,,,).
  M1[0].setv
    (0,
     {-1e0,
      1e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0,
      0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0,
      0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0,
      3.14e-3});
  std::cout << "\nM1:" << M1;

  auto desc2 = std::make_shared<gtpsa::desc>(nv+np, no);

  gtpsa::ss_vect<gtpsa::tpsa> M2(desc2, no);

  M1[0].print("\nM1[0]:", 1e-30, 0);
  M2 = param_to_ss_vect(M1);
  print_map("\nM2:", M2);

  // Validate inverse transformation.
  ss_vect_to_param(M2)[0].print("\nM1[0]:", 1e-30, 0);

  auto Id = gtpsa::ss_vect<gtpsa::tpsa>(desc2, no);
  Id.set_identity();
  for (auto k = 0; k < Id.size(); k++) 
    Id[k] = 2e0*Id[k];
  std::cout << "\nId:" << Id;
  Id[0].print("\nId[0]:", 1e-30, 0);

  compose(M2, Id)[0].print("\nM2*Id[0]:", 1e-30, 0);
  compose(Id, M2)[0].print("\nId*M2[0]:", 1e-30, 0);
  compose(M2, M2)[0].print("\nM2*M2[0]:", 1e-30, 0);
}

#endif

void test2(void)
{
  const int mo = 3, nv = 7;

  auto desc = std::make_shared<gtpsa::desc>(nv, mo);
  Config C;
  C.set<std::string>("name", "test");
  C.set<double>("K", 3e0);
  C.set<double>("L", 0e0);
  C.set<double>("N", 1);

  tsc::ConfigType calc_config;
  auto sext = tse::SextupoleType(C);

  auto M = gtpsa::ss_vect<gtpsa::tpsa>(desc, mo);
  M.set_identity();
  sext.propagate(calc_config, M);

  auto M_inv = gtpsa::ss_vect<gtpsa::tpsa>(desc, mo);
  M_inv = gtpsa::minv(M);
  std::cout << "\nM:\n" << M << "\n";
  print_map("\nM:", M);
  std::cout << "\nM^-1:\n" << M_inv << "\n";
  print_map("\nM^-1:", M_inv);

  auto M_M_inv = gtpsa::ss_vect<gtpsa::tpsa>(desc, mo);
  M_M_inv = gtpsa::compose(M, M_inv);
  std::cout << "\nM_M_inv:\n" << M_M_inv << "\n";
  print_map("\nM_M_inv:", M_M_inv);
  assert(false);

  auto h = gtpsa::M_to_h_DF(M);
  h.print("\nh:", 1e-30, 0);
}


gtpsa::ss_vect<gtpsa::tpsa> compute_sext_map
(const int no, const std::shared_ptr<gtpsa::mad::desc> &desc)
{
  Config C;
  C.set<std::string>("name", "test");
  C.set<double>("K", 3e0);
  C.set<double>("L", 0e0);
  C.set<double>("N", 1);

  tsc::ConfigType config;
  auto sext = tse::SextupoleType(C);
  auto M    = gtpsa::ss_vect<gtpsa::tpsa>(desc, no);

  M.set_identity();
  sext.propagate(config, M);
  print_map("\nM:", M);

  return M;
}


void test3(void)
{
  const int no = 3, nv = 6 + 1;

  auto desc = std::make_shared<gtpsa::desc>(nv, no);
  auto M    = gtpsa::ss_vect<gtpsa::tpsa>(desc, no);

  const int nm = desc->maxLen(no);

  M = compute_sext_map(no, desc);

  std::cout << "\nM[1].get(8)  = " << std::scientific << std::setprecision(3)
	    << M[1].get(8) << "\n";

  std::cout << "M[1].get()   = " << std::scientific << std::setprecision(3)
	    <<  M[1].get(std::vector<ord_t>{2, 0, 0, 0, 0, 0, 0}) << "\n";

  std::cout << "M[1].index() = "
	    <<  M[1].index(std::vector<ord_t>{2, 0, 0, 0, 0, 0, 0}) << "\n";

  ord_t              err;
  std::vector<ord_t> ind(nv);
  err = M[4].mono(11, &ind);
  std::cout << "M[4].mono(11):  \n  " << (int)err << "\n ";
  for (auto k: ind)
    std::cout << std::setw(2) << (int)k;
  std::cout << "\n";

  if (false) {
    std::vector<num_t> v(nm);
    M[0].getv(0, &v);
    print_vec("\nv:", v);
  }
}


void dragt_finn_fact(void)
{
  const int no = 3, nv = 6 + 1;

  auto desc = std::make_shared<gtpsa::desc>(nv, no);
  auto M    = gtpsa::ss_vect<gtpsa::tpsa>(desc, no);

  M = compute_sext_map(no, desc);

  auto h = gtpsa::M_to_h_DF(M);
  h.print("\nh:", 1e-30, 0);
}


void setvar(void)
{
  auto a_desc = std::make_shared<gtpsa::desc>(6, 1);
  auto a      = gtpsa::tpsa(a_desc, mad_tpsa_default);

  a.setVariable(0e0, 2, 0e0);
  a.print("\na:");
}


int main(int argc, char *argv[])
{
  // if (false)
  //   test1();

  // if (false)
  //   test2();

  if (false)
    setvar();

  if (false)
    test3();

  if (!false)
    dragt_finn_fact();
}
