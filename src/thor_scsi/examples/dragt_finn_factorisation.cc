#include <gtpsa/ss_vect.h>
#include <gtpsa/lielib.hpp>
#include <thor_scsi/elements/sextupole.h>
#include <iostream>
#include <assert.h>

namespace tsc = thor_scsi::core;
namespace tse = thor_scsi::elements;

void print_map(const std::string &str, gtpsa::ss_vect<gtpsa::tpsa> &M)
{
  for (auto k = 0; k < M.size(); k++) 
    M[k].print(str.c_str(), 1e-30, 0);
}


void print_vec(std::vector<num_t> &v)
{
  for (auto mn: v)
    std::cout << std::scientific << std::setprecision(3)
	      << std::setw(11) << mn;
  std::cout << "\n";
}


gtpsa::ss_vect<gtpsa::tpsa> param_to_ss_vect(gtpsa::ss_vect<gtpsa::tpsa> &A)
{
  const auto desc0 = A[0].getDescription();
  const auto info  = desc0->getInfo();
  const auto nv    = info.getNumberOfVariables();
  const auto no    = info.getVariablesMaximumOrder();
  const auto np    = info.getNumberOfParameters();
  const auto po    = info.getParametersMaximumOrder();

  auto desc1 = std::make_shared<gtpsa::desc>(nv+np, no);
  auto nm    = desc0->maxLen(no);

  gtpsa::ss_vect<gtpsa::tpsa> B(desc1, no);

  for (auto k = 0; k < A.size(); k++) {
    std::vector<num_t> v(nm);

    A[k].getv(0, &v);
    B[k].setv(0, v);
  }

  return B;
}


gtpsa::ss_vect<gtpsa::tpsa> ss_vect_to_param(gtpsa::ss_vect<gtpsa::tpsa> &A)
{
  const int ps_dim = 6;

  const auto desc0 = A[0].getDescription();
  const auto info  = desc0->getInfo();
  const auto nv    = info.getNumberOfVariables();
  const auto no    = info.getVariablesMaximumOrder();
  const auto np    = info.getNumberOfParameters();
  const auto po    = info.getParametersMaximumOrder();

  auto desc1 = std::make_shared<gtpsa::desc>(ps_dim, no, nv-ps_dim, no);
  auto nm    = desc0->maxLen(no);

  gtpsa::ss_vect<gtpsa::tpsa> B(desc1, no);

  for (auto k = 0; k < A.size(); k++) {
    std::vector<num_t> v(nm);

    A[k].getv(0, &v);
    B[k].setv(0, v);
  }

  return B;
}


void test(void)
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


void dragt_finn_fact(void)
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


int main(int argc, char *argv[])
{
  if (false)
    test();

  if (!false)
    dragt_finn_fact();
}
