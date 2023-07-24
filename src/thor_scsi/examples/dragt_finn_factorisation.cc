#include <gtpsa/ss_vect.h>
#include <gtpsa/lielib.hpp>
#include <thor_scsi/elements/sextupole.h>
#include <iostream>
#include <assert.h>

namespace tsc = thor_scsi::core;
namespace tse = thor_scsi::elements;

void print_map(const std::string &str, gtpsa::ss_vect<gtpsa::tpsa> &M)
{
  std::cout << str << "\n";
  for (auto k = 0; k < M.size(); k++) 
    M[k].print(str.c_str(), 1e-30, 0);
}


void dragt_finn_fact(void)
{
  const int mo = 4, nv = 6;

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

  std::cout << "\nAfter sextupole:\n" << M << std::endl;

  print_map("\nM:", M);

  auto h = gtpsa::M_to_h_DF(M);
  h.print("\nh:", 1e-30, 0);
}


gtpsa::ss_vect<gtpsa::tpsa> cct
(const gtpsa::ss_vect<gtpsa::tpsa> &A, const gtpsa::ss_vect<gtpsa::tpsa> &B)
{
  const auto desc = A[0].getDescription();
  const auto info  = desc->getInfo();

  const auto nv = info.getNumberOfVariables();
  const auto no = info.getVariablesMaximumOrder();
  const auto np = info.getNumberOfParameters();
  const auto po = info.getParametersMaximumOrder();

  auto nm = A[0].getDescription()->getNv(0,0,0);
  std::vector<num_t> v(nm);
  A[0].getv(0, &v);

  std::cout << info << "\n";
  std::cout << "\nnm = " << nm << "\n";
  for (num_t mn : v)
    std::cout << std::scientific << std::setprecision(3)
	      << std::setw(11) << mn;
  std::cout << "\n";

  auto C = gtpsa::ss_vect<gtpsa::tpsa>(desc, no);
  return C;
}


void test(void)
{
  const int
    mo = 4, nv = 7, np = 1, po = 4;

  auto desc = std::make_shared<gtpsa::desc>(nv, mo, np, po);
  auto Id   = gtpsa::ss_vect<gtpsa::tpsa>(desc, mo);
  auto M1   = gtpsa::ss_vect<gtpsa::tpsa>(desc, mo);
  auto M2   = gtpsa::ss_vect<gtpsa::tpsa>(desc, mo);

  Id.set_identity();
  std::cout << "\nId:" << Id;

  M1.set_zero();
  for (auto k = 0; k < M1.size(); k++)
    M1[k].set(k+1, 0e0, k+1e0);
  // (start location, cst, ,,,).
  M1[0].setv
    (0,
     {0e0,
      1e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0,
      0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0,
      0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0,
      3.14});
  std::cout << "\nM1:" << M1;
  M1[0].print("\nM1[0]:", 1e-30, 0);
  M2 = cct(M1, Id);
  assert(false);
  M2[0].print("\nM2[0]:", 1e-30, 0);
  M2 = compose_jb(Id, M1);
  M2[0].print("\nM2[0]:", 1e-30, 0);
  M2 = compose_jb(M1, M1);
  M2[0].print("\nM2[0]:", 1e-30, 0);
}


int main(int argc, char *argv[])
{
  // dragt_finn_fact();

  test();
}
