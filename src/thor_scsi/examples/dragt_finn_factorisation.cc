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


void dragt_finn_fact()
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


void test()
{
  const int
    mo = 4, nv = 6, np = 1, po = 4;

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
  M2 = compose(M1, Id);
  M2[0].print("\nM2[0]:", 1e-30, 0);
  M2 = compose(Id, M1);
  M2[0].print("\nM2[0]:", 1e-30, 0);
  M2 = compose(M1, M1);
  M2[0].print("\nM2[0]:", 1e-30, 0);
}


int main(int argc, char *argv[])
{
  // dragt_finn_fact();

  test();
}
