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

  const auto desc = std::make_shared<gtpsa::desc>(nv, no, np, po);

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

  const auto desc2 = std::make_shared<gtpsa::desc>(nv+np, no);

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

  const auto desc = std::make_shared<gtpsa::desc>(nv, mo);
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
(const std::shared_ptr<gtpsa::mad::desc> &desc, const int no)
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

  const auto desc = std::make_shared<gtpsa::desc>(nv, no);

  auto M    = gtpsa::ss_vect<gtpsa::tpsa>(desc, no);

  const int nm = desc->maxLen(no);

  M = compute_sext_map(desc, no);

  std::cout << "\nM[1].get(8)  = " << std::scientific << std::setprecision(3)
	    << M[1].get(8) << "\n";

  std::cout << "M[1].get()   = " << std::scientific << std::setprecision(3)
	    <<  M[1].get(std::vector<ord_t>{2, 0, 0, 0, 0, 0, 0}) << "\n";

  std::cout << "M[1].index() = "
	    <<  M[1].index(std::vector<ord_t>{2, 0, 0, 0, 0, 0, 0}) << "\n";

  M[4].print();

  auto n = 11;
  std::vector<ord_t> ind(nv);
  auto ord = M[4].mono(n, &ind);
  std::cout << std::scientific << std::setprecision(16)
	    << "\nM[4]:\n" << std::setw(6) << n
	    << std::setw(25) << M[4].get(ind)
	    << std::setw(5) << (int)ord << "   ";
  for (auto k = 0; k < ind.size(); k++)
    std::cout << ((k % 2 == 1)? std::setw(2) : std::setw(3)) << (int)ind[k];
  std::cout << "\n";

  if (false) {
    std::cout << "\n";
    std::vector<num_t> v(nm);
    M[0].getv(0, &v);
    print_vec("\nv:", v);
  }
}


void test4(void)
{
  const int no = 3, nv = 6 + 1;

  const auto desc = std::make_shared<gtpsa::desc>(nv, no);

  auto x    = gtpsa::tpsa(desc, no);

  x.setv(0, {-0.123, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0});
  x.print();

  std::vector<ord_t> exps = {0, 0, 0, 1, 0, 0, 0};

  auto ind = x.index(exps);
  printf("\n%2d\n", ind);
  printf("%10.3e\n", x.get(ind));
  printf("%10.3e\n", x.get(exps));
}


void dragt_finn_fact(void)
{
  const int no = 3, nv = 6 + 1;

  const auto desc = std::make_shared<gtpsa::desc>(nv, no);

  auto M    = gtpsa::ss_vect<gtpsa::tpsa>(desc, no);

  M = compute_sext_map(desc, no);

  auto h = gtpsa::M_to_h_DF(M);
  h.print("\nh:", 1e-30, 0);
}


void setvar(void)
{
  const auto a_desc = std::make_shared<gtpsa::desc>(6, 1);

  auto a      = gtpsa::tpsa(a_desc, mad_tpsa_default);

  // Parameters: (constant part, monomial index, value).
  a.setVariable(3.1415, 3, 2.0);
  a.print("\na:");
}

#define RET_VALUE 0
#if RET_VALUE
void test_map_norm(void)
{
  const int no = 3, nv = 6 + 1;

  const auto desc = std::make_shared<gtpsa::desc>(nv, no);

  auto M = gtpsa::ss_vect<gtpsa::tpsa>(desc, no);

  MNFType MNF(desc, no);

  M = compute_sext_map(desc, no);

  Map_Norm(M);
  // MNF = Map_Norm(M);
}
#endif


int main(int argc, char *argv[])
{
  // if (false)
  //   test1();

  // if (false)
  //   test2();

  if (false)
    test3();

  if (false)
    test4();

  if (false)
    setvar();

  if (false)
    dragt_finn_fact();

#if RET_VALUE
  if (!false)
    test_map_norm();
#endif
}
