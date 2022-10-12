#include <tps/ss_vect.h>
#include <iostream>

int main(int argc, char *argv[])
{

    ss_vect<tps> ss;
    std::cout << "state space vector start \n"
	      << ss.cst() << "\n"
	      << ss << std::endl;

    ss.set_identity();
    std::cout << "state space vector  identity\n"
	      << ss.cst() << "\n"
	      << ss << std::endl;

}
