#include <tps/tps_type.h>
#include <tps/ss_vect.h>
#include <cmath>
#include <iostream>

int main(int argc, char *argv[])
{

    auto a = tps(42e0);
    std::cout << "a 42\n" << a << std::endl;

    auto x = tps(0, 0+1);
    auto y = tps(0, 2+1);

    std::cout << "x\n" << x << std::endl;
    std::cout << "y\n" << y << "\n" << std::endl;

    std::cout << "a + x\n" << a + x << std::endl;
    std::cout << "a + y\n" << a + y << "\n" << std::endl;

    std::cout << "a * x\n" << a * x << std::endl;
    std::cout << "a * y\n" << a * y << "\n" << std::endl;

#if 0
    //t1 = 3e0;
    //t2 = 5e0;

    std::cout << "t1\n" << t1 << std::endl
	      << "t2\n" << t2 << std::endl;

    std::cout << "\nt1 * t2 * 42\n" << t1 * t2 + 42<< std::endl;
    std::cout << "t1 + t2 \n" << t1 + t2 << std::endl;


    // t1 = M_PI/6;
    std::cout << "sin(t1) \n" << sin(t1) << std::endl;
    std::cout << "sin(t2) \n" << sin(t2) << std::endl;

    std::cout << "\ncos(t1) \n" << cos(t1) << std::endl;
    std::cout << "cos(t2) \n" << cos(t2) << std::endl;

    ss_vect<tps> v1, v2;
    v1.set_identity();
    v2.set_identity();

    ss_vect<tps> v3 = v1 * v2;

    std::cout << "v1\n" << v1 << "v2\n" << v2
	      << "v1 * v2\n"<< v3 <<std::endl;

#endif
}
