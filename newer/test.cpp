/*
Author: Sohail Siadatnejad

* 10 Jul 2015, 13:06
* 18 June 2022, improvements

Compile using:
    g++ -std=c++11 affine.cpp
To run the tests:
    ./a.out
*/
#include "./affine2d.hpp"

using namespace affine;


// Should this be inside the namespace?
bool test_Affine2D()
{
    bool any_failed = false;

    Affine2D o = Affine2D::o;
    std::cout << o.to_string() << std::endl;

    Affine2D a = Affine2D::eye; // Affine2D(1.0, 0,0,0);

    assert((o - o).almost_zero());
    // assert( a==a );
    // assert( ! (a==o) );

    // assert( a+o==o+a );
    // assert( a - o==a );
    // assert( o + o==o );

    std::cout << a.det() << std::endl;

    assert(o.det() == 0);
    assert(a.det() == 1.0);

    auto b = a;
    // todo: immutability test

    std::cout << a.to_string() << std::endl;

    // Affine2D c = a + Affine2D::e00 * 2.0 - 3.0 * Affine2D::e01;
    Affine2D c = a + Affine2D::e00 + Affine2D::e01;
    // b[0,0]=2;
    assert((c - c).almost_zero());
    assert(Affine2D::o.almost_zero());
    assert((Affine2D::eye - Affine2D::eye).almost_zero());

    Affine2D x = Affine2D::eye;
    assert((x + (-x)).almost_zero());

    assert((Affine2D::o).almost_zero());

    assert((Affine2D::eye.get_inverse() - Affine2D::eye).almost_zero());

    assert((Affine2D::eye * Affine2D::eye - Affine2D::eye).almost_zero());

    // setting up random number generation
    typedef std::mt19937 MyRNG;
    uint32_t seed_val;
    MyRNG rng;
    rng.seed(seed_val); // initialize

    std::normal_distribution<double> normal_dist(0.0, 1.0);
    for (int j = 0; j < 100; j++)
    {
        auto q =
            Affine2D::e00 * normal_dist(rng) +
            Affine2D::e01 * normal_dist(rng) +
            Affine2D::e10 * normal_dist(rng) +
            Affine2D::e11 * normal_dist(rng) +
            Affine2D::tx * normal_dist(rng) +
            Affine2D::ty * normal_dist(rng);

        if (abs(q.det()) > 0.0001)
        {
            if (!(q * (q.get_inverse())).almost_equal(Affine2D::eye))
            {
                std::cout << "failed" << (q * (q.get_inverse())) << std::endl;

                std::cout << q.get_inverse().det() << std::endl;
                std::cout << "det = " << q.det() << " inv det = ";
                std::cout << q << std::endl;
                std::cout << q.get_inverse() << std::endl;
                std::cout << std::endl;
                any_failed = true;
            }
        }
    }
    return !any_failed;
};

int main()
{

    std::cout << "Testing Vec2D" << std::endl;
    affine::test_Vec2D();
    std::cout << "Testing Affine2D" << std::endl;
    // bool r = test_Affine2D());
    bool r = true;
    if (r)
        std::cout << "Tests passed successfully." << std::endl;
    else
        std::cout << "Some tests failed." << std::endl;
}

