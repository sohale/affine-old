#include <iostream>
#include <cassert>
#include <cmath>
#include <stdexcept>

#include <random> //for tests only

/*
Author: Sohail Siadatnejad

10 Jul 2015, 13:06

Compile using:
    g++ -std=c++11 affine.cpp
To run the tests:
    ./a.out
*/


//template <typename R> //todog++ -std=c++11 affine.cpp

namespace affine{

typedef double R;
bool almost_zero(R x)
{
    constexpr R TOLERANCE = 0.0000001;
    assert(TOLERANCE>0);
    return abs(x)<TOLERANCE;
}



/* A minimal but performant vector function using minimal accessors.
*/
class Vec2D
{
public:
    typedef unsigned short index_type;

private:
    R x[2];

public:
    const R& operator[](index_type i) const
    {
        if(i<0) throw std::invalid_argument( "index out of range" );
        if(i>1) throw std::invalid_argument( "index out of range" );
        return x[i];
    }

    static constexpr index_type dims()
    /* Keep for future extention to Vec3D, etc*/
    {
        return 2;
    }
    
    /*todo: Use of == is discouraged. Asserts may not work for real/float types. */
    bool operator==(const Vec2D& b) const
    {
        return (this->x[0] == b.x[0]) && (this->x[1] == b.x[1]);
    }
    Vec2D operator+(const Vec2D& b) const
    {
        return Vec2D( this->x[0] + b.x[0], this->x[1] + b.x[1] );
    }
    Vec2D operator-(const Vec2D& b) const
    {
        return Vec2D( this->x[0] - b.x[0], this->x[1] - b.x[1] );
    }
    Vec2D operator-() const
    {
        return Vec2D( -this->x[0], -this->x[1] );
    }

    constexpr Vec2D(R _x, R _y)
        : x{_x,_y}
    {
    }

    std::string to_string() const
    {
        //todo: dont show decimals
        return std::to_string(x[0])+","+std::to_string(x[1]);
    }
    R norm2_L2() const
    {
        return x[0]*x[0]+x[1]*x[1];
    }
    R inner(const Vec2D& b) const
    {
        return x[0]*b.x[0]+x[1]*b.x[1];
    }

    //a constexpr would have been preferred.
    static const Vec2D e(const index_type i){
        return Vec2D(i==0?1:0,i==1?1:0);
    }

    static const Vec2D& o;

    /*
        Notes:
        * Don't define a default constructor
        * The operator* was not used because it is ambiguous: inner, scalar, outer, etc
        Todo:
        * Replace the ::operator"==" with ::almost_equal()
    */
};

const Vec2D& Vec2D::o = Vec2D(0,0);

Vec2D operator "" _x (long double x)
{
    return Vec2D(x,0.0);
}
Vec2D operator "" _x (unsigned long long const x)
{
    return Vec2D(x,0.0);
}

Vec2D operator "" _y (long double y)
{
    return Vec2D(0.0,y);
}
Vec2D operator "" _y (unsigned long long y)
{
    return Vec2D(0.0,y);
}

std::ostream& operator<<(std::ostream& os, const Vec2D& v)
{
    os << v.to_string();
    return os;
}

/* end of class Vect2D  */

void test_Vec2D(){
    Vec2D o=Vec2D(0,0);
    std::cout << o.to_string() << std::endl;

    Vec2D a=Vec2D(1.0, 0);

    assert( o==o  );
    assert( a==a );
    assert( ! (a==o) );

    assert( a+o==o+a );
    assert( a - o==a );
    assert( o + o==o );

    auto b=a;  //Note that we dont a copy/move constructor or default constructor

    static_assert(std::is_copy_constructible<Vec2D>::value,
                  "copying required");
    static_assert(std::is_nothrow_move_constructible<Vec2D>::value
               && std::is_nothrow_move_assignable<Vec2D>::value,
                  "may throw");
    //http://en.cppreference.com/w/cpp/language/static_assert


    std::cout << (o+a).to_string() << std::endl;
    std::cout << (a+a).to_string() << std::endl;

    std::cout << a.to_string() << std::endl;

    std::cout << (3.0_x + 4.0_y).to_string() << std::endl;

    /*
    //copy (by-val) test not done. ==> Vec2D is immutable.
    auto c=b;
    c.x[0] = 5;
    assert(b.x[0]==1.0);
    */
    
    Vec2D c = 5_x+0_y;
    std::cout <<"b= "<< b.to_string() << std::endl;
    std::cout <<"c= "<< c.to_string() << std::endl;


    std::cout << b.norm2_L2() << std::endl;
    std::cout << c.norm2_L2() << std::endl;

    assert( a.inner(a)==a.norm2_L2() );
    assert( c.inner(c)==c.norm2_L2() );

    Vec2D d = Vec2D::o;
    std::cout << d[0] << " * ";
    for(int i=0;i<2;i++)
        std::cout << d[i] << " *** ";
    std::cout << std::endl;

    std::cout << d << " ";
    std::cout << d << " ";
    std::cout << d << " ";
    std::cout << std::endl;

    //The following cause error. (Which is good)
    //Vec2D f;
    //Vec2D g();
    //Vec2D h=Vec2D();
    //    d[0]=6.0; //lvalue required as left operand of assignment

    double u=2;
    double& z = u; //d[0];
    std::cout << u << "=u ";
    z=8;
    std::cout << u << "=u ";
    std::cout << std::endl;


    Vec2D v=Vec2D::o;
    v = 1.0_x;

    //todo: more tests
};


}
