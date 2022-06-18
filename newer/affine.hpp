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


//todo: Template<typename R,int size>
class Affine2D {
    /*
    Notes:
        * no default constructor
        * For performance in the 2D case, Matrix type was not used.
    */
    protected:
        R m[4];
        R d[2];
        bool _isNan;
        const unsigned short i00=0;
        const unsigned short i01=1;
        const unsigned short i10=2;
        const unsigned short i11=3;
    private:
        Affine2D(R m00, R m01, R m10, R m11, R dx, R dy, bool kk=false)
        {
            m[i00]=m00;
            m[i01]=m01;
            m[i10]=m10;
            m[i11]=m11;
            d[0]=dx;
            d[1]=dy;
            _isNan = false;
        }
        Affine2D()
        {
            _isNan = true;
        }
    public:
        Affine2D(Vec2D col0,Vec2D col1, Vec2D d)
           : Affine2D(col0[0],col1[0],col0[1],col1[1], d[0], d[1])
        {
        }
    public:
        Affine2D(R _m[], Vec2D _d)
           /* R[] not recommended because of range check */
           : Affine2D(_m[0],_m[1],_m[2],_m[3], _d[0], _d[1])
        {
        }
        Affine2D(R _m[], R _d[])
           /* R[] not recommended because of range check */
           : Affine2D(_m[0],_m[1],_m[2],_m[3], _d[0], _d[1])
        {
        }
    public:
        Affine2D(const Affine2D& a) //must be const&
        {
            m[0] = a.m[0];
            m[1] = a.m[1];
            m[2] = a.m[2];
            m[3] = a.m[3];
            d[0] = a.d[0];
            d[1] = a.d[1];
            _isNan = false;
            assert( ! a.isNan());
        }
        /* Move constructor not used.
        Affine2D(Affine2D&& a);
        Affine2D& operator(Affine2D&& a);
        */

        bool isNan() const
        {
            return _isNan;
        }
 
        /* applies the transformation */       
        Vec2D apply(Vec2D x) const
        {
            assert( ! isNan());
            return Vec2D(
                m[i00]*x[0]+m[i01]*x[1] + d[0],
                m[i10]*x[0]+m[i11]*x[1] + d[1]
                );
        };
        Vec2D operator*(Vec2D x) const
        {
            return apply(x);  //inline
        };

        Affine2D get_inverse() const
        {
            assert( ! isNan());
            R denom = m[i00] * m[i11] - m[i01]*m[i10];
            if (affine::almost_zero(denom))
                return Affine2D::NaN;

            R denom_inv = ((double)1.0)/denom;
            R m00_p =  m[i11] * denom_inv;
            R m01_p = -m[i01] * denom_inv;
            R m10_p = -m[i10] * denom_inv;
            R m11_p =  m[i00] * denom_inv;
            R dx_p = -(d[0] * m00_p + d[1] * m01_p);
            R dy_p = -(d[0] * m10_p + d[1] * m11_p);
            return Affine2D(
                m00_p,m01_p, m10_p, m11_p,
                dx_p, dy_p
                );
        };

        /* todo: Vec2D apply_inverse(Vec2D x) const {}*/

        /* operator"-" means if the (a-b).almost_zero() then a==b (appr.) */
        Affine2D operator-(const Affine2D& b) const
        {
            return Affine2D(
                m[i00]-b.m[i00],
                m[i01]-b.m[i01],
                m[i10]-b.m[i10],
                m[i11]-b.m[i11],
                d[0]  -b.d[0],
                d[1]  -b.d[1]
                );
        }
        Affine2D operator+(const Affine2D& b) const
        {
            return Affine2D(
                m[i00]+b.m[i00],
                m[i01]+b.m[i01],
                m[i10]+b.m[i10],
                m[i11]+b.m[i11],
                d[0]  +b.d[0],
                d[1]  +b.d[1]
                );
        }
        Affine2D operator-() const
        {
            return Affine2D(
                -m[i00],
                -m[i01],
                -m[i10],
                -m[i11],
                -d[0],
                -d[1]
                );
        }

        /* May not be the same value used in Vect2D */
        static constexpr R DEFAULT_TOLERANCE = 0.0000001;
        
        bool almost_zero(R tol=DEFAULT_TOLERANCE) const
        {
            assert(tol>0);
            return
                  abs(m[i00])<tol
               && abs(m[i01])<tol
               && abs(m[i10])<tol
               && abs(m[i11])<tol
               && abs(d[0])<tol
               && abs(d[1])<tol;
        }

        bool almost_equal(const Affine2D& b, R tol=DEFAULT_TOLERANCE) const {
            /*
            todo: Some parameters are more sensitive and need smaller tolerance value, such as m[i00].
            todo: check if the correct abs function is used for the given type R.
            */
            assert(tol>0);
            return
                  abs( m[i00] - b.m[i00] )<tol
               && abs( m[i01] - b.m[i01] )<tol
               && abs( m[i10] - b.m[i10] )<tol
               && abs( m[i11] - b.m[i11] )<tol
               && abs( d[0] - b.d[0] )<tol
               && abs( d[1] - b.d[1] )<tol ;            
        }
        

        Affine2D operator*(const Affine2D& b) const
        {
            return Affine2D(
                m[i00]*b.m[i00] + m[i01]*b.m[i10],
                m[i00]*b.m[i01] + m[i01]*b.m[i11],
                m[i10]*b.m[i00] + m[i11]*b.m[i10],
                m[i10]*b.m[i01] + m[i11]*b.m[i11],
                d[0]  + m[i00]*b.d[0] + m[i01]*b.d[1] ,
                d[1]  + m[i10]*b.d[0] + m[i11]*b.d[1]
                //d[0]  + m[i00]*b.d[0] + m[i01]*b.d[1] ,
                //d[1]  + m[i10]*b.d[0] + m[i11]*b.d[1]
                );
        }

        Affine2D operator*(const double r) const
        {
            return Affine2D(
                r*m[i00],
                r*m[i01],
                r*m[i10],
                r*m[i11],
                r*d[0],
                r*d[1]
                );
        }

        Affine2D& operator=(const Affine2D& rhs)
        {
            m[i00]=rhs.m[i00];
            m[i01]=rhs.m[i01];
            m[i10]=rhs.m[i10];
            m[i11]=rhs.m[i11];
            d[0]=rhs.d[1];
            d[0]=rhs.d[1];
            return *this;
        }

        std::string to_string() const
        {
            return "["
                + std::to_string(m[i00])+","+std::to_string(m[i01]) + ";"
                + std::to_string(m[i10])+","+std::to_string(m[i11]) + "]+("
                + std::to_string(d[0])+","+std::to_string(d[1]) + ")"
                ;
        }
        R det() const
        {
            return m[i00]*m[i11] - m[i01]*m[i10];
        }

        static const Affine2D& eye;
        static const Affine2D& o;
        static const Affine2D& e00;
        static const Affine2D& e01;
        static const Affine2D& e10;
        static const Affine2D& e11;
        static const Affine2D& tx;
        static const Affine2D& ty;
        static const Affine2D& NaN; //undef

    private:
        /*Used only for defining the Affine2D::Nan without the need for defining a new constructor.*/
        static Affine2D naNProvider()
        {
            Affine2D a = Affine2D::o;
            a._isNan = true;
        }
};

const Affine2D& Affine2D::eye = Affine2D(1,0,0,1, 0,0); //& or not & ?
const Affine2D& Affine2D::o = Affine2D(0,0,0,0, 0,0);
const Affine2D& Affine2D::e00 = Affine2D(1,0,0,0, 0,0);
const Affine2D& Affine2D::e01 = Affine2D(0,1,0,0, 0,0);
const Affine2D& Affine2D::e10 = Affine2D(0,0,1,0, 0,0);
const Affine2D& Affine2D::e11 = Affine2D(0,0,0,1, 0,0);
const Affine2D& Affine2D::tx = Affine2D(0,0,0,0, 1,0);
const Affine2D& Affine2D::ty = Affine2D(0,0,0,0, 0,1);
const Affine2D& Affine2D::NaN = Affine2D::naNProvider();

std::ostream& operator<<(std::ostream& os, const Affine2D& a)
{
    os << a.to_string();
    return os;
}


bool test_Affine2D(){
    bool any_failed = false;

    Affine2D o=Affine2D::o;
    std::cout << o.to_string() << std::endl;

    Affine2D a = Affine2D::eye; //Affine2D(1.0, 0,0,0);


    assert( (o-o).almost_zero()  );
    //assert( a==a );
    //assert( ! (a==o) );

    //assert( a+o==o+a );
    //assert( a - o==a );
    //assert( o + o==o );

    std::cout << a.det() << std::endl;
  
    assert( o.det()==0 );
    assert( a.det()==1.0 );

    auto b=a;
    //todo: immutability test

    std::cout << a.to_string() << std::endl;

    //Affine2D c = a + Affine2D::e00 * 2.0 - 3.0 * Affine2D::e01;
    Affine2D c = a + Affine2D::e00 + Affine2D::e01;
    //b[0,0]=2;
    assert( (c - c).almost_zero() );
    assert( Affine2D::o.almost_zero() );
    assert( (Affine2D::eye-Affine2D::eye).almost_zero() );

    Affine2D x = Affine2D::eye;
    assert( ( x + (-x) ).almost_zero() );

    assert( (Affine2D::o).almost_zero() );

    assert( (Affine2D::eye.get_inverse() - Affine2D::eye ).almost_zero() );

    assert( (Affine2D::eye * Affine2D::eye - Affine2D::eye ).almost_zero() );


    //setting up random number generation
    typedef std::mt19937 MyRNG;
    uint32_t seed_val;
    MyRNG rng;
    rng.seed(seed_val); //initialize

    std::normal_distribution<double> normal_dist(0.0, 1.0);
    for(int j=0;j<100;j++)
    {
        auto q =
            Affine2D::e00 * normal_dist(rng)+
            Affine2D::e01 * normal_dist(rng)+
            Affine2D::e10 * normal_dist(rng)+
            Affine2D::e11 * normal_dist(rng)+
            Affine2D::tx * normal_dist(rng)+
            Affine2D::ty * normal_dist(rng);

        if(abs(q.det()) > 0.0001)
        {
            if(!  (q * (q.get_inverse()) ).almost_equal( Affine2D::eye )){
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
    return ! any_failed;
};

}

int main(){

    std::cout << "Testing Vec2D" << std::endl;
    affine::test_Vec2D();
    std::cout << "Testing Affine2D" << std::endl;
    if (affine::test_Affine2D())
        std::cout << "Tests passed successfully." << std::endl;
    else
        std::cout << "Some tests failed." << std::endl;

}
