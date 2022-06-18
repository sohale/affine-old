#include "./affine.hpp"

namespace affine{

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
}
