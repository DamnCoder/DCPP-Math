//
//  matrix.h
//
//  Based upon Mathlib.h -- Copyright (c) 2005-2006 David Henry
//
//  This code is licenced under the MIT license.
//
//  Modified by Jorge López González on 17/07/12.
//

#ifndef bitthemall_matrix_h
#define bitthemall_matrix_h

#include <iostream>

#include "vector3.h"
#include "vector4.h"
#include "quaternion.h"

namespace dc
{
namespace math
{
    /////////////////////////////////////////////////////////////////////////////
	//
	// class Matrix4x4<Real> - Implement a 4x4 Matrix class that can represent
	// any 3D affine transformation.
	//
	/////////////////////////////////////////////////////////////////////////////
    
    /////////////////////////////////////////////////////////////////////////////
    //
    // class Matrix4x4<Real> implementation.
    //
    // --------------------------------------------------------------------------
    //
    // MATRIX ORGANIZATION
    //
    // The purpose of this class is so that a user might perform transformations
    // without fiddling with plus or minus signs or transposing the matrix
    // until the output "looks right".  But of course, the specifics of the
    // internal representation is important.  Not only for the implementation
    // in this file to be correct, but occasionally direct access to the
    // matrix variables is necessary, or beneficial for optimization.  Thus,
    // we document our matrix conventions here.
    //
    // Strict adherance to linear algebra rules dictates that the
    // multiplication of a 4x4 matrix by a 3D vector is actually undefined.
    // To circumvent this, we can consider the input and output vectors as
    // having an assumed fourth coordinate of 1.  Also, since the rightmost
    // column is [ 0 0 0 1 ], we can simplificate calculations ignoring
    // this last column. This is shown below:
    //
    //         | m11 m12 m13 0 | | x |   | x'|
    //         | m21 m22 m23 0 | | y | = | y'|
    //         | m31 m32 m33 0 | | z |   | z'|
    //         | tx  ty  tz  1 | | 1 |   | 1 |
    //
    // We use row vectors to represent the right, up and forward vectors
    // in the 4x4 matrix.  OpenGL uses column vectors, but the elements of
    // an OpenGL matrix are ordered in columns so that m[i][j] is in row j
    // and column i.  This is the reverse of the standard C convention in
    // which m[i][j] is in row i and column j.  The matrix should be
    // transposed before being sent to OpenGL.
    //
    //   | m11 m21 m31 tx |    |  m0  m4  m8 m12 |    |  m0  m1  m2  m3 |
    //   | m12 m22 m32 ty |    |  m1  m5  m9 m13 |    |  m4  m5  m6  m7 |
    //   | m13 m23 m33 tz |    |  m2  m6 m10 m14 |    |  m8  m9 m10 m11 |
    //   |  0   0   0  tw |    |  m3  m7 m11 m15 |    | m12 m13 m14 m15 |
    //
    //      OpenGL style          OpenGL matrix            standard C
    //                             arrangement             convention
    //
    // Fortunately, accessing to the raw matrix data via the m[] array gives
    // us the transpose matrix; i.e. in OpenGL form, so that we can directly use
    // it with glLoadMatrix{fd}() or glMultMatrix{fd}().
    //
    // Also, since the rightmost column (in standard C form) should always
    // be [ 0 0 0 1 ], and sice these values (h14, h24, h34 and tw) are
    // initialized in constructors, we don't need to modify them in our
    // matrix operations, so we don't perform useless calculations...
    //
    // The right-hand rule is used to define "positive" rotation.
    //
    //               +y                          +x'
    //                |                           |
    //                |                           |
    //                | +x                    +y' |
    //               /                            /
    //              /                            /
    //             +z                          +z'
    //
    //          initial position      Positive rotation of
    //                                 pi/2 around z-axis
    //
    /////////////////////////////////////////////////////////////////////////////
    
	template <typename Real>
	class Matrix4x4
    {
	public:

		// Constructor - Initialize the last (never used) row of the matrix
		// so that we can do any operation on matrices on the 3x4 portion
		// and forget that line which will (and should) never change.
		Matrix4x4 () :
            right(),        h1(0),
            up(),           h2(0),
            forward(),      h3(0),
            translation(),  w(1)
        {}
        //h14 (0.0f), h24 (0.0f), h34 (0.0f), tw (1.0f) { }
        
        Matrix4x4 (const Vector4<Real>& row1,
                   const Vector4<Real>& row2,
                   const Vector4<Real>& row3,
                   const Vector4<Real>& row4):
            row1(row1),
            row2(row2),
            row3(row3),
            row4(row4)
        {}
        
        Matrix4x4 (const Matrix4x4<Real>& copy):
            row1(copy.row1),
            row2(copy.row2),
            row3(copy.row3),
            row4(copy.row4)
        {}
        
        
		// Note that we don't define the copy constructor and let the compiler
		// doing it itself because such initialization is not necessary
		// since the source matrix has its last row already initialized...
        
	public:
        // Matrix operations
        
        // Matrix comparison
        const bool operator ==(const Matrix4x4<Real>& m) const;
        const bool operator !=(const Matrix4x4<Real>& m) const;
        
        // Combined assignment operators to conform to C notation convention
        Matrix4x4<Real>& operator+= (const Matrix4x4<Real>& m);
        Matrix4x4<Real>& operator-= (const Matrix4x4<Real>& m);
        
        Matrix4x4<Real> operator+ (const Matrix4x4<Real>& m) const;
        Matrix4x4<Real> operator- (const Matrix4x4<Real>& m) const;
        
		void Identity();
		void Transpose();
        
		void Translate(const Vector3<Real>& p);
        void Rotate(const Vector3<Real>& rotation);
        void Rotate(const Vector3<Real>& axis, const Real theta);
        void Scale(const Vector3<Real>& s);
        
        void FromQuaternion (const Quaternion<Real> &q);
        Quaternion<Real> ToQuaternion() const;
        
        void FromEulerAngles (const Vector3<Real>& rotation);
        void FromEulerAngles (Real x, Real y, Real z);
        Vector3<Real> EulerAngles () const;
        
        const bool IsSingular() const;
        const Real Determinant3x3() const;
        
		Vector3<Real> TransformPosition(const Vector3<Real>& position) const;
        Vector3<Real> InverseTransformPosition (const Vector3<Real>& position) const;
        
        Vector3<Real> TransformRotation (const Vector3<Real>& rotation) const;
		Vector3<Real> InverseTransformRotation (const Vector3<Real>& rotation) const;
        
    // Camera transformations
    public:
        
        Matrix4x4<Real> FrustumMatrix       (Real left, Real right, Real bottom, Real top, Real near, Real far);
        Matrix4x4<Real> OrthoMatrix         (Real left, Real right, Real bottom, Real top, Real near, Real far);
        Matrix4x4<Real> PerspectiveMatrix   (Real fovY, Real aspect, Real near, Real far);
        
		// Accessor.  This allows to use the matrix object
		// like an array of Real. For example:
		// Matrix4x4<float> mat;
		// float f = mat[4]; // access to m21
		operator const Real *() { return m; }
        
    public:
        // static operations for class methods
        static Matrix4x4<Real> Inverse(const Matrix4x4<Real>& mat);
        static Matrix4x4<Real> IdentityMatrix();
        
	public:
		// Member variables
        
		// The values of the matrix.  Basically the upper 3x3 portion
		// contains a linear transformation, and the last column is the
		// translation portion. Here data is transposed, see the Mathlib.inl
		// for more details.
		union
        {
			struct
            {
				Real m11, m12, m13, h14;
				Real m21, m22, m23, h24;
				Real m31, m32, m33, h34;
				Real tx,  ty,  tz,  tw;
			};
            
            // Basis vectors from matrix
            struct
            {
                Vector3<Real> right;        float h1;
                Vector3<Real> up;           float h2;
                Vector3<Real> forward;      float h3;
                Vector3<Real> translation;  float w;
            };
            
            struct
            {
                Vector4<Real> row1;
                Vector4<Real> row2;
                Vector4<Real> row3;
                Vector4<Real> row4;
            };
            
			// Access to raw packed matrix data (usefull for
			// glLoadMatrixf () and glMultMatrixf ())
			Real m[16];
		};
	};
    
    // Predefined Matrix4x4 types
	typedef Matrix4x4<float>    Matrix4x4f;
	typedef Matrix4x4<double>   Matrix4x4d;
    
    //
	// Nonmember Matrix4x4 functions
	//
    
    // Matrix concatenation
	template <typename Real>
	Matrix4x4<Real> operator* (const Matrix4x4<Real>& m1, const Matrix4x4<Real>& m2);
    
    template <typename Real>
	Matrix4x4<Real> operator*= (const Matrix4x4<Real>& m1, const Matrix4x4<Real>& m2);
    
    template <typename Real>
    std::ostream& operator << (std::ostream& output, const Matrix4x4<Real>& m);
    
    /////////////////////////////////////////////////////////////////////////////
    //
    // class Vector3<Real> member functions implementation.
    //
    /////////////////////////////////////////////////////////////////////////////
    
    template <typename Real>
    inline
    Matrix4x4<Real>&
    Matrix4x4<Real>::operator +=(const Matrix4x4<Real>& m)
    {
        m11 += m.m11;   m12 += m.m12;   m13 += m.m13;
        m21 += m.m21;   m22 += m.m22;   m23 += m.m23;
        m31 += m.m31;   m32 += m.m32;   m33 += m.m33;
        
        tx += m.tx; ty += m.ty; tz += m.tz;
        
        return *this;
    }
    
    template <typename Real>
    inline
    Matrix4x4<Real>&
    Matrix4x4<Real>::operator -=(const Matrix4x4<Real>& m)
    {
        m11 -= m.m11; m12 -= m.m12; m13 -= m.m13;
        m21 -= m.m21; m22 -= m.m22; m23 -= m.m23;
        m31 -= m.m31; m32 -= m.m32; m33 -= m.m33;
        
        tx -= m.tx; ty -= m.ty; tz -= m.tz;
        
        return *this;
    }
    
    template <typename Real>
    inline
    Matrix4x4<Real>
    Matrix4x4<Real>::operator+ (const Matrix4x4<Real>& m) const
    {
        Matrix4x4<Real> res;
        res.m11 = m11 + m.m11; res.m12 = m12 + m.m12; res.m13 = m13 + m.m13;
        res.m21 = m21 + m.m21; res.m22 = m22 + m.m22; res.m23 = m23 + m.m23;
        res.m31 = m31 + m.m31; res.m32 = m32 + m.m32; res.m33 = m33 + m.m33;
        
        res.tx = tx + m.tx; res.ty = ty + m.ty; res.tz = tz + m.tz;

        return res;
    }
    
    template <typename Real>
    inline
    Matrix4x4<Real>
    Matrix4x4<Real>::operator- (const Matrix4x4<Real>& m) const
    {
        Matrix4x4<Real> res;
        res.m11 = m11 - m.m11; res.m12 = m12 - m.m12; res.m13 = m13 - m.m13;
        res.m21 = m21 - m.m21; res.m22 = m22 - m.m22; res.m23 = m23 - m.m23;
        res.m31 = m31 - m.m31; res.m32 = m32 - m.m32; res.m33 = m33 - m.m33;
        
        res.tx = tx - m.tx; res.ty = ty - m.ty; res.tz = tz - m.tz;
        
        return res;
    }
     
    // --------------------------------------------------------------------------
    // Matrix4x4::operator ==
    //
    // Matrix equals comparison.
    // --------------------------------------------------------------------------
    
    template <typename Real>
    inline
    const bool
    Matrix4x4<Real>::operator==(const Matrix4x4<Real> &m) const
    {
        return row1 == m.row1 && row2 == m.row2 && row3 == m.row3 && row4 == m.row4;
    }
        
    // --------------------------------------------------------------------------
    // Matrix4x4::operator !=
    //
    // Matrix distinct comparison.
    // --------------------------------------------------------------------------
    
    template <typename Real>
    inline
    const bool
    Matrix4x4<Real>::operator!=(const Matrix4x4<Real> &m) const
    {
        return row1 != m.row1 || row2 != m.row2 || row3 != m.row3 || row4 != m.row4;
    }
        
        
    // --------------------------------------------------------------------------
    // Matrix4x4::identity
    //
    // Set matrix to identity.
    // --------------------------------------------------------------------------
    
    template <typename Real>
    inline void
    Matrix4x4<Real>::Identity()
    {
        m11 = 1.0; m12 = 0.0; m13 = 0.0; h14 = 0.0;
        m21 = 0.0; m22 = 1.0; m23 = 0.0; h24 = 0.0;
        m31 = 0.0; m32 = 0.0; m33 = 1.0; h34 = 0.0;
        
        tx = 0.0; ty = 0.0; tz = 0.0; tw = 1.0;
    }
    
    
    // --------------------------------------------------------------------------
    // Matrix4x4::transpose
    //
    // Transpose the current matrix.
    // --------------------------------------------------------------------------
    
    template <typename Real>
    inline
    void Matrix4x4<Real>::Transpose()
    {
        Matrix4x4<Real> res;
        
        res.m11 = m11; res.m12 = m21; res.m13 = m31; res.h14 = tx;
        res.m21 = m12; res.m22 = m22; res.m23 = m32; res.h24 = ty;
        res.m31 = m13; res.m32 = m23; res.m33 = m33; res.h34 = tz;
        
        res.tx = h14; res.ty = h24; res.tz = h34; res.tw = tw;
        *this = res;
    }
    
    // --------------------------------------------------------------------------
	// Translate
	//
	// Build a translation matrix given a translation vector.
	// --------------------------------------------------------------------------
    
	template <typename Real>
	inline
    void
    Matrix4x4<Real>::Translate (const Vector3<Real>& v)
	{
        translation += v;
	}
    
    // --------------------------------------------------------------------------
	// Rotate
	//
	// Build a rotation matrix given the axis and angle.
	// --------------------------------------------------------------------------
    
    template <typename Real>
	inline
    void
    Matrix4x4<Real>::Rotate(const Vector3<Real>& rotation)
    {
        FromEulerAngles(rotation.x, rotation.y, rotation.z);
    }
    
    template <typename Real>
	inline
    void
    Matrix4x4<Real>::Rotate(const Vector3<Real>& axis, const Real theta)
    {
        
        // Quick sanity check to make sure they passed in a unit vector
        // to specify the axis
        assert (axis.Length() <= 1.0);
        
        // Get sin and cosine of rotation angle
        Real s = std::sin (theta);
        Real c = std::cos (theta);
        
        // Compute 1 - cos(theta) and some common subexpressions
        Real a = 1.0 - c;
        Real ax = a * axis.x;
        Real ay = a * axis.y;
        Real az = a * axis.z;
        
        // Set the matrix elements.  There is still a little more
        // opportunity for optimization due to the many common
        // subexpressions.  We'll let the compiler handle that...
        m11 = (ax * axis.x) + c;
        m12 = (ax * axis.y) + (axis.z * s);
        m13 = (ax * axis.z) - (axis.y * s);
        
        m21 = (ay * axis.x) - (axis.z * s);
        m22 = (ay * axis.y) + c;
        m23 = (ay * axis.z) + (axis.x * s);
        
        m31 = (az * axis.x) + (axis.y * s);
        m32 = (az * axis.y) - (axis.x * s);
        m33 = (az * axis.z) + c;
        
        return ;
    }
    
	// --------------------------------------------------------------------------
	// Scale
	//
	// Setup the matrix to perform scale on each axis.  For uniform scale by k,
	// use a vector of the form Vector3( k, k, k ).
	//
	// --------------------------------------------------------------------------
    
	template <typename Real>
	inline
    void
    Matrix4x4<Real>::Scale(const Vector3<Real>& s)
	{
        // Set the matrix elements.  Pretty straightforward
        m11 = s.x;
        m22 = s.y;
        m33 = s.z;
	}
    
    template <typename Real>
    inline
    Quaternion<Real> Matrix4x4<Real>::ToQuaternion() const
    {
        return Quaternion<Real>(*this);
        /*
        Real q0 = ( m11 + m22 + m33 + 1.0) / 4.0;
        Real q1 = ( m11 - m22 - m33 + 1.0) / 4.0;
        Real q2 = (-m11 + m22 - m33 + 1.0) / 4.0;
        Real q3 = (-m11 - m22 + m33 + 1.0) / 4.0;
        
        if(q0 < 0.0) q0 = 0.0;
        if(q1 < 0.0) q1 = 0.0;
        if(q2 < 0.0) q2 = 0.0;
        if(q3 < 0.0) q3 = 0.0;
        
        q0 = sqrt(q0);
        q1 = sqrt(q1);
        q2 = sqrt(q2);
        q3 = sqrt(q3);
        
        if(q0 >= q1 && q0 >= q2 && q0 >= q3)
        {
            q0 *= +1.0f;
            q1 *= Sign(m32 - m23);
            q2 *= Sign(m13 - m31);
            q3 *= Sign(m21 - m12);
        }
        else if(q1 >= q0 && q1 >= q2 && q1 >= q3)
        {
            q0 *= Sign(m32 - m23);
            q1 *= +1.0f;
            q2 *= Sign(m21 + m12);
            q3 *= Sign(m13 + m31);
        }
        else if(q2 >= q0 && q2 >= q1 && q2 >= q3)
        {
            q0 *= Sign(m13 - m31);
            q1 *= Sign(m21 + m12);
            q2 *= +1.0f;
            q3 *= Sign(m32 + m23);
        }
        else if(q3 >= q0 && q3 >= q1 && q3 >= q2)
        {
            q0 *= Sign(m21 - m12);
            q1 *= Sign(m31 + m13);
            q2 *= Sign(m32 + m23);
            q3 *= +1.0f;
        }
        
        Real length = Normal(q0, q1, q2, q3);
        
        return Quaternion<Real>(q0/length, q1/length, q2/length, q3/length);
         */
        
        // Algorithm in Ken Shoemake's article in 1987 SIGGRAPH course notes
        // article "Quaternion Calculus and Fast Animation".
        
    }
    
    // --------------------------------------------------------------------------
    // Matrix4x4::fromQuaternion
    //
    // Convert a quaternion to a matrix.
    // --------------------------------------------------------------------------
    
    template <typename Real>
    inline void
    Matrix4x4<Real>::FromQuaternion (const Quaternion<Real> &q)
    {
        // Compute a few values to optimize common subexpressions
        Real ww = 2.0 * q.w;
        Real xx = 2.0 * q.x;
        Real yy = 2.0 * q.y;
        Real zz = 2.0 * q.z;
        
        // Set the matrix elements.  There is still a little more
        // opportunity for optimization due to the many common
        // subexpressions.  We'll let the compiler handle that...
        m11 = 1.0 - (yy * q.y) - (zz * q.z);
        m12 = (xx * q.y) + (ww * q.z);
        m13 = (xx * q.z) - (ww * q.y);
        
        m21 = (xx * q.y) - (ww * q.z);
        m22 = 1.0 - (xx * q.x) - (zz * q.z);
        m23 = (yy * q.z) + (ww * q.x);
        
        m31 = (xx * q.z) + (ww * q.y);
        m32 = (yy * q.z) - (ww * q.x);
        m33 = 1.0 - (xx * q.x) - (yy * q.y);
        
        // Reset the translation portion
        tx = ty = tz = 0.0;
    }
    
    // --------------------------------------------------------------------------
	// Matrix4x4::fromEulerAngles
	//
	// Setup a rotation matrix, given three X-Y-Z rotation angles. The
	// rotations are performed first on x-axis, then y-axis and finaly z-axis.
	// --------------------------------------------------------------------------
    
    template <typename Real>
	inline void
	Matrix4x4<Real>::FromEulerAngles (const Vector3<Real>& rotation)
	{
        FromEulerAngles(rotation.x, rotation.y, rotation.z);
    }
    
	template <typename Real>
	inline void
	Matrix4x4<Real>::FromEulerAngles (Real x, Real y, Real z)
	{
        /*
        // Fetch sine and cosine of angles
        const Real cr = std::cos (x);
        const Real sr = std::sin (x);
        
        const Real cp = std::cos (y);
        const Real sp = std::sin (y);
        
        const Real cy = std::cos (z);
        const Real sy = std::sin (z);
        
        const Real srsp = sr * sp;
        const Real crsp = cr * sp;
        
        // Fill in the matrix elements
        m11 =  (cp * cy);   m12 =  (srsp * cy) - (cr * sy); m13 =  (crsp * cy) + (sr * sy);
        m21 =  (cp * sy);   m22 =  (srsp * sy) + (cr * cy); m23 =  (crsp * sy) - (sr * cy);
        m31 = -(sp);        m32 =  (sr * cp);               m33 =  (cr * cp);
        //*/
        
        /*
        // Fetch sine and cosine of angles
        const Real cr = std::cos (x);
        const Real sr = std::sin (x);
        
        const Real cp = std::cos (y);
        const Real sp = std::sin (y);
        
        const Real cy = std::cos (z);
        const Real sy = std::sin (z);
        
        const Real srsp = sr * sp;
        const Real crsp = cr * sp;
        
        // Fill in the matrix elements
        m11 =  (cp * cy);   m12 =  (srsp * cy) + (cr * sy); m13 = -(crsp * cy) + (sr * sy);
        m21 = -(cp * sy);   m22 = -(srsp * sy) + (cr * cy); m23 =  (crsp * sy) + (sr * cy);
        m31 =  (sp);        m32 = -(sr * cp);               m33 =  (cr * cp);
        //*/
        //*
        // Fetch sine and cosine of angles
        // Bank
        Real cb = std::cos (x);
        Real sb = std::sin (x);
        
        // Heading
        Real ch = std::cos (y);
        Real sh = std::sin (y);
        
        // Attitude
        Real ca = std::cos (z);
        Real sa = std::sin (z);
        
        Real sacb = sa * cb;
        Real sasb = sa * sb;
        
        // Fill in the matrix elements
        m11 =  (ch * ca);   m12 = -(ch * sacb) + (sh * sb); m13 =  (ch * sasb) + (sh * cb);
        m21 =  (sa);        m22 =  (ca * cb);               m23 = -(ca * sb);
        m31 = -(sh * ca);   m32 = -(sh * sacb) + (ch * sb); m33 = -(sh * sasb) + (ch * cb);
        //*/
	}

    // --------------------------------------------------------------------------
	// Matrix4x4::toEulerAngles
	//
	// Setup the euler angles in radians, given a rotation matrix. The rotation
	// matrix could have been obtained from euler angles given the expression:
	//   M = Y * X * Z
	// where X, Y and Z are rotation matrices about X, Y and Z axes.
	// This is the "opposite" of the fromEulerAngles function.
	// --------------------------------------------------------------------------
    
	template <typename Real>
	inline
    Vector3<Real>
	Matrix4x4<Real>::EulerAngles() const
	{
        /*
         if (m.m10 > 0.998) { // singularity at north pole
         heading = Math.atan2(m.m02,m.m22);
         attitude = Math.PI/2;
         bank = 0;
         return;
         }
         if (m.m10 < -0.998) { // singularity at south pole
         heading = Math.atan2(m.m02,m.m22);
         attitude = -Math.PI/2;
         bank = 0;
         return;
         }
         
         bank = Math.atan2(-m.m12,m.m11);       // rotX
         heading = Math.atan2(-m.m20,m.m00);    // rotY
         attitude = Math.asin(m.m10);           // rotZ
         
         */
        //*
        float bank, heading, attitude;
        
        if(0.998 < m21)
        {
            bank  = 0;
            heading = std::atan2( m13, m33);
            attitude =  MO_PI_OVER_2;
            return Vector3<Real>(bank, heading, attitude);
        }
        
        if(m21 < -0.998)
        {
            bank  = 0;
            heading = std::atan2( m13, m33);
            attitude =  -MO_PI_OVER_2;
            return Vector3<Real>(bank, heading, attitude);
        }
        
        bank  = std::atan2( -m23, m22);
        heading = std::atan2( -m31, m11);
        attitude = std::asin ( m21);
        
        return Vector3<Real>(bank, heading, attitude);
        //*/
        /*
        // Compute Y-axis angle
        const Real y = std::asin (m31);
        
        // Compute cos and one over cos for optimization
        const Real cy = std::cos (y);
        const Real oneOverCosY = 1.0 / cy;
        
        // No gimball lock
        const Real x = std::atan2 (-m32 * oneOverCosY, m33 * oneOverCosY);
        const Real z = std::atan2 (-m21 * oneOverCosY, m11 * oneOverCosY);
        
        return Vector3<Real>(x, y, z);
        //*/
	}
    
    // --------------------------------------------------------------------------
	// Matrix4x4::TransformPosition
	//
	// Transform a point by the matrix.
	// --------------------------------------------------------------------------
    
	template <typename Real>
	inline
    Vector3<Real>
	Matrix4x4<Real>::TransformPosition(const Vector3<Real>& position) const
	{
        // Grind through the linear algebra.
        return Vector3<Real>((position.x * m11) + (position.y * m21) + (position.z * m31) + tx,
                             (position.x * m12) + (position.y * m22) + (position.z * m32) + ty,
                             (position.x * m13) + (position.y * m23) + (position.z * m33) + tz );
	}
    
    // --------------------------------------------------------------------------
	// Matrix4x4::InverseTransformPosition
	//
	// Translate a point by the inverse matrix.
	// --------------------------------------------------------------------------
    
	template <typename Real>
	inline
    Vector3<Real>
	Matrix4x4<Real>::InverseTransformPosition (const Vector3<Real> &v) const
	{
        return Vector3<Real>(v.x - tx, v.y - ty, v.z - tz);
	}
    
    // --------------------------------------------------------------------------
	// Matrix4x4::TransformRotation
	//
	// Rotate a point by the 3x3 upper left portion of the matrix.
	// --------------------------------------------------------------------------
    
	template <typename Real>
	inline
    Vector3<Real>
	Matrix4x4<Real>::TransformRotation(const Vector3<Real>& rotation) const
	{
        return Vector3<Real>((rotation.x * m11) + (rotation.y * m21) + (rotation.z * m31),
                             (rotation.x * m12) + (rotation.y * m22) + (rotation.z * m32),
                             (rotation.x * m13) + (rotation.y * m23) + (rotation.z * m33));
	}
    
    // --------------------------------------------------------------------------
	// Matrix4x4::inverseRotate
	//
	// Rotate a point by the inverse 3x3 upper left portion of the matrix.
	// --------------------------------------------------------------------------
    
	template <typename Real>
	inline
    Vector3<Real>
	Matrix4x4<Real>::InverseTransformRotation (const Vector3<Real> &v) const
	{
        return Vector3<Real> ((v.x * m11) + (v.y * m12) + (v.z * m13),
                              (v.x * m21) + (v.y * m22) + (v.z * m23),
                              (v.x * m31) + (v.y * m32) + (v.z * m33));
	}
    
    // --------------------------------------------------------------------------
	// FrustumMatrix
	//
	// Setup a frustum matrix given the left, right, bottom, top, near, and far
	// values for the frustum boundaries.
	// --------------------------------------------------------------------------
    
	template <typename Real>
	inline
    Matrix4x4<Real>
	Matrix4x4<Real>::FrustumMatrix (Real l, Real r, Real b, Real t, Real n, Real f)
	{
        assert (n >= 0.0);
        assert (f >= 0.0);
        
        Matrix4x4<Real> res;
        
        Real width  = r - l;
        Real height = t - b;
        Real depth  = f - n;
        
        res.m[0] = (2 * n) / width;
        res.m[1] = 0.0;
        res.m[2] = 0.0;
        res.m[3] = 0.0;
        
        res.m[4] = 0.0;
        res.m[5] = (2 * n) / height;
        res.m[6] = 0.0;
        res.m[7] = 0.0;
        
        res.m[8] = (r + l) / width;
        res.m[9] = (t + b) / height;
        res.m[10]= -(f + n) / depth;
        res.m[11]= -1.0;
        
        res.m[12]= 0.0;
        res.m[13]= 0.0;
        res.m[14]= -(2 * f * n) / depth;
        res.m[15]= 0.0;
        
        return res;
	}
    
    // --------------------------------------------------------------------------
	// OrthoMatrix
	//
	// Setup an orthographic Matrix4x4 given the left, right, bottom, top, near,
	// and far values for the frustum boundaries.
	// --------------------------------------------------------------------------
    
	template <typename Real>
	inline
    Matrix4x4<Real>
	Matrix4x4<Real>::OrthoMatrix (Real l, Real r, Real b, Real t, Real n, Real f)
	{
        Matrix4x4<Real> res;
        
        Real width  = r - l;
        Real height = t - b;
        Real depth  = f - n;
        
        res.m[0] =  2.0 / width;
        res.m[1] =  0.0;
        res.m[2] =  0.0;
        res.m[3] =  0.0;
        
        res.m[4] =  0.0;
        res.m[5] =  2.0 / height;
        res.m[6] =  0.0;
        res.m[7] =  0.0;
        
        res.m[8] =  0.0;
        res.m[9] =  0.0;
        res.m[10]= -2.0 / depth;
        res.m[11]=  0.0;
        
        res.m[12]= -(r + l) / width;
        res.m[13]= -(t + b) / height;
        res.m[14]= -(f + n) / depth;
        res.m[15]=  1.0;
        
        return res;
	}
    
    // --------------------------------------------------------------------------
	// PerspectiveMatrix
	//
	// Setup a perspective matrix given the field-of-view in the Y direction
	// in degrees, the aspect ratio of Y/X, and near and far plane distances.
	// --------------------------------------------------------------------------
    
	template <typename Real>
	inline
    Matrix4x4<Real>
	Matrix4x4<Real>::PerspectiveMatrix (Real fovY, Real aspect, Real n, Real f)
	{
        Matrix4x4<Real> res;
        
        Real angle;
        Real cot;
        
        angle = fovY / 2.0;
        angle = DegToRad(angle);
        
        cot = std::cos (angle) / std::sin (angle);
        
        res.m[0] = cot / aspect;
        res.m[1] = 0.0;
        res.m[2] = 0.0;
        res.m[3] = 0.0;
        
        res.m[4] = 0.0;
        res.m[5] = cot;
        res.m[6] = 0.0;
        res.m[7] = 0.0;
        
        res.m[8] = 0.0;
        res.m[9] = 0.0;
        res.m[10]= -(f + n) / (f - n);
        res.m[11]= -1.0;
        
        res.m[12]= 0.0;
        res.m[13]= 0.0;
        res.m[14]= -(2 * f * n) / (f - n);
        res.m[15]= 0.0;
        
        return res;
	}
    
    // --------------------------------------------------------------------------
	// Matrix4x4::Determinant
	//
	// Compute the determinant of the 3x3 portion of the matrix.
	// --------------------------------------------------------------------------
    
	template <typename Real>
	inline
    const Real
    Matrix4x4<Real>::Determinant3x3 () const
	{
        return
        m11 * ((m22 * m33) - (m23 * m32))
        +
        m12 * ((m23 * m31) - (m21 * m33))
        +
        m13 * ((m21 * m32) - (m22 * m31));
	}
    
    // --------------------------------------------------------------------------
	// Matrix4x4::IsSingular
	//
	// If the determinant of the 3x3 matrix part is 0, then is singular, which
    // means there are no Inverse of the 3x3 matrix part.
	// --------------------------------------------------------------------------
    
	template <typename Real>
	inline
    const bool
    Matrix4x4<Real>::IsSingular() const
	{
        return Approximately(Determinant3x3(), 0.0f);
	}
    
    /////////////////////////////////////////////////////////////////////////////
    //
    // class Vector3<Real> static functions implementation.
    //
    /////////////////////////////////////////////////////////////////////////////
    
    // --------------------------------------------------------------------------
	// Inverse
	//
	// Compute the inverse of a matrix.  We use the classical adjoint divided
	// by the determinant method.
    //
    // Caution! This function did not check if the matrix is singular
	// --------------------------------------------------------------------------
    
	template <typename Real>
	inline
    Matrix4x4<Real>
	Matrix4x4<Real>::Inverse(const Matrix4x4<Real>& mat)
	{
        // Compute the determinant of the 3x3 portion
        const Real det = mat.Determinant3x3 ();
        
        // Compute one over the determinant, so we divide once and
        // can *multiply* per element
        const Real oneOverDet = 1.0 / det;
        
        // Compute the 3x3 portion of the inverse, by
        // dividing the adjoint by the determinant
        Matrix4x4<Real> res;
        
        res.m11 = ((mat.m22 * mat.m33) - (mat.m23 * mat.m32)) * oneOverDet;
        res.m12 = ((mat.m13 * mat.m32) - (mat.m12 * mat.m33)) * oneOverDet;
        res.m13 = ((mat.m12 * mat.m23) - (mat.m13 * mat.m22)) * oneOverDet;
        
        res.m21 = ((mat.m23 * mat.m31) - (mat.m21 * mat.m33)) * oneOverDet;
        res.m22 = ((mat.m11 * mat.m33) - (mat.m13 * mat.m31)) * oneOverDet;
        res.m23 = ((mat.m13 * mat.m21) - (mat.m11 * mat.m23)) * oneOverDet;
        
        res.m31 = ((mat.m21 * mat.m32) - (mat.m22 * mat.m31)) * oneOverDet;
        res.m32 = ((mat.m12 * mat.m31) - (mat.m11 * mat.m32)) * oneOverDet;
        res.m33 = ((mat.m11 * mat.m22) - (mat.m12 * mat.m21)) * oneOverDet;
        
        // Compute the translation portion of the inverse
        res.tx = -((mat.tx * res.m11) + (mat.ty * res.m21) + (mat.tz * res.m31));
        res.ty = -((mat.tx * res.m12) + (mat.ty * res.m22) + (mat.tz * res.m32));
        res.tz = -((mat.tx * res.m13) + (mat.ty * res.m23) + (mat.tz * res.m33));
        
        // Return it.
        return res;
	}
    
    template <typename Real>
	inline
    Matrix4x4<Real>
	Matrix4x4<Real>::IdentityMatrix()
    {
        Matrix4x4<Real> res;
        
        res.m11 = 1.0; res.m12 = 0.0; res.m13 = 0.0; res.h14 = 0.0;
        res.m21 = 0.0; res.m22 = 1.0; res.m23 = 0.0; res.h24 = 0.0;
        res.m31 = 0.0; res.m32 = 0.0; res.m33 = 1.0; res.h34 = 0.0;
        
        res.tx = 0.0; res.ty = 0.0; res.tz = 0.0; res.tw = 1.0;
        
        return res;
    }
    
    /////////////////////////////////////////////////////////////////////////////
    //
    // class Vector3<Real> non - member functions implementation.
    //
    /////////////////////////////////////////////////////////////////////////////
    
    // --------------------------------------------------------------------------
	// Matrix4x4 * Matrix4x4
	//
	// Matrix concatenation.
	//
	// We also provide a *= operator, as per C convention.
	// --------------------------------------------------------------------------
    
	template <typename Real>
	inline Matrix4x4<Real>
	operator* (const Matrix4x4<Real>& m1, const Matrix4x4<Real>& m2)
	{
        Matrix4x4<Real> res;
        
        // Compute the left 4x3 (linear transformation) portion
        res.m11 = (m1.m11 * m2.m11) + (m1.m12 * m2.m21) + (m1.m13 * m2.m31);
        res.m12 = (m1.m11 * m2.m12) + (m1.m12 * m2.m22) + (m1.m13 * m2.m32);
        res.m13 = (m1.m11 * m2.m13) + (m1.m12 * m2.m23) + (m1.m13 * m2.m33);
        
        res.m21 = (m1.m21 * m2.m11) + (m1.m22 * m2.m21) + (m1.m23 * m2.m31);
        res.m22 = (m1.m21 * m2.m12) + (m1.m22 * m2.m22) + (m1.m23 * m2.m32);
        res.m23 = (m1.m21 * m2.m13) + (m1.m22 * m2.m23) + (m1.m23 * m2.m33);
        
        res.m31 = (m1.m31 * m2.m11) + (m1.m32 * m2.m21) + (m1.m33 * m2.m31);
        res.m32 = (m1.m31 * m2.m12) + (m1.m32 * m2.m22) + (m1.m33 * m2.m32);
        res.m33 = (m1.m31 * m2.m13) + (m1.m32 * m2.m23) + (m1.m33 * m2.m33);
        
        // Compute the translation portion
        res.tx = (m1.tx * m2.m11) + (m1.ty * m2.m21) + (m1.tz * m2.m31) + m2.tx;
        res.ty = (m1.tx * m2.m12) + (m1.ty * m2.m22) + (m1.tz * m2.m32) + m2.ty;
        res.tz = (m1.tx * m2.m13) + (m1.ty * m2.m23) + (m1.tz * m2.m33) + m2.tz;
        
        return res;
	}
    
    template <typename Real>
	inline Matrix4x4<Real> &
	operator*= (Matrix4x4<Real>& m1, const Matrix4x4<Real>& m2)
	{
        m1 = m1 * m2;
        return m1;
	}
    
    // --------------------------------------------------------------------------
	// ostream << Matrix4x4
	//
	// Standard output of a Matrix4x4.
	// --------------------------------------------------------------------------
    
    template <typename Real>
    inline std::ostream&
    operator<< (std::ostream& output, const Matrix4x4<Real>& m)
    {
        output << m.row1 << "\n";
        output << m.row2 << "\n";
        output << m.row3 << "\n";
        output << m.row4 << "\n";
        return output;
    }

}
}

#endif
