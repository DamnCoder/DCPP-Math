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
	//	 Column - Row								  Row - Column
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
	//                |__+x					 +y'__|
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
		// ===========================================================
		// Constant / Enums / Typedefs internal usage
		// ===========================================================
		
		// ===========================================================
		// Static fields / methods
		// ===========================================================
	public:
		// Matrix-builder functions
		static Matrix4x4<Real> Frustum		(Real left, Real right, Real bottom, Real top, Real near, Real far);
		
		static Matrix4x4<Real> Orthographic	(Real left, Real right, Real bottom, Real top, Real near, Real far);
		static Matrix4x4<Real> Perspective	(Real fovY, Real aspect, Real near, Real far);
		
		static Matrix4x4<Real> LookAt		(const Vector3<Real>& eye, const Vector3<Real>& target, const Vector3<Real>& camUp);
		
		static Matrix4x4<Real> Inverse		(const Matrix4x4<Real>& mat);
		static Matrix4x4<Real> Identity		();
		
		// ===========================================================
		// Inner and Anonymous Classes
		// ===========================================================
		
		// ===========================================================
		// Getter & Setter
		// ===========================================================
	public:
		// Basis vectors from matrix
		const Vector3<Real>	Right()			const	{ return Vector3<Real>(m11, m12, m13); }
		const Vector3<Real>	Up()			const	{ return Vector3<Real>(m21, m22, m23); }
		const Vector3<Real>	Forward()		const	{ return Vector3<Real>(m31, m32, m33); }
		
		Vector3<Real>		Translation()			{ return Vector3<Real>(tx, ty, tz); }
		const Vector3<Real>	Translation()	const	{ return Vector3<Real>(tx, ty, tz); }
		
		const Vector4<Real>	Row1()			const	{ return Vector4<Real>(m11, m12, m13, h14); }
		const Vector4<Real>	Row2()			const	{ return Vector4<Real>(m21, m22, m23, h24); }
		const Vector4<Real>	Row3()			const	{ return Vector4<Real>(m31, m32, m33, h34); }
		const Vector4<Real>	Row4()			const	{ return Vector4<Real>(tx, ty, tz, tw); }

		void				FromQuaternion (const Quaternion<Real> &q);
		Quaternion<Real>	ToQuaternion() const;
		
		void				FromEulerAngles (const Vector3<Real>& rotation);
		void				FromEulerAngles (Real x, Real y, Real z);
		Vector3<Real>		EulerAngles () const;
		
		const bool			IsSingular() const;
		const Real			Determinant3x3() const;
		
		// ===========================================================
		// Constructors
		// ===========================================================
	public:
		
		Matrix4x4(): h14(0), h24(0), h34(0), tw(1)
        {}
		
		Matrix4x4(const Matrix4x4& copy):
			m11(copy.m11), m12(copy.m12), m13(copy.m13), h14(copy.h14),
			m21(copy.m21), m22(copy.m22), m23(copy.m23), h24(copy.h24),
			m31(copy.m31), m32(copy.m32), m33(copy.m33), h34(copy.h34),
			tx(copy.tx), ty(copy.ty), tz(copy.tz), tw(copy.tw)
		{}

		// ===========================================================
		// Methods for/from SuperClass/Interfaces
		// ===========================================================
	public:
		void				operator= (const Matrix4x4<Real>& m);
		
        // Matrix comparison
        const bool			operator==(const Matrix4x4<Real>& m) const;
        const bool			operator!=(const Matrix4x4<Real>& m) const;
        
        // Combined assignment operators to conform to C notation convention
        Matrix4x4<Real>&	operator+= (const Matrix4x4<Real>& m);
        Matrix4x4<Real>&	operator-= (const Matrix4x4<Real>& m);
        
        Matrix4x4<Real>		operator+ (const Matrix4x4<Real>& m) const;
        Matrix4x4<Real>		operator- (const Matrix4x4<Real>& m) const;
		
		// ===========================================================
		// Methods
		// ===========================================================
	public:
		void			Identify();
		void			Transpose();
		
		void			Translate	(const Vector3<Real>& position);
		void			Rotate		(const Vector3<Real>& rotation);
		void			Rotate		(const Vector3<Real>& axis, const Real theta);
		void			Scale		(const Vector3<Real>& scale);
		
		Vector3<Real>	TransformPosition			(const Vector3<Real>& position) const;
		Vector3<Real>	InverseTransformPosition	(const Vector3<Real>& position) const;
		
		Vector3<Real>	TransformRotation			(const Vector3<Real>& rotation) const;
		Vector3<Real>	InverseTransformRotation	(const Vector3<Real>& rotation) const;
		
	public:
		// Accessor.  This allows to use the matrix object
		// like an array of Real. For example:
		// Matrix4x4<float> mat;
		// float f = mat[4]; // access to m21
		operator const Real *()			{ return m; }
		operator const Real *() const	{ return m; }

		// ===========================================================
		// Fields
		// ===========================================================
	public:
		// Member variables
		
		
		// The values of the matrix.  Basically the upper 3x3 portion
		// contains a linear transformation, and the last column is the
		// translation portion. Here data is transposed, check the header
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
			
			// Access to raw packed matrix data (usefull for
			// glLoadMatrixf () and glMultMatrixf ())
			Real m[16];
		};
	};
	// ===========================================================
	// Class typedefs
	// ===========================================================
	// Predefined Matrix4x4 types
	typedef Matrix4x4<float>    Matrix4x4f;
	typedef Matrix4x4<double>   Matrix4x4d;

	inline void PrintMatrix(const math::Matrix4x4f& matrix)
	{
		printf("%f %f %f %f\n", matrix[0], matrix[1], matrix[2], matrix[3]);
		printf("%f %f %f %f\n", matrix[4], matrix[5], matrix[6], matrix[7]);
		printf("%f %f %f %f\n", matrix[8], matrix[9], matrix[10], matrix[11]);
		printf("%f %f %f %f\n", matrix[12], matrix[13], matrix[14], matrix[15]);
	}
	
	// Nonmember Matrix4x4 functions
    
    // Matrix concatenation
	template <typename Real>
	Matrix4x4<Real> operator* (const Matrix4x4<Real>& m1, const Matrix4x4<Real>& m2);
    
    template <typename Real>
	Matrix4x4<Real> operator*= (const Matrix4x4<Real>& m1, const Matrix4x4<Real>& m2);
    
    template <typename Real>
    std::ostream& operator << (std::ostream& output, const Matrix4x4<Real>& m);
	
	// ===========================================================
	// Template/Inline implementation
	// ===========================================================
	template <typename Real>
	inline
	void Matrix4x4<Real>::operator= (const Matrix4x4<Real>& copy)
	{
		for(int i=0; i<16; ++i)
		{
			m[i] = copy.m[i];
		}
	}
	
	#include "matrix.inl"
}
}

#endif
