//
//  vector3.h
//
//  Based upon Mathlib.h -- Copyright (c) 2005-2006 David Henry
//
//  This code is licenced under the MIT license.
//
//  Modified by Jorge López González on 14/07/12.
//

#ifndef bitthemall_vector_h
#define bitthemall_vector_h

#include <iostream>
#include "math.h"

namespace dc
{
namespace math
{
    //////////////////////////////////////////////////
	//
	// class Vector3<Real> - A simple 3D vector class.
	//
	//////////////////////////////////////////////////
    
	template <typename Real>
	class Vector3
	{
    // Static operations for class methods
    public:
        static const Vector3<Real> Zero();
        static const Vector3<Real> One();
        
        static const Vector3<Real> Up();
        static const Vector3<Real> Down();
        
        static const Vector3<Real> Left();
        static const Vector3<Real> Right();
        
        static const Vector3<Real> Forward();
        static const Vector3<Real> Backward();
        
        static const Vector3<Real> PositiveEulerAngles(const Vector3<Real>& euler);
        
        static Vector3<Real> FromPolarAngles(Real azimuth, Real zenith);
        static Vector3<Real> UnitVectorFromOrientation(Real degrees);
        
        static Real DotProduct (const Vector3<Real>& a, const Vector3<Real>& b);
        static Vector3<Real> CrossProduct (const Vector3<Real>& a, const Vector3<Real>& b);
        static Vector3<Real> ComputeNormal (const Vector3<Real>& p1, const Vector3<Real>& p2, const Vector3<Real>& p3);
        
        static Real Distance (const Vector3<Real>& from, const Vector3<Real>& to);
        static Real DistanceSquared (const Vector3<Real>& from, const Vector3<Real>& to);
        
        static Vector3<Real> Lerp(const Vector3<Real>& from, const Vector3<Real>& to, const Real perc);
        
    // Getters / Setters
    public:
        /*
         Accessor.  This allows to use the vector object
         like an array of Real. For example:
         Vector3<float> v (...);
         float f = v[1]; // access to _y
         */
        operator const Real *() const { return vec; }
        
        void SetCoords(const Real x, const Real y, const Real z);
        
        const bool IsZero () const;
        const bool IsNearZero() const;
        const bool IsNearZero(const Real epsilon) const;
        
        const bool IsApprox(const Vector3<Real>& v) const;
        const bool IsApprox(const Vector3<Real>& v, const Real epsilon) const;
        
        const Real Length() const;
        
        void Normalize ();


    // Constructors
	public:
        Vector3 () : x(0), y(0), z(0) {}
        Vector3 (const Real x, const Real y,const Real z = 0):x(x), y(y), z(z) {}
    
        // Copy constructor
        Vector3 (const Vector3<Real>& copy):x(copy.x), y(copy.y), z(copy.z) {}
        
    // Vector operations
	public:
        // Vector comparison
        const bool operator== (const Vector3<Real>& v) const;
        const bool operator!= (const Vector3<Real>& v) const;
        
        // Vector negation
        Vector3<Real> operator- () const;

        // Vector transformations
        Vector3<Real> operator+ (const Vector3<Real>& v) const;
        Vector3<Real> operator- (const Vector3<Real>& v) const;
        Vector3<Real> operator* (const Vector3<Real>& v) const;
        Vector3<Real> operator^ (const Vector3<Real>& v) const;
        
        Vector3<Real> operator* (const Real s) const;
        Vector3<Real> operator/ (const Real s) const;
        
        // Combined assignment operators to conform to C notation convention
        Vector3<Real>& operator+= (const Vector3<Real>& v);
        Vector3<Real>& operator-= (const Vector3<Real>& v);
        Vector3<Real>& operator*= (const Vector3<Real>& v);
        Vector3<Real>& operator^= (const Vector3<Real>& v);

        Vector3<Real>& operator*= (const Real s);
        Vector3<Real>& operator/= (const Real s);

    // Methods
    public:
        void RotateAndMoveXZ (Real degrees, Real distance);
        void RotateAndMoveXY (Real degrees, Real distance);
        void RotateAndMoveYZ (Real degrees, Real distance);

	public:
        // Member variables
        union
        {
            struct
            {
                Real x, y, z;
            };
            
            struct
            {
                Real r, g, b;
            };
            
            Real vec[3];
        };
	};

    // Predefined Vector3 types
	typedef Vector3<float> Vector3f;
	typedef Vector3<float> Position3f;
    typedef Vector3<float> ColorRGBf;
    
	typedef Vector3<double> Vector3d;
	typedef Vector3<double> Position3d;
    typedef Vector3<double> ColorRGBd;
    
	//
	// Nonmember Vector3 functions
	//
    
	template <typename Real>
	Vector3<Real>
    operator* (const Real k, const Vector3<Real>& v);
    
    template <typename Real>
    std::ostream&
    operator << (std::ostream &output, const Vector3<Real>& v);
    
    #include "vector3.inl"
}
}

#endif
