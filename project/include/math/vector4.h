//
//  vector4.h
//
//  Based upon Mathlib.h -- Copyright (c) 2005-2006 David Henry
//
//  This code is licenced under the MIT license.
//
//  Modified by Jorge López González on 17/07/12.
//

#ifndef bitthemall_vector4_h
#define bitthemall_vector4_h

#include "vector3.h"
#include "math.h"

/////////////////////////////////////////////////////////////////////////////
//
// class Vector4<Real> - A simple 4D vector class.
//
/////////////////////////////////////////////////////////////////////////////

namespace dc
{
namespace math
{

    template <typename Real>
    class Vector4
    {
	public:
		static Vector4<Real> Red()		{ return Vector4<Real>(1.0, 0.0, 0.0, 1.0); }
		static Vector4<Real> Green()	{ return Vector4<Real>(0.0, 1.0, 0.0, 1.0); }
		static Vector4<Real> Blue()		{ return Vector4<Real>(0.0, 0.0, 1.0, 1.0); }
		static Vector4<Real> White()	{ return Vector4<Real>(1.0, 1.0, 1.0, 1.0); }
		static Vector4<Real> Black()	{ return Vector4<Real>(0.0, 0.0, 0.0, 1.0); }
		
    public:
        Vector4() : x (0), y (0), z(0), w(0) {}
        Vector4(Real x, Real y, Real z, Real w)
            : x(x), y(y), z(z), w(w) {}
        
        // Copy constructor
        Vector4(const Vector4<Real>& copy): x(copy.x), y(copy.y), z(copy.z), w(copy.w) {}
        
    public:
        // Vector operations
        
        // Vector comparison
        const bool operator== (const Vector4<Real>& v) const;
        const bool operator!= (const Vector4<Real>& v) const;
        
        // Accessor.  This allows to use the vector object
        // like an array of Real. For example:
        // Vector3<float> v (...);
        // float f = v[1]; // access to _y
        operator const Real *() const { return vec; }
        
    public:
        void Normalize ();
        
        union
        {
            struct
            {
                Real x, y, z, w;
            };
            
            struct
            {
                Real r, g, b, a;
            };
            
            Real vec[4];
        };
    };

    // Predefined Vector3 types
    typedef Vector4<float> Vector4f;
    typedef Vector4<float> Planef;
    typedef Vector4<float> ColorRGBAf;

    typedef Vector4<double> Vector4d;
    typedef Vector4<double> Planed;
    typedef Vector4<double> ColorRGBAd;
    //
    // Nonmember Vector4 functions
    //
    
    template <typename Real>
    bool IsFacingTo (const Vector4<Real> plane, const Vector3<Real> &v);

    template <typename Real>
    bool IsFacingTo (const Vector4<Real> plane, const Vector3<Real> &v, Real radius);
    
    /////////////////////////////////////////////////////////////////////////////
    //
    // class Vector4<Real> implementation.
    //
    /////////////////////////////////////////////////////////////////////////////
    
    // --------------------------------------------------------------------------
    // Vector4 operators
    //
    // Operator overloading for basic vector operations.
    // --------------------------------------------------------------------------
    
    // Vector comparison
    
    template <typename Real>
    inline
    const bool
    Vector4<Real>::operator== (const Vector4<Real> &v) const
    {
        return Approximately(x, v.x) && Approximately(y, v.y) && Approximately(z, v.z) && Approximately(w, v.w);
    }
    
    template <typename Real>
    inline
    const bool
    Vector4<Real>::operator!= (const Vector4<Real> &v) const
    {
        return (!Approximately(x, v.x) || !Approximately(y, v.y) || !Approximately(z, v.z) || !Approximately(w, v.w));
    }
        
    // Decide if a plane is facing a point
    template <typename Real>
    inline
    bool
    IsFacingTo (const Vector4<Real> plane, const Vector3<Real> &v)
    {
        return ((plane.x*v.x) + (plane.y*v.y) + (plane.z*v.z) + (plane.w))>0;
    }
    
    // Decide if a plane is facing a sphere
    template <typename Real>
    inline
    bool
    IsFacingTo (const Vector4<Real> plane, const Vector3<Real> &v, Real radius)
    {
        return ((plane.x*v.x) + (plane.y*v.y) + (plane.z*v.z) + (plane.w))>-radius;
    }

}
}

#endif
