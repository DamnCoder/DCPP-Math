//
//  dcmath.h
//
//  A small group of utility functions to be used with real numbers.
//
//  Created by Jorge López González on 15/07/12.
//

#ifndef bitthemall_real_h
#define bitthemall_real_h

#include <cmath>

namespace dc
{
namespace math
{
    const double MO_EPSILON = 0.0001;
    
    // Declare a global constant for pi and a few multiples.
    const double MO_PI = 3.14159265358979323846264338327950288;
	const double MO_2PI = MO_PI * 2.0;
	const double MO_PI_OVER_2 = MO_PI / 2.0;
	const double MO_ONE_OVER_PI = 1.0 / MO_PI;
	const double MO_ONE_OVER_2PI = 1.0 / MO_2PI;
    const double MO_ONE_OVER_180 = 1.0 / 180.0;
	const double MO_PI_OVER_180 = MO_PI / 180.0;	// Grados a radianes
	const double MO_180_OVER_PI = 180.0 / MO_PI;	// Radianes a grados
    
    template <typename Real>
    inline
    const Real
    Sin (const Real rad)
    {
        return std::sin (rad);
    }
    
    template <typename Real>
    inline
    const Real
    Cos (const Real rad)
    {
        return std::cos (rad);
    }
    
    template <typename Real>
    inline
    const Real
    Atan2 (const Real y, const Real x)
    {
        return std::atan2 (y, x);
    }
    
    template <typename Real>
    inline
    const Real
    Asin (const Real x)
    {
        return std::asin (x);
    }
    
	// Convert between degrees and radians
	template <typename Real>
	inline 
    const Real 
    DegToRad (const Real deg) 
    { 
        return deg * MO_PI_OVER_180;
    }
    
    // Convert between radians and degrees
	template <typename Real>
	inline 
    const Real 
    RadToDeg (const Real rad)
    { 
        return rad * MO_180_OVER_PI;
    }

    // Clamps a value between min and max
    template<typename Real>
    inline
    const Real 
    Clamp(const Real value, const Real min, const Real max)
    {
        if(value > 1.0) return max;
        if(value < 0.0) return min;
        return value;
    }
    
    template<typename Real>
    inline
    const Real 
    LerpClamped(const Real from, const Real to, const Real perc)
    {
        return (from + Clamp<Real>(perc, 0.0, 1.0)*(to - from));
    }
    
    template<typename Real>
    inline
    const Real 
    Lerp(const Real from, const Real to, const Real perc)
    {
        return (from + perc*(to - from));
    }
    
    template<typename Real>
    inline
    Real Abs(Real x)
    {
        return (x > 0) ? x : -x;
    }
    
    template<typename Real>
    inline
    const bool
    Approximately(const Real a, const Real b, const Real epsilon)
    {
        return Abs<Real>(a - b) < epsilon;
    }
    
    template<typename Real>
    inline
    const bool
    Approximately(const Real a, const Real b)
    {
        return Abs<Real>(a - b) < MO_EPSILON;
    }
    
    template<typename Real>
    inline
    const Real Sign(const Real x)
    {
        if (0 < x) return 1.0;
        if (x < 0) return -1.0;
        return 0.0;
    }
    
    template<typename Real>
    inline
    const Real Normal(const Real a, const Real b, const Real c, const Real d)
    {
        return Sqrt(a * a + b * b + c * c + d * d);
    }
    
    template<typename Real>
    inline
    const Real Sqrt(const Real value)
    {
        return std::sqrt(value);
    }
    
    inline
    const int
    Half(const int value)
    {
        return value >> 1;
    }
    
    inline
    const int
    Double(const int value)
    {
        return 1 << value;
    }
    
    /*
     inline
     int Abs(int a) {
     int mask = (a >> (sizeof(int) * CHAR_BIT - 1));
     return (a + mask) ^ mask;
     }
     
    inline
    float Abs(float f) 
    {
        int i=((*(int*)&f)&0x7fffffff);
        return (*(float*)&i);
    }
    
    inline
    double Abs(double g)
    {
        unsigned long int *gg;
        gg=(unsigned long int*)&g;
        *(gg)&=9223372036854775807llu;
        return g;
    }
    
    inline
    float Neg(float f) 
    {
        int i=((*(int*)&f)^0x80000000);
        return (*(float*)&i);
    }
    
    inline
    int Sign(float f) 
    {
        return 1+(((*(int*)&f)>>31)<<1);
    }
    */
}
}

#endif
