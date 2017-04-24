//
//  vector2.h
//
//  Based upon Mathlib.h -- Copyright (c) 2005-2006 David Henry
//
//  This code is licenced under the MIT license.
//
//  Modified by Jorge López González on 17/07/12.
//

#ifndef bitthemall_vector2_h
#define bitthemall_vector2_h

namespace dc
{
namespace math
{
    /////////////////////////////////////////////////////////////////////////////
    //
    // class Vector2<Real> - A simple 2D vector class.
    //
    /////////////////////////////////////////////////////////////////////////////
    
    template <typename Number>
    class Vector2
    {
    // Static operations for class methods
    public:
        static const Vector2<Number> Zero();
        static const Vector2<Number> One();
        
        static const Vector2<Number> Up();
        static const Vector2<Number> Down();
        
        static const Vector2<Number> Left();
        static const Vector2<Number> Right();

    // Constructors
    public:
        
        Vector2() : x (0), y (0) {}
        Vector2(Number x, Number y): x(x), y(y) {}
        
        // Copy constructor
        Vector2(const Vector2<Number>& copy): x(copy.x), y(copy.y) {}
        
    // Accessors
    public:
        // Accessor.  This allows to use the vector object
        // like an array of Real. For example:
        // Vector3<float> v (...);
        // float f = v[1]; // access to _y
        operator const Number *() const { return vec; }
       
    // Methods
    public:
        // Vector comparison
        const bool operator== (const Vector2<Number>& v) const;
        const bool operator!= (const Vector2<Number>& v) const;
        
        // Vector negation
        Vector2<Number> operator- () const;
        
        // Vector transformations
        Vector2<Number> operator+ (const Vector2<Number>& v) const;
        Vector2<Number> operator- (const Vector2<Number>& v) const;
        Vector2<Number> operator* (const Vector2<Number>& v) const;
        
        Vector2<Number> operator* (const Number s) const;
        Vector2<Number> operator/ (const Number s) const;
        
        // Combined assignment operators to conform to C notation convention
        Vector2<Number>& operator+= (const Vector2<Number>& v);
        Vector2<Number>& operator-= (const Vector2<Number>& v);
        Vector2<Number>& operator*= (const Vector2<Number>& v);
        
        Vector2<Number>& operator*= (const Number s);
        Vector2<Number>& operator/= (const Number s);
        
    public:
        Vector2<Number>& Half() const;
        Vector2<Number>& Double() const;

    // Fields
    public:
        union
        {
            struct
            {
                Number x, y;
            };
            
            struct
            {
                Number u, v;
            };
            
            struct
            {
                Number width, height;
            };
            
            Number vec[2];
        };
    };
    
    // Predefined Vector3 types
    typedef Vector2<float> Vector2f;
    typedef Vector2<float> UVf;
    
    typedef Vector2<double> Vector2d;
    typedef Vector2<double> UVd;
    
    typedef Vector2<int> Vector2i;

    #include "vector2.inl"
}
}

#endif
