//
//  rect.h
//  Grid
//
//  Created by Jorge López on 10/10/15.
//  Copyright © 2015 Jorge López. All rights reserved.
//

#ifndef rect_h
#define rect_h

#include "vector2.h"

namespace dc
{
namespace math
{
    /**
     * Rect
     * 
     * The coordinate origin is set in upper left corner
     */
    template <typename Number>
    class Rect
    {
    // Static / Constants
    public:
        static const Vector2f W_LEFT_H_UP;
        static const Vector2f W_CENTER_H_UP;
        static const Vector2f W_RIGHT_H_UP;
        
        static const Vector2f W_LEFT_H_CENTER;
        static const Vector2f W_CENTER_H_CENTER;
        static const Vector2f W_RIGHT_H_CENTER;
        
        static const Vector2f W_LEFT_H_DOWN;
        static const Vector2f W_CENTER_H_DOWN;
        static const Vector2f W_RIGHT_H_DOWN;
        
        
    // Getters / Setters
    public:
        
    // Constructors
    public:
        Rect() {}
        
        Rect(const Vector2<Number>& position,
             const Vector2<Number>& dimensions):
            pivot(Vector2<Number>::Zero()),
            position(position),
            dimensions(dimensions) {}
        
        Rect(const Vector2<Number>& pivot,
             const Vector2<Number>& position,
             const Vector2<Number>& dimensions):
            pivot(pivot),
            position(position),
            dimensions(dimensions) {}
        
        Rect(Number x, Number y, Number width, Number height):
            pivot(Vector2<Number>::Zero()),
            position(Vector2<Number>(x, y)),
            dimensions(Vector2<Number>(width, height)) {}
        
    // Copy Constructor
    public:
        Rect(const Rect<Number>& copy): position(copy.position), dimensions(copy.dimensions) {}
        
        // Methods
    public:
        // Vector comparison
        const bool operator== (const Rect<Number>& r) const
        {
            return Center() == r.Center() & Dimensions() == r.Dimensions() & Pivot() == r.Pivot();
        }
        
        const bool operator!= (const Rect<Number>& r) const
        {
            return Center() != r.Center() | Dimensions() != r.Dimensions() | Pivot() != r.Pivot();
        }
        
    // Methods
    public:
        
        const Number Area() const { return dimensions.width * dimensions.height; }
        
        const Number Width() const { return dimensions.width; }
        const Number Height() const { return dimensions.height; }
        
        const Vector2<Number>& Pivot() const { return pivot; }
        
        void SetPivot(const Vector2<Number>& pivot)
        {
            this->pivot = pivot;
        }
        
        const Vector2<Number> Center() const { return position + (dimensions * pivot); }
        
        void SetPosition(const Vector2<Number>& position)
        {
            this->position = position;
        }
        
        const Vector2<Number>& Dimensions() const { return dimensions; }
        
        void SetDimensions(const Vector2<Number>& dimensions)
        {
            this->dimensions = dimensions;
        }
        
        void SetDimensions(const Number width, const Number height)
        {
            this->dimensions = Vector2<Number>(width, height);
        }
        
        const Vector2<Number> LeftUpperCorner() const
        {
            return position;
        }
        
        const Vector2<Number> LeftLowerCorner() const
        {
            return Vector2<Number>(position.x, position.y + dimensions.height);
        }
        
        const Vector2<Number> RightUpperCorner() const
        {
            return Vector2<Number>(position.x + dimensions.width, position.y);
        }
        
        const Vector2<Number> RightLowerCorner() const
        {
            return position + dimensions;
        }
        
    // Fields
    private:
        Vector2<Number> pivot;
        Vector2<Number> position;
        Vector2<Number> dimensions;
    };
    
    static const Vector2f W_LEFT_H_UP = Vector2f(0.0f, 0.0f);
    static const Vector2f W_CENTER_H_UP = Vector2f(0.5f, 0.0f);
    static const Vector2f W_RIGHT_H_UP = Vector2f(1.0f, 0.0f);
    
    static const Vector2f W_LEFT_H_CENTER = Vector2f(0.0f, 0.5f);
    static const Vector2f W_CENTER_H_CENTER = Vector2f(0.5f, 0.5f);
    static const Vector2f W_RIGHT_H_CENTER = Vector2f(1.0f, 0.5f);
    
    static const Vector2f W_LEFT_H_DOWN = Vector2f(0.0f, 1.0f);
    static const Vector2f W_CENTER_H_DOWN = Vector2f(0.5f, 1.0f);
    static const Vector2f W_RIGHT_H_DOWN = Vector2f(1.0f, 1.0f);
    
    // Predefined types
    typedef Rect<float>     Rectf;
    typedef Rect<double>    Rectd;
    typedef Rect<int>       Recti;
}
}

#endif /* rect_h */
