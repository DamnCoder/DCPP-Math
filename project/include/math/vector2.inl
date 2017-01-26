/////////////////////////////////////////////////////////////////////////////
//
// class Vector2<Number> implementation.
//
/////////////////////////////////////////////////////////////////////////////

// --------------------------------------------------------------------------
// Vector2 operators
//
// Operator overloading for basic vector operations.
// --------------------------------------------------------------------------

// Vector comparison

template <typename Number>
inline
const bool
Vector2<Number>::operator== (const Vector2<Number> &v) const
{
    return Approximately(x, v.x) && Approximately(y, v.y);
}

template <typename Number>
inline
const bool
Vector2<Number>::operator!= (const Vector2<Number> &v) const
{
    return (!Approximately(x, v.x) || !Approximately(y, v.y));
}

// Vector negation
            
template <typename Number>
inline
Vector2<Number>
Vector2<Number>::operator- () const
{
    return Vector2<Number> (-x, -y);
}
            
// Vector transformations

template <typename Number>
inline
Vector2<Number>
Vector2<Number>::operator+ (const Vector2<Number> &v) const
{
    return Vector2<Number> (x + v.x, y + v.y);
}
            
template <typename Number>
inline
Vector2<Number>
Vector2<Number>::operator- (const Vector2<Number> &v) const
{
    return Vector2<Number> (x - v.x, y - v.y);
}
            
template <typename Number>
inline
Vector2<Number>
Vector2<Number>::operator* (const Vector2<Number>& v) const
{
    return Vector2<Number>(x * v.x, y * v.y);
}
            
template <typename Number>
inline
Vector2<Number>
Vector2<Number>::operator* (const Number s) const
{
    return Vector2<Number> (x * s, y * s);
}
            
template <typename Number>
inline
Vector2<Number>
Vector2<Number>::operator/ (const Number s) const
{
    const Number oneOverS = 1.0 / s; // Note: no check for divide by zero
    return Vector2<Number> (x * oneOverS, y * oneOverS);
}

// Combined assignment operators to conform to C notation convention

template <typename Number>
inline
Vector2<Number>&
Vector2<Number>::operator+= (const Vector2<Number> &v)
{
    x += v.x;
    y += v.y;
    return *this;
}

template <typename Number>
inline
Vector2<Number>&
Vector2<Number>::operator-= (const Vector2<Number> &v)
{
    x -= v.x;
    y -= v.y;
    return *this;
}

template <typename Number>
inline
Vector2<Number>&
Vector2<Number>::operator*= (const Vector2<Number> &v)
{
    x *= v.x;
    y *= v.y;
    return *this;
}

template <typename Number>
inline
Vector2<Number>&
Vector2<Number>::operator*= (const Number s)
{
    x *= s;
    y *= s;
    return *this;
}

template <typename Number>
inline
Vector2<Number>&
Vector2<Number>::operator/= (const Number s)
{
    Number oneOverS = 1.0 / s; // Note: no check for divide by zero!
    x *= oneOverS;
    y *= oneOverS;
    return *this;
}


// --------------------------------------------------------------------------
// Static operations
// --------------------------------------------------------------------------


template <typename Number>
inline
const Vector2<Number>
Vector2<Number>::Zero()
{
    return Vector2<Number>();
}

template <typename Number>
inline
const Vector2<Number>
Vector2<Number>::One()
{
    return Vector2<Number>(1, 1);
}

template <typename Number>
inline
const Vector2<Number>
Vector2<Number>::Up()
{
    return Vector2<Number>(0, 1);
}

template <typename Number>
inline
const Vector2<Number>
Vector2<Number>::Down()
{
    return Vector2<Number>(0, -1);
}

template <typename Number>
inline
const Vector2<Number>
Vector2<Number>::Left()
{
    return Vector2<Number>(-1, 0);
}

template <typename Number>
inline
const Vector2<Number>
Vector2<Number>::Right()
{
    return Vector2<Number>(1, 0);
}

template <typename Number>
inline
Vector2<Number>&
Vector2<Number>::Half() const
{
    return this * 0.5;
}

template <typename Number>
inline
Vector2<Number>&
Vector2<Number>::Double() const
{
    return this * 2.0;
}
