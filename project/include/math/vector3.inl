/////////////////////////////////////////////////////////////////////////////
//
// class Vector3<Real> implementation.
//
/////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//
// Static operations
//
//////////////////////////////////////////////////////////////////////////////

template <typename Real>
inline
const Vector3<Real>
Vector3<Real>::PositiveEulerAngles(const Vector3<Real>& euler)
{
    //return Vector3<Real>(remainder(euler.x + 360.0, 360.0), remainder(euler.y + 360.0, 360.0), remainder(euler.z + 360.0, 360.0));
    Vector3<Real> positiveEuler;
    
    positiveEuler.x = euler.x;
    if (positiveEuler.x < 0.0)
        positiveEuler.x += 360.0;
   
    positiveEuler.y = euler.y;
    if (positiveEuler.y < 0.0)
        positiveEuler.y += 360.0;
    
    positiveEuler.z = euler.z;
    if (positiveEuler.z < 0.0)
        positiveEuler.z += 360.0;
    
    return positiveEuler;
}


// Decodifica las normales en forma polar y devuelve un vector normal
template <typename Real>
inline
Vector3<Real>
Vector3<Real>::FromPolarAngles(Real zenith, Real azimuth)
{
    Real lat = zenith * (DC_PI2)/255.0;
    Real lng = azimuth * (DC_PI2)/255.0;
    
    return Vector3<Real>(std::cos(lng) * std::sin(lat), std::sin(lng) * std::sin(lat), std::cos(lat));
}

// Devuelve el vector unidad que apunta en la orientacion pasada como argumento, en el plano X,Y
template <typename Real>
inline
Vector3<Real>
Vector3<Real>::UnitVectorFromOrientation(Real degrees)
{
    const float rads = DegToRad(degrees);
    return Vector3<Real>(std::sin(rads), -std::cos(rads), 0);
}

// Vector3 dot product
template <typename Real>
inline
Real
Vector3<Real>::DotProduct (const Vector3<Real> &a, const Vector3<Real> &b)
{
    return ((a.x * b.x) +  (a.y * b.y) +  (a.z * b.z));
}

// Vector3 cross product
template <typename Real>
inline
Vector3<Real>
Vector3<Real>::CrossProduct (const Vector3<Real> &a, const Vector3<Real> &b)
{
    return Vector3<Real> ((a.y * b.z) - (a.z * b.y),
                          (a.z * b.x) - (a.x * b.z),
                          (a.x * b.y) - (a.y * b.x));
}

// Compute normal plane given three points
template <typename Real>
inline
Vector3<Real>
Vector3<Real>::ComputeNormal (const Vector3<Real> &p1, const Vector3<Real> &p2, const Vector3<Real> &p3)
{
    Vector3<Real> vec1 (p1 - p2);
    Vector3<Real> vec2 (p1 - p3);
    
    Vector3<Real> result (CrossProduct (vec1, vec2));
    result.Normalize ();
    
    return result;
}

// Compute distance between two points
template <typename Real>
inline
Real
Vector3<Real>::Distance (const Vector3<Real>& from, const Vector3<Real>& to)
{
    Real dx = to.x - from.x;
    Real dy = to.y - from.y;
    Real dz = to.z - from.z;
    return std::sqrt ((dx * dx) + (dy * dy) + (dz * dz));
}

// Compute squared distance between two points.
// Useful when comparing distances, since we don't need
// to square the result.
template <typename Real>
inline
Real
Vector3<Real>::DistanceSquared (const Vector3<Real>& from, const Vector3<Real>& to)
{
    Real dx = to.x - from.x;
    Real dy = to.y - from.y;
    Real dz = to.z - from.z;
    return ((dx * dx) + (dy * dy) + (dz * dz));
}

/**
 * Interpolacion lineal de vectores en funcion de un parametro
 * porcentual t.
 * Siendo t el porcentaje de interpolacion entre un vector y el otro.
 * Si t = 0, no hay interpolacion, si t = 1 la interpolacion seria
 * el vector final
 */

template <typename Real>
inline
Vector3<Real>
Vector3<Real>::Lerp(const Vector3<Real>& from, const Vector3<Real>& to, const Real t)
{
    return (from + Clamp<float>(t, 0.0f, 1.0f)*(to - from));
}

//////////////////////////////////////////////////////////////////////////////
//
// Getters / Setters
//
//////////////////////////////////////////////////////////////////////////////

template <typename Real>
inline
void Vector3<Real>::SetCoords(const Real x, const Real y, const Real z)
{
    this->x = x; this->y = y; this->z = z;
}

// --------------------------------------------------------------------------
// Vector3::IsZero
//
// Return true if is zero vector.
// --------------------------------------------------------------------------
template <typename Real>
inline
const bool
Vector3<Real>::IsZero () const
{
    return (x == 0.0) && (y == 0.0) && (z == 0.0);
}

template <typename Real>
inline
const bool
Vector3<Real>::IsNearZero() const
{
    return Approximately(x, 0.0) && Approximately(y, 0.0) && Approximately(z, 0.0);
}

template <typename Real>
inline
const bool
Vector3<Real>::IsNearZero(const Real epsilon) const
{
    return Approximately(x, 0.0, epsilon) && Approximately(y, 0.0, epsilon) && Approximately(z, 0.0, epsilon);
}

template <typename Real>
inline
const bool
Vector3<Real>::IsApprox(const Vector3<Real>& v) const
{
    return Approximately(x, v.x) && Approximately(y, v.y) && Approximately(z, v.z);
}

template <typename Real>
inline
const bool
Vector3<Real>::IsApprox(const Vector3<Real>& v, const Real epsilon) const
{
    return Approximately(x, v.x, epsilon) && Approximately(y, v.y, epsilon) && Approximately(z, v.z, epsilon);
}

// --------------------------------------------------------------------------
// Vector3::Length
//
// Gets the vector magnitude
// --------------------------------------------------------------------------

template <typename Real>
inline
const Real
Vector3<Real>::Length() const
{
    return std::sqrt((x * x) + (y * y) + (z * z));
}

// --------------------------------------------------------------------------
// Vector3::normalize
//
// Set vector length to 1.
// --------------------------------------------------------------------------

template <typename Real>
inline
void
Vector3<Real>::Normalize()
{
    Real magSq = (x * x) + (y * y) + (z * z);
    
    if (magSq > 0.0)
    {
        Real oneOverMag = 1.0 / std::sqrt (magSq);
        x *= oneOverMag;
        y *= oneOverMag;
        z *= oneOverMag;
    }
}

//////////////////////////////////////////////////////////////////////////////
//
// Vector transformations
//
// Operator overloading for basic vector operations.
//
//////////////////////////////////////////////////////////////////////////////

// Vector comparison

template <typename Real>
inline
const bool
Vector3<Real>::operator== (const Vector3<Real> &v) const
{
    return Approximately(x, v.x) && Approximately(y, v.y) && Approximately(z, v.z);
}

template <typename Real>
inline
const bool
Vector3<Real>::operator!= (const Vector3<Real> &v) const
{
    return (!Approximately(x, v.x) || !Approximately(y, v.y) || !Approximately(z, v.z));
}

// Vector negation

template <typename Real>
inline
Vector3<Real>
Vector3<Real>::operator- () const
{
    return Vector3<Real> (-x, -y, -z);
}


template <typename Real>
inline
Vector3<Real>
Vector3<Real>::operator+ (const Vector3<Real> &v) const
{
    return Vector3<Real> (x + v.x, y + v.y, z + v.z);
}

template <typename Real>
inline
Vector3<Real>
Vector3<Real>::operator- (const Vector3<Real> &v) const
{
    return Vector3<Real> (x - v.x, y - v.y, z - v.z);
}

template <typename Real>
inline
Vector3<Real>
Vector3<Real>::operator* (const Vector3<Real>& v) const
{
    return Vector3<Real>(x * v.x, y * v.y, z * v.z);
}

template <typename Real>
inline
Vector3<Real>
Vector3<Real>::operator^ (const Vector3<Real>& v) const
{
    return Vector3<Real> ((y * v.z) - (z * v.y),
                          (z * v.x) - (x * v.z),
                          (x * v.y) - (y * v.x));
}

template <typename Real>
inline
Vector3<Real>
Vector3<Real>::operator* (const Real s) const
{
    return Vector3<Real> (x * s, y * s, z * s);
}

template <typename Real>
inline
Vector3<Real>
Vector3<Real>::operator/ (const Real s) const
{
    Real oneOverS = 1.0 / s; // Note: no check for divide by zero
    return Vector3<Real> (x * oneOverS, y * oneOverS, z * oneOverS);
}

// Combined assignment operators to conform to C notation convention

template <typename Real>
inline
Vector3<Real>&
Vector3<Real>::operator+= (const Vector3<Real> &v)
{
    x += v.x;
    y += v.y;
    z += v.z;
    return *this;
}

template <typename Real>
inline
Vector3<Real>&
Vector3<Real>::operator-= (const Vector3<Real> &v)
{
    x -= v.x;
    y -= v.y;
    z -= v.z;
    return *this;
}

template <typename Real>
inline
Vector3<Real>&
Vector3<Real>::operator*= (const Vector3<Real> &v)
{
    x *= v.x;
    y *= v.y;
    z *= v.z;
    return *this;
}

template <typename Real>
inline
Vector3<Real>&
Vector3<Real>::operator^= (const Vector3<Real> &v)
{
    // We can't modify the member atributes because with this operation we modify them and the result is affected
    *this = *this ^ v;
    return *this;
}

template <typename Real>
inline
Vector3<Real>&
Vector3<Real>::operator*= (const Real s)
{
    x *= s;
    y *= s;
    z *= s;
    return *this;
}

template <typename Real>
inline
Vector3<Real>&
Vector3<Real>::operator/= (const Real s)
{
    Real oneOverS = 1.0 / s; // Note: no check for divide by zero!
    x *= oneOverS;
    y *= oneOverS ;
    z *= oneOverS;
    return *this;
}

//////////////////////////////////////////////////////////////////////////////
//
// Methods
//
//////////////////////////////////////////////////////////////////////////////

/**
 * Situa el vector en el vertice más alejado
 * de la hipotenusa de un triangulo rectangulo
 * con catetos x = cos(degrees) y z = sin (degrees)
 *
 *		  Eje Z				 * (vector resultado)
 *			^			    /|
 *			|	          /  |
 *			|	  		/	 |
 *			|	      /		 |
 *			|	    /		 |	sin(degrees)
 *			|     /			 |
 *			|   /			 |
 *			| /	] degrees 	 |
 * (origen) *----------------|-----> Eje X
 *				cos(degrees)
 */

template <typename Real>
inline void
Vector3<Real>::RotateAndMoveXZ(const Real degrees, const Real distance)
{
    const float rads = DegToRad(degrees);
    x += std::cos(rads) * distance;
    z -= std::sin(rads) * distance;
}

template <typename Real>
inline void
Vector3<Real>::RotateAndMoveXY(const Real degrees, const Real distance)
{
    const float rads = DegToRad(degrees);
    x += std::cos(rads) * distance;
    y -= std::sin(rads) * distance;
}

template <typename Real>
inline void
Vector3<Real>::RotateAndMoveYZ(const Real degrees, const Real distance)
{
    const float rads = DegToRad(degrees);
    y += std::cos(rads) * distance;
    z -= std::sin(rads) * distance;
}

//////////////////////////////////////////////////////////////////////////////
//
// Nonmember methods
//
//////////////////////////////////////////////////////////////////////////////

// Scalar on the left multiplication, for symmetry
template <typename Real>
inline
Vector3<Real>
operator* (const Real k, const Vector3<Real>& v)
{
    return Vector3<Real> (k * v.x, k * v.y, k * v.z);
}
