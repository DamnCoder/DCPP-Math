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
//                | +x        +y' |
//               /                            /
//              /                            /
//             +z                          +z'
//
//          initial position      Positive rotation of
//                                 pi/2 around z-axis
//
/////////////////////////////////////////////////////////////////////////////


template <typename Real>
inline Matrix4x4<Real>
Matrix4x4<Real>::operator +=(const Matrix4x4<Real>& m)
{
    m11 += m.m11;   m12 += m.m12;   m13 += m.m13;
    m21 += m.m21;   m22 += m.m22;   m23 += m.m23;
    m31 += m.m31;   m32 += m.m32;   m33 += m.m33;
    
    tx += m.tx; ty += m.ty; tz += m.tz;
    
    return *this;
}

template <typename Real>
inline Matrix4x4<Real>
Matrix4x4<Real>::operator -=(const Matrix4x4<Real>& m)
{
    m11 -= m.m11; m12 -= m.m12; m13 -= m.m13;
    m21 -= m.m21; m22 -= m.m22; m23 -= m.m23;
    m31 -= m.m31; m32 -= m.m32; m33 -= m.m33;
    
    tx -= m.tx; ty -= m.ty; tz -= m.tz;
    
    return *this;
}

// --------------------------------------------------------------------------
// Matrix4x4::multiplication / concatenation
//
// Concatenación de dos matrices.
//
// Debido a que estamos usando una version traspuesta de la matriz
// que se usa en OpenGL la multiplicación la hacemos columna por fila,
// en lugar de fila por columna.
// --------------------------------------------------------------------------

template <typename Real>
inline Matrix4x4<Real>
Matrix4x4<Real>::operator *=(const Matrix4x4<Real> m)
{
    // Compute the left 4x3 (linear transformation) portion
    m11 = (m11 * m.m11) + (m21 * m.m12) + (m31 * m.m13);
    m12 = (m12 * m.m11) + (m22 * m.m12) + (m32 * m.m13);
    m13 = (m13 * m.m11) + (m23 * m.m12) + (m33 * m.m13);
    
    m21 = (m11 * m.m21) + (m21 * m.m22) + (m31 * m.m23);
    m22 = (m12 * m.m21) + (m22 * m.m22) + (m32 * m.m23);
    m23 = (m13 * m.m21) + (m23 * m.m22) + (m33 * m.m23);
    
    m31 = (m11 * m.m31) + (m21 * m.m32) + (m31 * m.m33);
    m32 = (m12 * m.m31) + (m22 * m.m32) + (m32 * m.m33);
    m33 = (m13 * m.m31) + (m23 * m.m32) + (m33 * m.m33);
    
    // Compute the translation portion
    tx = (m11 * b.tx) + (m21 * b.ty) + (m31 * b.tz) + tx;
    ty = (m12 * b.tx) + (m22 * b.ty) + (m32 * b.tz) + ty;
    tz = (m13 * b.tx) + (m23 * b.ty) + (m33 * b.tz) + tz;
    
    return *this;
}

template <typename Real>
inline void
Vector3<Real>::operator *= (const Vector3<Real>& v)
{
    // Grind through the linear algebra.
    v = Vector3<Real> ((v.x * m11) + (v.y * m21) + (v.z * m31) + tx,
                       (v.x * m12) + (v.y * m22) + (v.z * m32) + ty,
                       (v.x * m13) + (v.y * m23) + (v.z * m33) + tz );
}

// --------------------------------------------------------------------------
// Matrix4x4::identity
//
// Set matrix to identity.
// --------------------------------------------------------------------------

template <typename Real>
inline void
Matrix4x4<Real>::identity ()
{
    m11 = 1.0; m21 = 0.0; m31 = 0.0; tx = 0.0;
    m12 = 0.0; m22 = 1.0; m32 = 0.0; ty = 0.0;
    m13 = 0.0; m23 = 0.0; m33 = 1.0; tz = 0.0;
    h14 = 0.0; h24 = 0.0; h34 = 0.0; tw = 1.0;
}


// --------------------------------------------------------------------------
// Matrix4x4::transpose
//
// Transpose the current matrix.
// --------------------------------------------------------------------------

template <typename Real>
inline void
Matrix4x4<Real>::transpose ()
{
    *this = Transpose (*this);
}


// --------------------------------------------------------------------------
// Matrix4x4::invert
//
// Invert the current matrix.
// --------------------------------------------------------------------------

template <typename Real>
inline void
Matrix4x4<Real>::invert ()
{
    *this = Invert (*this);
}


// --------------------------------------------------------------------------
// Matrix4x4::setTranslation
//
// Set the translation portion of the matrix.
// --------------------------------------------------------------------------

template <typename Real>
inline void
Matrix4x4<Real>::setTranslation (const Vector3<Real> &v)
{
    tx = v.x; ty = v.y; tz = v.z;
}


// --------------------------------------------------------------------------
// Matrix4x4::transform
//
// Transform a point by the matrix.
// --------------------------------------------------------------------------

template <typename Real>
inline void
Matrix4x4<Real>::transform (Vector3<Real> &v) const
{
    // Grind through the linear algebra.
    v = Vector3<Real> (
                       (v.x * m11) + (v.y * m21) + (v.z * m31) + tx,
                       (v.x * m12) + (v.y * m22) + (v.z * m32) + ty,
                       (v.x * m13) + (v.y * m23) + (v.z * m33) + tz
                       );
}


// --------------------------------------------------------------------------
// Matrix4x4::rotate
//
// Rotate a point by the 3x3 upper left portion of the matrix.
// --------------------------------------------------------------------------

template <typename Real>
inline void
Matrix4x4<Real>::rotate (Vector3<Real> &v) const
{
    v = Vector3<Real> (
                       (v.x * m11) + (v.y * m21) + (v.z * m31),
                       (v.x * m12) + (v.y * m22) + (v.z * m32),
                       (v.x * m13) + (v.y * m23) + (v.z * m33)
                       );
}


// --------------------------------------------------------------------------
// Matrix4x4::inverseRotate
//
// Rotate a point by the inverse 3x3 upper left portion of the matrix.
// --------------------------------------------------------------------------

template <typename Real>
inline void
Matrix4x4<Real>::inverseRotate (Vector3<Real> &v) const
{
    v = Vector3<Real> (
                       (v.x * m11) + (v.y * m12) + (v.z * m13),
                       (v.x * m21) + (v.y * m22) + (v.z * m23),
                       (v.x * m31) + (v.y * m32) + (v.z * m33)
                       );
}


// --------------------------------------------------------------------------
// Matrix4x4::inverseRotate
//
// Translate a point by the inverse matrix.
// --------------------------------------------------------------------------

template <typename Real>
inline void
Matrix4x4<Real>::inverseTranslate (Vector3<Real> &v) const
{
    v.x -= tx;
    v.y -= ty;
    v.z -= tz;
}


// --------------------------------------------------------------------------
// Matrix4x4::fromQuaternion
//
// Convert a quaternion to a matrix.
// --------------------------------------------------------------------------

template <typename Real>
inline void
Matrix4x4<Real>::fromQuaternion (const Quaternion<Real> &q)
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
Matrix4x4<Real>::fromEulerAngles (Real x, Real y, Real z)
{
    // Fetch sine and cosine of angles
    Real cx = std::cos (x);
    Real sx = std::sin (x);
    Real cy = std::cos (y);
    Real sy = std::sin (y);
    Real cz = std::cos (z);
    Real sz = std::sin (z);
    
    Real sxsy = sx * sy;
    Real cxsy = cx * sy;
    
    // Fill in the matrix elements
    m11 =  (cy * cz);
    m12 =  (sxsy * cz) + (cx * sz);
    m13 = -(cxsy * cz) + (sx * sz);
    
    m21 = -(cy * sz);
    m22 = -(sxsy * sz) + (cx * cz);
    m23 =  (cxsy * sz) + (sx * cz);
    
    m31 =  (sy);
    m32 = -(sx * cy);
    m33 =  (cx * cy);
    
    // Reset the translation portion
    tx = ty = tz = 0.0;
}


// --------------------------------------------------------------------------
// Matrix4x4::toEulerAngles
//
// Setup the euler angles in radians, given a rotation matrix. The rotation
// matrix could have been obtained from euler angles given the expression:
//   M = X.Y.Z
// where X, Y and Z are rotation matrices about X, Y and Z axes.
// This is the "opposite" of the fromEulerAngles function.
// --------------------------------------------------------------------------

template <typename Real>
inline void
Matrix4x4<Real>::toEulerAngles (Real &x, Real &y, Real &z) const
{
    // Compute Y-axis angle
    y = std::asin (m31);
    
    // Compute cos and one over cos for optimization
    Real cy = std::cos (y);
    Real oneOverCosY = 1.0 / cy;
    
    if (std::fabs (cy) > 0.001)
    {
        // No gimball lock
        x = std::atan2 (-m32 * oneOverCosY, m33 * oneOverCosY);
        z = std::atan2 (-m21 * oneOverCosY, m11 * oneOverCosY);
    }
    else
    {
        // Gimbal lock case
        x = 0.0;
        z = std::atan2 (m12, m22);
    }
}


// --------------------------------------------------------------------------
// Matrix4x4::rightVector
// Matrix4x4::upVector
// Matrix4x4::forwardVector
// Matrix4x4::translationVector
//
// Return a base vector from the matrix.
// --------------------------------------------------------------------------

template <typename Real>
inline Vector3<Real>
Matrix4x4<Real>::rightVector () const
{
    return Vector3<Real> (m11, m12, m13);
}

template <typename Real>
inline Vector3<Real>
Matrix4x4<Real>::upVector () const
{
    return Vector3<Real> (m21, m22, m23);
}

template <typename Real>
inline Vector3<Real>
Matrix4x4<Real>::forwardVector () const
{
    return Vector3<Real> (m31, m32, m33);
}

template <typename Real>
inline Vector3<Real>
Matrix4x4<Real>::translationVector () const
{
    return Vector3<Real> (tx, ty, tz);
}


// --------------------------------------------------------------------------
//
// Nonmember Matrix4x4<Real> functions
//
// --------------------------------------------------------------------------


// --------------------------------------------------------------------------
// Matrix4x4 * Matrix4x4
//
// Matrix concatenation.
//
// We also provide a *= operator, as per C convention.
// --------------------------------------------------------------------------

template <typename Real>
inline Matrix4x4<Real>
operator* (const Matrix4x4<Real> &a, const Matrix4x4<Real> &b)
{
    Matrix4x4<Real> res;
    
    // Compute the left 4x3 (linear transformation) portion
    res.m11 = (a.m11 * b.m11) + (a.m21 * b.m12) + (a.m31 * b.m13);
    res.m12 = (a.m12 * b.m11) + (a.m22 * b.m12) + (a.m32 * b.m13);
    res.m13 = (a.m13 * b.m11) + (a.m23 * b.m12) + (a.m33 * b.m13);
    
    res.m21 = (a.m11 * b.m21) + (a.m21 * b.m22) + (a.m31 * b.m23);
    res.m22 = (a.m12 * b.m21) + (a.m22 * b.m22) + (a.m32 * b.m23);
    res.m23 = (a.m13 * b.m21) + (a.m23 * b.m22) + (a.m33 * b.m23);
    
    res.m31 = (a.m11 * b.m31) + (a.m21 * b.m32) + (a.m31 * b.m33);
    res.m32 = (a.m12 * b.m31) + (a.m22 * b.m32) + (a.m32 * b.m33);
    res.m33 = (a.m13 * b.m31) + (a.m23 * b.m32) + (a.m33 * b.m33);
    
    // Compute the translation portion
    res.tx = (a.m11 * b.tx) + (a.m21 * b.ty) + (a.m31 * b.tz) + a.tx;
    res.ty = (a.m12 * b.tx) + (a.m22 * b.ty) + (a.m32 * b.tz) + a.ty;
    res.tz = (a.m13 * b.tx) + (a.m23 * b.ty) + (a.m33 * b.tz) + a.tz;
    
    return res;
}

template <typename Real>
inline Matrix4x4<Real> &
operator*= (Matrix4x4<Real> &a, const Matrix4x4<Real> &b)
{
    a = a * b;
    return a;
}


// --------------------------------------------------------------------------
// Matrix4x4 * Vector3
//
// Transform a point by a matrix.  This makes using the vector class look
// like it does with linear algebra notation on paper.
// --------------------------------------------------------------------------

template <typename Real>
inline Vector3<Real>
operator* (const Matrix4x4<Real> &m, const Vector3<Real> &p)
{
    Vector3<Real> res (p);
    m.transform (res);
    return res;
}


// --------------------------------------------------------------------------
// Transpose
//
// Return the transpose matrix.
// --------------------------------------------------------------------------

template <typename Real>
inline Matrix4x4<Real>
Transpose (const Matrix4x4<Real> &m)
{
    Matrix4x4<Real> res;
    
    res.m11 = m.m11; res.m21 = m.m12; res.m31 = m.m13; res.tx = m.h14;
    res.m12 = m.m21; res.m22 = m.m22; res.m32 = m.m23; res.ty = m.h24;
    res.m13 = m.m31; res.m23 = m.m32; res.m33 = m.m33; res.tz = m.h34;
    res.h14 = m.tx;  res.h24 = m.ty;  res.h34 = m.tz;  res.tw = m.tw;
    
    return res;
}


// --------------------------------------------------------------------------
// Determinant
//
// Compute the determinant of the 3x3 portion of the matrix.
// --------------------------------------------------------------------------

template <typename Real>
inline static Real
Determinant3x3 (const Matrix4x4<Real> &m)
{
    return m.m11 * ((m.m22 * m.m33) - (m.m23 * m.m32))
    + m.m12 * ((m.m23 * m.m31) - (m.m21 * m.m33))
    + m.m13 * ((m.m21 * m.m32) - (m.m22 * m.m31));
}


// --------------------------------------------------------------------------
// Invert
//
// Compute the inverse of a matrix.  We use the classical adjoint divided
// by the determinant method.
// --------------------------------------------------------------------------

template <typename Real>
inline Matrix4x4<Real>
Invert (const Matrix4x4<Real> &m)
{
    // Compute the determinant of the 3x3 portion
    Real det = Determinant3x3 (m);
    
    // If we're singular, then the determinant is zero and there's
    // no inverse
    assert (std::fabs (det) > 0.000001);
    
    // Compute one over the determinant, so we divide once and
    // can *multiply* per element
    Real oneOverDet = 1.0 / det;
    
    // Compute the 3x3 portion of the inverse, by
    // dividing the adjoint by the determinant
    Matrix4x4<Real> res;
    
    res.m11 = ((m.m22 * m.m33) - (m.m23 * m.m32)) * oneOverDet;
    res.m12 = ((m.m13 * m.m32) - (m.m12 * m.m33)) * oneOverDet;
    res.m13 = ((m.m12 * m.m23) - (m.m13 * m.m22)) * oneOverDet;
    
    res.m21 = ((m.m23 * m.m31) - (m.m21 * m.m33)) * oneOverDet;
    res.m22 = ((m.m11 * m.m33) - (m.m13 * m.m31)) * oneOverDet;
    res.m23 = ((m.m13 * m.m21) - (m.m11 * m.m23)) * oneOverDet;
    
    res.m31 = ((m.m21 * m.m32) - (m.m22 * m.m31)) * oneOverDet;
    res.m32 = ((m.m12 * m.m31) - (m.m11 * m.m32)) * oneOverDet;
    res.m33 = ((m.m11 * m.m22) - (m.m12 * m.m21)) * oneOverDet;
    
    // Compute the translation portion of the inverse
    res.tx = -((m.tx * res.m11) + (m.ty * res.m21) + (m.tz * res.m31));
    res.ty = -((m.tx * res.m12) + (m.ty * res.m22) + (m.tz * res.m32));
    res.tz = -((m.tx * res.m13) + (m.ty * res.m23) + (m.tz * res.m33));
    
    // Return it.
    return res;
}


// --------------------------------------------------------------------------
// RotationMatrix
//
// Setup the matrix to perform a rotation about one of the three cardinal
// X-Y-Z axes.
//
// The axis of rotation is specified by the 1-based "axis" index.
//
// theta is the amount of rotation, in radians.  The right-hand rule is
// used to define "positive" rotation.
//
// The translation portion is reset.
// --------------------------------------------------------------------------

template <typename Real>
inline Matrix4x4<Real>
RotationMatrix (Axis axis, Real theta)
{
    Matrix4x4<Real> res;
    
    // Get sin and cosine of rotation angle
    Real s = std::sin (theta);
    Real c = std::cos (theta);
    
    // Check which axis they are rotating about
    switch (axis)
    {
	    case kXaxis: // Rotate about the x-axis
            res.m11 = 1.0; res.m21 = 0.0; res.m31 = 0.0;
            res.m12 = 0.0; res.m22 = c;   res.m32 = -s;
            res.m13 = 0.0; res.m23 = s;   res.m33 =  c;
            break;
            
	    case kYaxis: // Rotate about the y-axis
            res.m11 = c;   res.m21 = 0.0; res.m31 = s;
            res.m12 = 0.0; res.m22 = 1.0; res.m32 = 0.0;
            res.m13 = -s;  res.m23 = 0.0; res.m33 = c;
            break;
            
	    case kZaxis: // Rotate about the z-axis
            res.m11 = c;   res.m21 = -s;  res.m31 = 0.0;
            res.m12 = s;   res.m22 =  c;  res.m32 = 0.0;
            res.m13 = 0.0; res.m23 = 0.0; res.m33 = 1.0;
            break;
            
	    default:
            // bogus axis index
            assert (false);
    }
    
    // Reset the translation portion
    res.tx = res.ty = res.tz = 0.0;
    
    return res;
}


//---------------------------------------------------------------------------
// AxisRotationMatrix
//
// Setup the matrix to perform a rotation about an arbitrary axis.
// The axis of rotation must pass through the origin.
//
// axis defines the axis of rotation, and must be a unit vector.
//
// theta is the amount of rotation, in radians.  The right-hand rule is
// used to define "positive" rotation.
//
// The translation portion is reset.
// --------------------------------------------------------------------------

template <typename Real>
inline Matrix4x4<Real>
RotationMatrix (const Vector3<Real> &axis, Real theta)
{
    Matrix4x4<Real> res;
    
    // Quick sanity check to make sure they passed in a unit vector
    // to specify the axis
    assert (std::fabs (DotProduct (axis, axis) - 1.0) < 0.001);
    
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
    res.m11 = (ax * axis.x) + c;
    res.m12 = (ax * axis.y) + (axis.z * s);
    res.m13 = (ax * axis.z) - (axis.y * s);
    
    res.m21 = (ay * axis.x) - (axis.z * s);
    res.m22 = (ay * axis.y) + c;
    res.m23 = (ay * axis.z) + (axis.x * s);
    
    res.m31 = (az * axis.x) + (axis.y * s);
    res.m32 = (az * axis.y) - (axis.x * s);
    res.m33 = (az * axis.z) + c;
    
    // Reset the translation portion
    res.tx = res.ty = res.tz = 0.0;
    
    return res;
}


// --------------------------------------------------------------------------
// TranslationMatrix
//
// Build a translation matrix given a translation vector.
// --------------------------------------------------------------------------

template <typename Real>
inline Matrix4x4<Real>
TranslationMatrix (Real x, Real y, Real z)
{
    return TranslationMatrix (Vector3<Real> (x, y, z));
}

template <typename Real>
inline Matrix4x4<Real>
TranslationMatrix (const Vector3<Real> &v)
{
    Matrix4x4<Real> res;
    
    res.m11 = 1.0; res.m21 = 0.0; res.m31 = 0.0; res.tx = v.x;
    res.m12 = 0.0; res.m22 = 1.0; res.m32 = 0.0; res.ty = v.y;
    res.m13 = 0.0; res.m23 = 0.0; res.m33 = 1.0; res.tz = v.z;
    
    return res;
}


// --------------------------------------------------------------------------
// ScaleMatrix
//
// Setup the matrix to perform scale on each axis.  For uniform scale by k,
// use a vector of the form Vector3( k, k, k ).
//
// The translation portion is reset.
// --------------------------------------------------------------------------

template <typename Real>
inline Matrix4x4<Real>
ScaleMatrix (const Vector3<Real> &s)
{
    Matrix4x4<Real> res;
    
    // Set the matrix elements.  Pretty straightforward
    res.m11 = s.x; res.m21 = 0.0;  res.m31 = 0.0;
    res.m12 = 0.0;  res.m22 = s.y; res.m32 = 0.0;
    res.m13 = 0.0;  res.m23 = 0.0;  res.m33 = s.z;
    
    // Reset the translation portion
    res.tx = res.ty = res.tz = 0.0;
    
    return res;
}


// --------------------------------------------------------------------------
// ScaleAlongAxisMatrix
//
// Setup the matrix to perform scale along an arbitrary axis.
//
// The axis is specified using a unit vector.
//
// The translation portion is reset.
// --------------------------------------------------------------------------

template <typename Real>
inline Matrix4x4<Real>
ScaleAlongAxisMatrix (const Vector3<Real> &axis, Real k)
{
    Matrix4x4<Real> res;
    
    // Quick sanity check to make sure they passed in a unit vector
    // to specify the axis
    assert (std::fabs (DotProduct (axis, axis) - 1.0) < 0.001);
    
    // Compute k-1 and some common subexpressions
    Real a = k - 1.0;
    Real ax = a * axis.x;
    Real ay = a * axis.y;
    Real az = a * axis.z;
    
    // Fill in the matrix elements.  We'll do the common
    // subexpression optimization ourselves here, since diagonally
    // opposite matrix elements are equal
    res.m11 = (ax * axis.x) + 1.0;
    res.m22 = (ay * axis.y) + 1.0;
    res.m32 = (az * axis.z) + 1.0;
    
    res.m12 = res.m21 = (ax * axis.y);
    res.m13 = res.m31 = (ax * axis.z);
    res.m23 = res.m32 = (ay * axis.z);
    
    // Reset the translation portion
    res.tx = res.ty = res.tz = 0.0;
    
    return res;
}


// --------------------------------------------------------------------------
// ShearMatrix
//
// Setup the matrix to perform a shear
//
// The type of shear is specified by the 1-based "axis" index.  The effect
// of transforming a point by the matrix is described by the pseudocode
// below:
//
//	xAxis  =>  y += s * x, z += t * x
//	yAxis  =>  x += s * y, z += t * y
//	zAxis  =>  x += s * z, y += t * z
//
// The translation portion is reset.
// --------------------------------------------------------------------------

template <typename Real>
inline Matrix4x4<Real>
ShearMatrix (Axis axis, Real s, Real t)
{
    Matrix4x4<Real> res;
    
    // Check which type of shear they want
    switch (axis)
    {
	    case kXaxis: // Shear y and z using x
            res.m11 = 1.0; res.m21 = 0.0; res.m31 = 0.0;
            res.m12 = s;   res.m22 = 1.0; res.m32 = 0.0;
            res.m13 = t;   res.m23 = 0.0; res.m33 = 1.0;
            break;
            
	    case kYaxis: // Shear x and z using y
            res.m11 = 1.0; res.m21 = s;   res.m31 = 0.0;
            res.m12 = 0.0; res.m22 = 1.0; res.m32 = 0.0;
            res.m13 = 0.0; res.m23 = t;   res.m33 = 1.0;
            break;
            
	    case kZaxis: // Shear x and y using z
            res.m11 = 1.0; res.m21 = 0.0; res.m31 = s;
            res.m12 = 0.0; res.m22 = 1.0; res.m32 = t;
            res.m13 = 0.0; res.m23 = 0.0; res.m33 = 1.0;
            break;
            
	    default:
            // bogus axis index
            assert (false);
    }
    
    // Reset the translation portion
    res.tx = res.ty = res.tz = 0.0;
    
    return res;
}


// --------------------------------------------------------------------------
// ProjectionMatrix
//
// Setup the matrix to perform a projection onto a plane passing
// through the origin.  The plane is perpendicular to the
// unit vector n.
// --------------------------------------------------------------------------

template <typename Real>
inline Matrix4x4<Real>
ProjectionMatrix (const Vector3<Real> &n)
{
    Matrix4x4<Real> res;
    
    // Quick sanity check to make sure they passed in a unit vector
    // to specify the axis
    assert (std::fabs (DotProduct (n, n) - 1.0) < 0.001);
    
    // Fill in the matrix elements.  We'll do the common
    // subexpression optimization ourselves here, since diagonally
    // opposite matrix elements are equal
    res.m11 = 1.0 - (n.x * n.x);
    res.m22 = 1.0 - (n.y * n.y);
    res.m33 = 1.0 - (n.z * n.z);
    
    res.m12 = res.m21 = -(n.x * n.y);
    res.m13 = res.m31 = -(n.x * n.z);
    res.m23 = res.m32 = -(n.y * n.z);
    
    // Reset the translation portion
    res.tx = res.ty = res.tz = 0.0;
    
    return res;
}


// --------------------------------------------------------------------------
// ReflectionMatrix
//
// Setup the matrix to perform a reflection about a plane parallel
// to a cardinal plane.
//
// axis is a 1-based index which specifies the plane to project about:
//
//	xAxis => reflect about the plane x=k
//	yAxis => reflect about the plane y=k
//	zAxis => reflect about the plane z=k
//
// The translation is set appropriately, since translation must occur if
// k != 0
// --------------------------------------------------------------------------

template <typename Real>
inline Matrix4x4<Real>
ReflectionMatrix (Axis axis, Real k)
{
    Matrix4x4<Real> res;
    
    // Check which plane they want to reflect about
    switch (axis)
    {
	    case kXaxis: // Reflect about the plane x=k
            res.m11 = -1.0; res.m21 =  0.0; res.m31 =  0.0; res.tx = 2.0 * k;
            res.m12 =  0.0; res.m22 =  1.0; res.m32 =  0.0; res.ty = 0.0;
            res.m13 =  0.0; res.m23 =  0.0; res.m33 =  1.0; res.tz = 0.0;
            break;
            
	    case kYaxis: // Reflect about the plane y=k
            res.m11 =  1.0; res.m21 =  0.0; res.m31 =  0.0; res.tx = 0.0;
            res.m12 =  0.0; res.m22 = -1.0; res.m32 =  0.0; res.ty = 2.0 * k;
            res.m13 =  0.0; res.m23 =  0.0; res.m33 =  1.0; res.tz = 0.0;
            break;
            
	    case kZaxis: // Reflect about the plane z=k
            res.m11 =  1.0; res.m21 =  0.0; res.m31 =  0.0; res.tx = 0.0;
            res.m12 =  0.0; res.m22 =  1.0; res.m32 =  0.0; res.ty = 0.0;
            res.m13 =  0.0; res.m23 =  0.0; res.m33 = -1.0; res.tz = 2.0 * k;
            break;
            
	    default:
            // bogus axis index
            assert (false);
    }
    
    return res;
}


// --------------------------------------------------------------------------
// AxisReflectionMatrix
//
// Setup the matrix to perform a reflection about an arbitrary plane
// through the origin.  The unit vector n is perpendicular to the plane.
//
// The translation portion is reset.
// --------------------------------------------------------------------------

template <typename Real>
inline Matrix4x4<Real>
AxisReflectionMatrix (const Vector3<Real> &n)
{
    Matrix4x4<Real> res;
    
    // Quick sanity check to make sure they passed in a unit vector
    // to specify the axis
    assert (std::fabs (DotProduct (n, n) - 1.0) < 0.001);
    
    // Compute common subexpressions
    Real ax = -2.0 * n.x;
    Real ay = -2.0 * n.y;
    Real az = -2.0 * n.z;
    
    // Fill in the matrix elements.  We'll do the common
    // subexpression optimization ourselves here, since diagonally
    // opposite matrix elements are equal
    res.m11 = 1.0 + (ax * n.x);
    res.m22 = 1.0 + (ay * n.y);
    res.m32 = 1.0 + (az * n.z);
    
    res.m12 = res.m21 = (ax * n.y);
    res.m13 = res.m31 = (ax * n.z);
    res.m23 = res.m32 = (ay * n.z);
    
    // Reset the translation portion
    res.tx = res.ty = res.tz = 0.00;
    
    return res;
}


// --------------------------------------------------------------------------
// LookAtMatrix
//
// Setup the matrix to perform a "Look At" transformation like a first
// person camera.
// --------------------------------------------------------------------------
template <typename Real>
inline Matrix4x4<Real>
LookAtMatrix (const Vector3<Real> &camPos, const Vector3<Real> &target,
		      const Vector3<Real> &camUp)
{
    Matrix4x4<Real> rot, trans;
    
    Vector3<Real> forward (camPos - target);
    forward.normalize ();
    
    Vector3<Real> right (CrossProduct (camUp, forward));
    Vector3<Real> up (CrossProduct (forward, right));
    
    right.normalize ();
    up.normalize ();
    
    rot.m11 = right.x;
    rot.m21 = right.y;
    rot.m31 = right.z;
    
    rot.m12 = up.x;
    rot.m22 = up.y;
    rot.m32 = up.z;
    
    rot.m13 = forward.x;
    rot.m23 = forward.y;
    rot.m33 = forward.z;
    
    rot.tx  = 0.0;
    rot.ty  = 0.0;
    rot.tz  = 0.0;
    
    trans = TranslationMatrix (-camPos);
    
    return (rot * trans);
}


// --------------------------------------------------------------------------
// FrustumMatrix
//
// Setup a frustum matrix given the left, right, bottom, top, near, and far
// values for the frustum boundaries.
// --------------------------------------------------------------------------

template <typename Real>
inline Matrix4x4<Real>
FrustumMatrix (Real l, Real r, Real b, Real t, Real n, Real f)
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
// PerspectiveMatrix
//
// Setup a perspective matrix given the field-of-view in the Y direction
// in degrees, the aspect ratio of Y/X, and near and far plane distances.
// --------------------------------------------------------------------------

template <typename Real>
inline Matrix4x4<Real>
PerspectiveMatrix (Real fovY, Real aspect, Real n, Real f)
{
    Matrix4x4<Real> res;
    
    Real angle;
    Real cot;
    
    angle = fovY / 2.0;
    angle = degToRad (angle);
    
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
// OrthoMatrix
//
// Setup an orthographic Matrix4x4 given the left, right, bottom, top, near,
// and far values for the frustum boundaries.
// --------------------------------------------------------------------------

template <typename Real>
inline Matrix4x4<Real>
OrthoMatrix (Real l, Real r, Real b, Real t, Real n, Real f)
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
// OrthoNormalMatrix
//
// Setup an orientation matrix using 3 basis normalized vectors.
// --------------------------------------------------------------------------

template <typename Real>
inline Matrix4x4<Real>
OrthoNormalMatrix (const Vector3<Real> &xdir, const Vector3<Real> &ydir,
                   const Vector3<Real> &zdir)
{
    Matrix4x4<Real> res;
    
    res.m[0] = xdir.x; res.m[4] = ydir.x; res.m[8] = zdir.x; res.m[12] = 0.0;
    res.m[1] = xdir.y; res.m[5] = ydir.y; res.m[9] = zdir.y; res.m[13] = 0.0;
    res.m[2] = xdir.z; res.m[6] = ydir.z; res.m[10]= zdir.z; res.m[14] = 0.0;
    res.m[3] = 0.0;    res.m[7] = 0.0;    res.m[11]= 0.0;    res.m[15] = 1.0;
    
    return res;
}
