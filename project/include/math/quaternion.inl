//
//  quaternion.inl.h
//  bitthemall
//
//  Created by Jorge López González on 17/07/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//


/////////////////////////////////////////////////////////////////////////////
//
// class Quaternion<Real> implementation.
//
/////////////////////////////////////////////////////////////////////////////

// --------------------------------------------------------------------------
// Quaternion::identity
//
// Set to identity
// --------------------------------------------------------------------------

template <typename Real>
inline void
Quaternion<Real>::Identity ()
{
    w = 1.0; 
    x = y = z = 0.0;
}


// --------------------------------------------------------------------------
// Quaternion::normalize
//
// "Normalize" a quaternion.  Note that normally, quaternions
// are always normalized (within limits of numerical precision).
//
// This function is provided primarily to combat floating point "error
// creep", which can occur when many successive quaternion operations
// are applied.
// --------------------------------------------------------------------------

template <typename Real>
inline void
Quaternion<Real>::Normalize ()
{
    // Compute magnitude of the quaternion
    Real mag = std::sqrt ((w * w) + (x * x) + (y * y) + (z * z));
    
    // Check for bogus length, to protect against divide by zero
    if (mag > 0.0)
    {
        // Normalize it
        Real oneOverMag = 1.0 / mag;
        
        w *= oneOverMag;
        x *= oneOverMag;
        y *= oneOverMag;
        z *= oneOverMag;
    }
}


// --------------------------------------------------------------------------
// Quaternion<Real>::computeW
//
// Compute the W component of a unit quaternion given its x,y,z components.
// --------------------------------------------------------------------------

template <typename Real>
inline void
Quaternion<Real>::ComputeW ()
{
    Real t = 1.0 - (x * x) - (y * y) - (z * z);
    
    if (t < 0.0)
	    w = 0.0;
    else
	    w = -std::sqrt (t);
}


// --------------------------------------------------------------------------
// Quaternion<Real>::rotate
//
// Rotate a point by quaternion.  v' = q.p.q*, where p = <0, v>.
// --------------------------------------------------------------------------

template <typename Real>
inline void
Quaternion<Real>::Rotate (Vector3<Real>& v) const
{
    Quaternion<Real> qf = *this * v * ~(*this);
    v.x = qf.x;
    v.y = qf.y;
    v.z = qf.z;
}


// --------------------------------------------------------------------------
// Quaternion::fromMatrix
//
// Setup the quaternion to perform a rotation, given the angular displacement
// in matrix form.
// --------------------------------------------------------------------------

template <typename Real>
inline void
Quaternion<Real>::FromMatrix (const Matrix4x4<Real>& m)
{
    Real trace = m.m11 + m.m22 + m.m33 + 1.0;
    
    if (trace > 0.0001)
    {
        Real s = 0.5 / std::sqrt (trace);
        w = 0.25 / s;
        x = (m.m23 - m.m32) * s;
        y = (m.m31 - m.m13) * s;
        z = (m.m12 - m.m21) * s;
    }
    else
    {
        if ((m.m11 > m.m22) && (m.m11 > m.m33))
		{
            Real s = 0.5 / std::sqrt (1.0 + m.m11 - m.m22 - m.m33);
            x = 0.25 / s;
            y = (m.m21 + m.m12) * s;
            z = (m.m31 + m.m13) * s;
            w = (m.m32 - m.m23) * s;
		}
        else if (m.m22 > m.m33)
		{
            Real s = 0.5 / std::sqrt (1.0 + m.m22 - m.m11 - m.m33);
            x = (m.m21 + m.m12) * s;
            y = 0.25 / s;
            z = (m.m32 + m.m23) * s;
            w = (m.m31 - m.m13) * s;
		}
        else
		{
            Real s = 0.5 / std::sqrt (1.0 + m.m33 - m.m11 - m.m22);
            x = (m.m31 + m.m13) * s;
            y = (m.m32 + m.m23) * s;
            z = 0.25 / s;
            w = (m.m21 - m.m12) * s;
		}
    }
}


// --------------------------------------------------------------------------
// Quaternion::fromEulerAngles
//
// Setup the quaternion to perform an object->inertial rotation, given the
// orientation in XYZ-Euler angles format.  x,y,z parameters must be in
// radians.
// --------------------------------------------------------------------------

template <typename Real>
static inline
Quaternion<Real>
Quaternion<Real>::FromEulerAngles (const Real roll, const Real pitch, const Real yaw) const
{
    // Compute sine and cosine of the half angles
    Real sr = std::sin (roll * 0.5);
    Real cr = std::cos (roll * 0.5);
    Real sp = std::sin (pitch * 0.5);
    Real cp = std::cos (pitch * 0.5);
    Real sy = std::sin (yaw * 0.5);
    Real cy = std::cos (yaw * 0.5);
    
    // Compute values
    return Quaternion<Real>( (cy * cp * cr) + (sy * sp * sr),
                            -(sy * sp * cr) + (cy * cp * sr),
                             (cy * sp * cr) + (sy * cp * sr),
                            -(cy * sp * sr) + (sy * cp * cr));
}


// --------------------------------------------------------------------------
// Quaternion::toEulerAngles
//
// Setup the Euler angles, given an object->inertial rotation quaternion.
// Returned x,y,z are in radians.
// --------------------------------------------------------------------------

template <typename Real>
static inline 
Vector3<Real>
Quaternion<Real>::ToEulerAngles (const Quaternion<Real>& q) const
{
    // Compute square products for optimization
    float sqw = q.w * q.w;
    float sqx = q.x * q.x;
    float sqy = q.y * q.y;
    float sqz = q.z * q.z;
    
    // Compute Y-axis angle
    y = std::asin (2.0 * ((x * z) + (w * y)));
    
    // Compute cos and one over cos for optimization
    Real cy = std::cos (y);
    Real oneOverCosY = 1.0 / cy;
    
    if (std::fabs (cy) > 0.001)
    {
        // No gimball lock
        x = std::atan2 (2.0 * ((w * x) - (y * z)) * oneOverCosY, (1.0 - 2.0 * (sqx + sqy)) * oneOverCosY);
        z = std::atan2 (2.0 * ((w * z) - (x * y)) * oneOverCosY, (1.0 - 2.0 * (sqy + sqz)) * oneOverCosY);
    }
    else
    {
        // Gimbal lock case
        x = 0.0;
        z = std::atan2 (2.0 * ((x * y) + (w * z)),
                        1.0 - 2.0 * (sqx + sqz));
    }
    
    double sqw;
    double sqx;
    double sqy;
    double sqz;
    
    double rotxrad;
    double rotyrad;
    double rotzrad;
    
    sqw = _w * _w;
    sqx = _x * _x;
    sqy = _y * _y;
    sqz = _z * _z;
    
    rotxrad = (double)atan2l(2.0 * ( _y * _z + _x * _w ) , ( -sqx - sqy + sqz + sqw ));
    
    rotyrad = (double)asinl(-2.0 * ( _x * _z - _y * _w ));
    
    rotzrad = (double)atan2l(2.0 * ( _x * _y + _z * _w ) , (  sqx - sqy - sqz + sqw ));
    
    float pi = 3.141592;
    return maths::CVector3((rotxrad*180)/pi,(rotyrad*180)/pi,(rotzrad*180)/pi);
}


// --------------------------------------------------------------------------
// Quaternion::getRotationAngle
//
// Return the rotation angle theta (in radians).
// --------------------------------------------------------------------------

template <typename Real>
inline Real
Quaternion<Real>::rotationAngle () const
{
    // Compute the half angle.  Remember that w = cos(theta / 2)
    Real thetaOver2 = safeAcos (w);
    
    // Return the rotation angle
    return thetaOver2 * 2.0;
}


// --------------------------------------------------------------------------
// Quaternion::getRotationAxis
//
// Return the rotation axis.
// --------------------------------------------------------------------------

template <typename Real>
inline Vector3<Real>
Quaternion<Real>::rotationAxis () const
{
    // Compute sin^2(theta/2).  Remember that w = cos(theta/2),
    // and sin^2(x) + cos^2(x) = 1
    Real sinThetaOver2Sq = 1.0 - (w * w);
    
    // Protect against numerical imprecision
    if (sinThetaOver2Sq <= 0.0)
    {
        // Identity quaternion, or numerical imprecision.  Just
        // return any valid vector, since it doesn't matter
        
        return Vector3<Real> (1.0, 0.0, 0.0);
    }
    
    // Compute 1 / sin(theta/2)
    Real oneOverSinThetaOver2 = 1.0 / std::sqrt (sinThetaOver2Sq);
    
    // Return axis of rotation
    return Vector3<Real> (
                          x * oneOverSinThetaOver2,
                          y * oneOverSinThetaOver2,
                          z * oneOverSinThetaOver2
                          );
}


// --------------------------------------------------------------------------
// Quaternion operators
//
// Operator overloading for basic quaternion operations.
// --------------------------------------------------------------------------

template <typename Real>
inline Quaternion<Real>
Quaternion<Real>::operator+ (const Quaternion<Real> &q) const
{
    return Quaternion<Real> (w + q.w, x + q.x, y + q.y, z + q.z);
}

template <typename Real>
inline Quaternion<Real> &
Quaternion<Real>::operator+= (const Quaternion<Real> &q)
{
    w += q.w; x += q.x; y += q.y; z += q.z;
    return *this;
}

template <typename Real>
inline Quaternion<Real>
Quaternion<Real>::operator- (const Quaternion<Real> &q) const
{
    return Quaternion<Real> (w - q.w, x - q.x, y - q.y, z - q.z);
}

template <typename Real>
inline Quaternion<Real> &
Quaternion<Real>::operator-= (const Quaternion<Real> &q)
{
    w -= q.w; x -= q.x; y -= q.y; z -= q.z;
    return *this;
}

// Quaternion multiplication.  The order of multiplication, from left
// to right, corresponds to the order of concatenated rotations.
// NOTE: Quaternion multiplication is NOT commutative, p * q != q * p
template <typename Real>
inline Quaternion<Real>
Quaternion<Real>::operator* (const Quaternion<Real> &q) const
{
    // We use the Grassman product formula:
    // pq = <w1w2 - dot(v1, v2), w1v2 + w2v1 + cross(v1, v2)>
    return Quaternion<Real> (
                             (w * q.w) - (x * q.x) - (y * q.y) - (z * q.z),
                             (x * q.w) + (w * q.x) + (y * q.z) - (z * q.y),
                             (y * q.w) + (w * q.y) + (z * q.x) - (x * q.z),
                             (z * q.w) + (w * q.z) + (x * q.y) - (y * q.x)
                             );
}

template <typename Real>
inline Quaternion<Real> &
Quaternion<Real>::operator*= (const Quaternion<Real> &q)
{
    *this = *this * q;
    return *this;
}

template <typename Real>
inline Quaternion<Real>
Quaternion<Real>::operator* (const Vector3<Real> &v) const
{
    // q * v = q * p where p = <0,v>
    // Thus, we can simplify the operations.
    return Quaternion<Real> (
                             - (x * v.x) - (y * v.y) - (z * v.z),
                             (w * v.x) + (y * v.z) - (z * v.y),
                             (w * v.y) + (z * v.x) - (x * v.z),
                             (w * v.z) + (x * v.y) - (y * v.x)
                             );
}

template <typename Real>
inline Quaternion<Real> &
Quaternion<Real>::operator*= (const Vector3<Real> &v)
{
    *this = *this * v;
    return *this;
}

template <typename Real>
inline Quaternion<Real>
Quaternion<Real>::operator* (Real k) const
{
    return Quaternion<Real> (w * k, x * k, y * k, z * k);
}

template <typename Real>
inline Quaternion<Real> &
Quaternion<Real>::operator*= (Real k)
{
    w *= k; x *= k; y *= k; z *= k;
    return *this;
}

template <typename Real>
inline Quaternion<Real>
Quaternion<Real>::operator/ (Real k) const
{
    Real oneOverK = 1.0 / k;
    return Quaternion<Real> (w * oneOverK, x * oneOverK, y * oneOverK, z * oneOverK);
}

template <typename Real>
inline Quaternion<Real> &
Quaternion<Real>::operator/= (Real k)
{
    Real oneOverK = 1.0 / k;
    w *= oneOverK; x *= oneOverK; y *= oneOverK; z *= oneOverK;
    return *this;
}

// Quaternion conjugate
template <typename Real>
inline Quaternion<Real>
Quaternion<Real>::operator~ () const
{
    return Quaternion<Real> (w, -x, -y, -z);
}


// Quaternion negation
template <typename Real>
inline Quaternion<Real>
Quaternion<Real>::operator- () const
{
    return Quaternion<Real> (-w, -x, -y, -z);
}


// --------------------------------------------------------------------------
//
// Nonmember Quaternion functions
//
// --------------------------------------------------------------------------

// Scalar on left multiplication
template <typename Real>
inline Quaternion<Real>
operator* (Real k, const Quaternion<Real> &q)
{
    return Quaternion<Real> (q.w * k, q.x * k, q.y * k, q.z * k);
}

// Quaternion dot product
template <typename Real>
inline Real
DotProduct (const Quaternion<Real> &a, const Quaternion<Real> &b)
{
    return ((a.w * b.w) + (a.x * b.x) + (a.y * b.y) + (a.z * b.z));
}

// Compute the quaternion conjugate.  This is the quaternian
// with the opposite rotation as the original quaternion.
template <typename Real>
inline Quaternion<Real>
Conjugate (const Quaternion<Real> &q)
{
    return Quaternion<Real> (q.w, -q.x, -q.y, -q.z);
}


// Compute the inverse quaternion (for unit quaternion only).
template <typename Real>
inline Quaternion<Real>
Inverse (const Quaternion<Real> &q)
{
    // Assume this is a unit quaternion! No check for this!
    Quaternion<Real> res (q.w, -q.x, -q.y, -q.z);
    res.normalize ();
    return res;
}


// --------------------------------------------------------------------------
// RotationQuaternion
//
// Setup the quaternion to rotate about the specified axis.  theta must
// be in radians.
// --------------------------------------------------------------------------

template <typename Real>
inline Quaternion<Real>
RotationQuaternion (Axis axis, Real theta)
{
    Quaternion<Real> res;
    
    // Compute the half angle
    Real thetaOver2 = theta * 0.5;
    
    // Set the values
    switch (axis)
    {
	    case kXaxis:
            res.w = std::cos (thetaOver2);
            res.x = std::sin (thetaOver2);
            res.y = 0.0;
            res.z = 0.0;
            break;
            
	    case kYaxis:
            res.w = std::cos (thetaOver2);
            res.x = 0.0;
            res.y = std::sin (thetaOver2);
            res.z = 0.0;
            break;
            
	    case kZaxis:
            res.w = std::cos (thetaOver2);
            res.x = 0.0;
            res.y = 0.0;
            res.z = std::sin (thetaOver2);
            break;
            
	    default:
            // Bad axis
            assert (false);
    }
    
    return res;
}

template <typename Real>
inline Quaternion<Real>
RotationQuaternion (const Vector3<Real> &axis, Real theta)
{
    Quaternion<Real> res;
    
    // The axis of rotation must be normalized
    assert (std::fabs (DotProduct (axis, axis) - 1.0) < 0.001);
    
    // Compute the half angle and its sin
    Real thetaOver2 = theta * 0.5;
    Real sinThetaOver2 = std::sin (thetaOver2);
    
    // Set the values
    res.w = std::cos (thetaOver2);
    res.x = axis.x * sinThetaOver2;
    res.y = axis.y * sinThetaOver2;
    res.z = axis.z * sinThetaOver2;
    
    return res;
}


// --------------------------------------------------------------------------
// Log
//
// Unit quaternion logarithm. log(q) = log(cos(theta) + n*sin(theta))
// --------------------------------------------------------------------------

template <typename Real>
inline Quaternion<Real>
Log (const Quaternion<Real> &q)
{
    Quaternion<Real> res;
    res.w = 0.0;
    
    if (std::fabs (q.w) < 1.0)
    {
        Real theta = std::acos (q.w);
        Real sintheta = std::sin (theta);
        
        if (std::fabs (sintheta) > 0.00001)
		{
            Real thetaOverSinTheta = theta / sintheta;
            res.x = q.x * thetaOverSinTheta;
            res.y = q.y * thetaOverSinTheta;
            res.z = q.z * thetaOverSinTheta;
            return res;
		}
    }
    
    res.x = q.x;
    res.y = q.y;
    res.z = q.z;
    
    return res;
}


// --------------------------------------------------------------------------
// Exp
//
// Quaternion exponential.
// --------------------------------------------------------------------------

template <typename Real>
inline Quaternion<Real>
Exp (const Quaternion<Real> &q)
{
    Real theta = std::sqrt (DotProduct (q, q));
    Real sintheta = std::sin (theta);
    
    Quaternion<Real> res;
    res.w = std::cos (theta);
    
    if (std::fabs (sintheta) > 0.00001)
    {
        Real sinThetaOverTheta = sintheta / theta;
        res.x = q.x * sinThetaOverTheta;
        res.y = q.y * sinThetaOverTheta;
        res.z = q.z * sinThetaOverTheta;
    }
    else
    {
        res.x = q.x;
        res.y = q.y;
        res.z = q.z;
    }
    
    return res;
}


// --------------------------------------------------------------------------
// Pow
//
// Quaternion exponentiation.
// --------------------------------------------------------------------------

template <typename Real>
inline Quaternion<Real>
Pow (const Quaternion<Real> &q, Real exponent)
{
    // Check for the case of an identity quaternion.
    // This will protect against divide by zero
    if (std::fabs (q.w) > 0.9999)
	    return q;
    
    // Extract the half angle alpha (alpha = theta/2)
    Real alpha = std::acos (q.w);
    
    // Compute new alpha value
    Real newAlpha = alpha * exponent;
    
    // Compute new quaternion
    Vector3<Real> n (q.x, q.y, q.z);
    n *= std::sin (newAlpha) / std::sin (alpha);
    
    return Quaternion<Real> (std::cos (newAlpha), n);
}


// --------------------------------------------------------------------------
// Slerp
//
// Spherical linear interpolation.
// --------------------------------------------------------------------------

template <typename Real>
inline Quaternion<Real>
Slerp (const Quaternion<Real> &q0, const Quaternion<Real> &q1, Real t)
{
    // Check for out-of range parameter and return edge points if so
    if (t <= 0.0) return q0;
    if (t >= 1.0) return q1;
    
    // Compute "cosine of angle between quaternions" using dot product
    Real cosOmega = DotProduct (q0, q1);
    
    // If negative dot, use -q1.  Two quaternions q and -q
    // represent the same rotation, but may produce
    // different slerp.  We chose q or -q to rotate using
    // the acute angle.
    Real q1w = q1.w;
    Real q1x = q1.x;
    Real q1y = q1.y;
    Real q1z = q1.z;
    
    if (cosOmega < 0.0)
    {
        q1w = -q1w;
        q1x = -q1x;
        q1y = -q1y;
        q1z = -q1z;
        cosOmega = -cosOmega;
    }
    
    // We should have two unit quaternions, so dot should be <= 1.0
    assert (cosOmega < 1.1);
    
    // Compute interpolation fraction, checking for quaternions
    // almost exactly the same
    Real k0, k1;
    
    if (cosOmega > 0.9999)
    {
        // Very close - just use linear interpolation,
        // which will protect againt a divide by zero
        
        k0 = 1.0 - t;
        k1 = t;
    }
    else
    {
        // Compute the sin of the angle using the
        // trig identity sin^2(omega) + cos^2(omega) = 1
        Real sinOmega = std::sqrt (1.0 - (cosOmega * cosOmega));
        
        // Compute the angle from its sin and cosine
        Real omega = std::atan2 (sinOmega, cosOmega);
        
        // Compute inverse of denominator, so we only have
        // to divide once
        Real oneOverSinOmega = 1.0 / sinOmega;
        
        // Compute interpolation parameters
        k0 = std::sin ((1.0 - t) * omega) * oneOverSinOmega;
        k1 = std::sin (t * omega) * oneOverSinOmega;
    }
    
    // Interpolate and return new quaternion
    return Quaternion<Real> (
                             (k0 * q0.w) + (k1 * q1w),
                             (k0 * q0.x) + (k1 * q1x),
                             (k0 * q0.y) + (k1 * q1y),
                             (k0 * q0.z) + (k1 * q1z)
                             );
}


// --------------------------------------------------------------------------
// Squad
//
// Spherical cubic interpolation.
// squad = slerp (slerp (q0, q1, t), slerp (qa, qb, t), 2t(1 - t)).
// --------------------------------------------------------------------------

template <typename Real>
inline Quaternion<Real>
Squad (const Quaternion<Real> &q0, const Quaternion<Real> &qa,
       const Quaternion<Real> &qb, const Quaternion<Real> &q1, Real t)
{
    Real slerpt = 2.0 * t * (1.0 - t);
    
    Quaternion<Real> slerpq0 = Slerp (q0, q1, t);
    Quaternion<Real> slerpq1 = Slerp (qa, qb, t);
    
    return Slerp (slerpq0, slerpq1, slerpt);
}


// --------------------------------------------------------------------------
// Intermediate
//
// Compute intermediate quaternions for building spline segments.
// --------------------------------------------------------------------------

template <typename Real>
inline void
Intermediate (const Quaternion<Real> &qprev, const Quaternion<Real> &qcurr,
		      const Quaternion<Real> &qnext, Quaternion<Real> &qa, Quaternion<Real> &qb)
{
    // We should have unit quaternions
    assert (DotProduct (qprev, qprev) <= 1.0001);
    assert (DotProduct (qcurr, qcurr) <= 1.0001);
    
    Quaternion<Real> invprev = Conjugate (qprev);
    Quaternion<Real> invcurr = Conjugate (qcurr);
    
    Quaternion<Real> p0 = invprev * qcurr;
    Quaternion<Real> p1 = invcurr * qnext;
    
    Quaternion<Real> arg = (Log (p0) - Log (p1)) * 0.25;
    
    qa = qcurr * Exp ( arg);
    qb = qcurr * Exp (-arg);
}

