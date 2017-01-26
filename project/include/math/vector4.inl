//
//  vector4.inl.h
//  bitthemall
//
//  Created by Jorge López González on 17/07/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

// Decide if a plane is facing a point
template <typename Real>
inline bool
IsFacingTo (const Vector4<Real> plane, const Vector3<Real> &v)
{
    return ((plane.x*v.x) + (plane.y*v.y) + (plane.z*v.z) + (plane.w))>0;
}

// Decide if a plane is facing a sphere
template <typename Real>
inline bool
IsFacingTo (const Vector4<Real> plane, const Vector3<Real> &v, Real radius)
{
    return ((plane.x*v.x) + (plane.y*v.y) + (plane.z*v.z) + (plane.w))>-radius;
}
