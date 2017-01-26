//
//  matrixstacks.cpp
//  BitThemAll
//
//  Created by Jorge López González on 05/11/12.
//  Copyright (c) 2012 Jorge López González. All rights reserved.
//

#include "matrixstacks.h"

#include "matrix.h"

#include <cassert>

namespace dc
{
namespace math
{
    void MatrixStacks::Init()
    {
        // Hay tantas pilas como unidades de textura y las pilas de modelo/vista y proyeccion
        ml_matrixStacks = TMatrixStackList(mn_textureUnits + 2);
        
        // Inicializamos todas las pilas con una matriz vacia
        TMatrixStackList::iterator it, end;
        it = ml_matrixStacks.begin();
        end = ml_matrixStacks.end();
        for (; it!=end; ++it)
        {
            (*it).push(Matrix4x4f());
        }
        
        me_mode = MODELVIEW;
    }
    
    void MatrixStacks::End()
    {
        ml_matrixStacks.clear();
    }
    
    void MatrixStacks::SetMatrixMode(const EMatrixMode mode)
    {
        me_mode = MODELVIEW;
    }
    
    Matrix4x4f& MatrixStacks::ModelViewMatrix()
    {
        return ml_matrixStacks[MODELVIEW].back();
    }
    
    Matrix4x4f& MatrixStacks::ProjectionMatrix()
    {
        return ml_matrixStacks[PROJECTION].back();
    }
    
    Matrix4x4f& MatrixStacks::TextureMatrix(const int index)
    {
        assert(0<=index && index<=mn_textureUnits);
        return ml_matrixStacks[index + 2].back();
    }
    
    void MatrixStacks::Push()
    {
        ml_matrixStacks[me_mode].push(Matrix4x4f());
    }
    
    void MatrixStacks::Pop()
    {
        ml_matrixStacks[me_mode].pop();
    }
    
    void MatrixStacks::LoadIdentity()
    {
        ml_matrixStacks[me_mode].back().Identity();
    }
    
    void MatrixStacks::LoadMatrix(const Matrix4x4f& m)
    {
        ml_matrixStacks[me_mode].back() = m;
    }
    
    void MatrixStacks::Translate(const Vector3f& vec)
    {
        ml_matrixStacks[me_mode].back().Translate(vec);
    }
    
    void MatrixStacks::Rotate(const Vector3f& vec)
    {
        ml_matrixStacks[me_mode].back().TransformRotation(vec);
    }
    
    void MatrixStacks::Scale(const Vector3f& vec)
    {
        ml_matrixStacks[me_mode].back().Scale(vec);
    }
    
    void MatrixStacks::Multiply(const Matrix4x4f& m)
    {
        ml_matrixStacks[me_mode].back() *= m;
    }
    
    void MatrixStacks::Perspective( float fovy, float aspect,float near, float far )
    {
        ml_matrixStacks[me_mode].back().PerspectiveMatrix(fovy, aspect, near, far);
    }
    
    void MatrixStacks::Ortho(float left, float right, float bottom, float top, float near, float far)
    {
        ml_matrixStacks[me_mode].back().OrthoMatrix(left, right, bottom, top, near, far);
    }
}
}