//
//  matrixstacks.h
//  BitThemAll
//
//  Created by Jorge López González on 05/11/12.
//  Copyright (c) 2012 Jorge López González. All rights reserved.
//

#ifndef __BitThemAll__matrixstacks__
#define __BitThemAll__matrixstacks__

#include <queue>
#include <vector>

#include "vector.h"
#include "matrix.h"

namespace dc
{
namespace math
{
    
    enum EMatrixMode
    {
        TEXTURE0 = 0,
        TEXTURE1 = 1,
        TEXTURE2 = 2,
        TEXTURE3 = 3,
        TEXTURE4 = 4,
        TEXTURE5 = 5,
        TEXTURE6 = 6,
        TEXTURE7 = 7,
        MODELVIEW = 8,
        PROJECTION = 9,
    };
    
    typedef std::queue<Matrix4x4f>      TMatrixStack;
    typedef std::vector<TMatrixStack>   TMatrixStackList;
    
    class MatrixStacks
    {
    public:
        MatrixStacks(const unsigned int textureunits);
        ~MatrixStacks();
        
        void Init();
        void End();
        
        void SetMatrixMode(const EMatrixMode mode);
        
    public:
        Matrix4x4f& ModelViewMatrix();
        Matrix4x4f& ProjectionMatrix();
        Matrix4x4f& TextureMatrix(const int index);
        
    public:
        void Push();
        void Pop();
        
        void LoadIdentity();
        void LoadMatrix(const Matrix4x4f& m);
        
        void Translate(const Vector3f& vec );
        void Rotate(const Vector3f& vec );
        void Scale(const Vector3f& vec );
        void Multiply(const Matrix4x4f& m);
        
        void Perspective( float fovy, float aspect,float near, float far);
        void Ortho(float left, float right, float bottom, float top, float near, float far);
        
    private:
        unsigned int        mn_textureUnits;
        EMatrixMode         me_mode;
        
        TMatrixStackList    ml_matrixStacks;
    };
}
}

#endif /* defined(__BitThemAll__matrixstacks__) */
