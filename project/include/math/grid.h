//
//  Grid.h
//  Grid
//
//  Created by Jorge López on 08/07/15.
//  Copyright (c) 2015 Jorge López. All rights reserved.
//

#ifndef __Grid__Grid__
#define __Grid__Grid__

#include "math.h"
#include "vector2.h"
#include "rect.h"

namespace dc
{
namespace math
{
    struct Cell
    {
        unsigned int column;
        unsigned int row;
        
        Cell (const unsigned int column, const unsigned int row): column(column), row (row){}
    };
    
    class Grid
    {
    // Constructors / Destructors
    public:
        Grid (const unsigned int columns, const unsigned int rows,
              const float x, const float y,
              const float width, const float height,
              const Vector2f& cellPivot):
            columns (columns),
            rows (rows),
            position(Vector2f(x, y)),
            dimensions(Vector2f(width, height)),
            cellDimensions(Vector2f(height / (float) rows, width / (float) columns)),
            cellPivot(cellPivot)
        {}
        
        Grid (const unsigned int columns, const unsigned int rows,
              const Vector2f& position,
              const Vector2f& dimensions,
              const Vector2f& cellPivot):
        columns (columns),
        rows (rows),
        position(position),
        dimensions(dimensions),
        cellDimensions(Vector2f(dimensions.height / (float) rows, dimensions.width / (float) columns)),
        cellPivot(cellPivot)
        {}
        
        virtual ~Grid() {}
        
    // Accessors
    public:
        const unsigned int TotalCells() const { return columns * rows; }

        const unsigned int Columns() const { return columns; }
        const unsigned int Rows() const { return rows; }
        
        const Vector2f& Position() const { return position; }
        const float X() const { return position.x; }
        const float Y() const { return position.y; }
        
        void SetPosition(const Vector2f& position)
        {
            this->position = position;
        }
        
        const float     Width() const { return dimensions.width; }
        const float     Height() const { return dimensions.height; }
        const Vector2f& Dimensions() const { return dimensions; }
        
        void SetDimensions(const Vector2f& dimensions)
        {
            this->dimensions = dimensions;
        }
        
        const float     CellWidth() const { return cellDimensions.width; }
        const float     CellHeight() const { return cellDimensions.height; }
        const Vector2f& CellDimensions() const { return cellDimensions; }
        
        void SetCellDimensions(const Vector2f& cellDimensions)
        {
            this->cellDimensions = cellDimensions;
        }
        
    private:
        const float GetXPos(const unsigned int column) const { return X() + (column * CellWidth()); }
        const float GetYPos(const unsigned int row) const { return Y() + (row * CellHeight()); }
    
    // Functions
    public:
        const Cell CellFromPosition(const Vector2f& position) const
        {
            unsigned int column = (unsigned int) (position.x / CellWidth());
            unsigned int row = (unsigned int) (position.y / CellHeight());
            return Cell (column, row);
        }
        
        const Vector2f CenterPos(const Cell& cell) const
        {
            return CenterPos(cell.row, cell.column);
        }
        
        const Vector2f CenterPos(const unsigned int row, const unsigned int column) const
        {
            return Vector2f(GetXPos(column) + Half(CellWidth()), GetYPos(row) + Half(CellHeight()));
        }
        
        const Vector2f UpperLeftCornerPos(const Cell& cell) const
        {
            return UpperLeftCornerPos(cell.row, cell.column);
        }
        
        const Vector2f UpperLeftCornerPos(const unsigned int row, const unsigned int column) const
        {
            return Vector2f(GetXPos(column), GetYPos(row));
        }
        
        const Vector2f LowerLeftCornerPos(const Cell& cell) const
        {
            return LowerLeftCornerPos(cell.row, cell.column);
        }
        
        const Vector2f LowerLeftCornerPos(const unsigned int row, const unsigned int column) const
        {
            return Vector2f(GetXPos(column), GetYPos(row) + CellHeight());
        }
        
        const Vector2f UpperRightCornerPos(const Cell& cell) const
        {
            return UpperRightCornerPos(cell.row, cell.column);
        }
        
        const Vector2f UpperRightCornerPos(const unsigned int row, const unsigned int column) const
        {
            return Vector2f(GetXPos(column) + CellWidth(), GetYPos(row));
        }
        
        const Vector2f LowerRightCornerPos(const Cell& cell) const
        {
            return LowerRightCornerPos(cell.row, cell.column);
        }
        
        const Vector2f LowerRightCornerPos(const unsigned int row, const unsigned int column) const
        {
            return Vector2f(GetXPos(column) + CellWidth(), GetYPos(row) + CellHeight());
        }
        
        const Rectf CellRect (const Cell& cell) const
        {
            return CellRect(cell.row, cell.column);
        }
        
        const Rectf CellRect (const unsigned int row, const unsigned int column) const
        {
            return Rectf(cellPivot, UpperLeftCornerPos(row, column), Vector2f(CellWidth(), CellHeight()));
        }

    // Fields
    private:
        
        unsigned int    columns;
        unsigned int    rows;
        
        Vector2f        position;
        Vector2f        dimensions;
        
        Vector2f        cellDimensions;
        Vector2f        cellPivot;
    };
}
}
#endif /* defined(__Grid__Grid__) */
