#include <vector>

namespace lflr {

class CSRMatrix 
{
public:
   CSRMatrix( int rows, int cols, int nonz );
   CSRMatrix( );
   CSRMatrix( int rows, int cols, int nonz, int *row_offs, int *col_inds, double *vals );
   void setRowSize( int );
   void setColSize( int );
   void setNumNonz( int );
   void setValMTX( int, int , double, int , bool);
   void spmv( std::vector<double> &vin, std::vector<double> &vout);
   int getRowSize();
   int getColSize();
   int getNumNonz();
   const std::vector<int> & get_col_index() const { return col_index;};
   const std::vector<int> & get_row_offset() const { return row_offset;};
   const std::vector<double> & get_values() const { return values;};
private:
   int nRows, nCols;
   int nNonz;
   std::vector<int>  row_offset; 
   std::vector<int>  col_index;
   std::vector<double>  values;

};

CSRMatrix::CSRMatrix () 
{
  nRows = 1;
  nCols = 1;
  nNonz = 1;
  row_offset.resize(2);
  col_index.resize(1);
  values.resize(1);

}

CSRMatrix::CSRMatrix( int rows, int cols, int nonz ) 
{
  nRows = rows;
  nCols = cols;
  nNonz = nonz;
  row_offset.resize(nRows+1);
  col_index.resize(nonz );
  values.resize(nonz );
  row_offset[nRows] = nonz ;
  row_offset[0] = 0;
}

CSRMatrix::CSRMatrix( int rows, int cols, int nonz, int *row_offs, int *col_inds, double *vals )
{
  nRows = rows;
  nCols = cols;
  nNonz = nonz;
  row_offset.assign( row_offs, row_offs+(rows+1) );
  col_index.assign( col_inds, col_inds+nonz );
  values.assign( vals, vals+nonz );
}

void CSRMatrix::setRowSize( int size )
{
  nRows = size;
  row_offset.resize(nRows+1);
  row_offset[0] = 0;
}

void CSRMatrix::setColSize( int size )
{
  nCols = size;
}

void CSRMatrix::setNumNonz( int size )
{
  nNonz = size;
  col_index.resize(size);
  values.resize(size);
  row_offset[nRows] = size;
  row_offset[0] = 0;
}

void CSRMatrix::setValMTX( int col, int row, double val, int loc, bool isnewrow )
{
  col_index[loc] = col;
  values[loc] = val;
  if( isnewrow )
  {
    row_offset[row] = loc;
  }
}

int CSRMatrix::getNumNonz()
{
  return nNonz;
}

int CSRMatrix::getRowSize()
{
  return nRows;
}

int CSRMatrix::getColSize()
{
  return nCols;
}

void CSRMatrix::spmv( std::vector<double> &vin, std::vector<double> &vout)
{
  for( int i = 0; i < nRows; ++i )
  {
    double tmp_dot = 0;
    for( int j = row_offset[i]; j < row_offset[i+1]; ++j ) 
    {
      tmp_dot += values[j]*vin[col_index[j]];
    }
    vout[i] = tmp_dot;
  }
}
}
