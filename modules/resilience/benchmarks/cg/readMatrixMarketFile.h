#include <vector>
#include <string>
#include <iostream>
#include <fstream>  
#include "CSRMatrix.h"
#include "formfull_CSR.h"


std::vector<std::string> split(const char *str, char c = ' ')
{
    std::vector<std::string> result;

    do
    {
        const char *begin = str;

        while(*str != c && *str)
            str++;

        result.push_back(std::string(begin, str));
    } while (0 != *str++);

    return result;
}


lflr::CSRMatrix readMatrixMarketFile( const std::string &sparseMatrixFile )
{
  // Read Sparse Matrix
  std::string line;
  std::vector<std::string> sizes;
  std::ifstream myfile;
  myfile.open (sparseMatrixFile);
  bool gotSize = false;
  int nRows, nCols, nNonz, iNonz = 0, prevRow  = 0;
  std::string::size_type sz;
  //Dynamic memory to store dat from half matrix, which will be extended
  //to the full matrix later on
  int *colind, *rowptr;
  double *nzval;

  if( myfile.is_open() )
  {
    while( std::getline ( myfile, line ) )
    {
      if( line[0] != '%' )  // Process the line
      {
        if( gotSize == false )
        {
          sizes = split( line.c_str(), ' ' );
          nRows = std::stoi( sizes[0] );
          nCols = std::stoi( sizes[1] );
          nNonz = std::stoi( sizes[2] );

	  //malloc buffers. Remember, these are only for the half matrix yet.
	  nzval = (double*)malloc(sizeof(double)*nNonz);
	  colind = (int*)malloc(sizeof(int)*nNonz);
	  rowptr = (int*)malloc(sizeof(int)*(nRows+1));
	  

	  gotSize = true;
        }
        else // Insert sparse elements
        {
          sizes = split( line.c_str(), ' ' );
          int col = std::stoi(sizes[0]);
          int row = std::stoi(sizes[1]);
          double val = std::stod(sizes[2],&sz);

          if( prevRow == row-1 )
          {
	    colind[iNonz] = col-1;
	    nzval[iNonz] = val;
          }
          else
          {
	    colind[iNonz] = col-1;
	    nzval[iNonz] = val;
	    rowptr[row-1] = iNonz;
            prevRow = row-1;
          }
          iNonz ++ ;
        }
        rowptr[nRows] = iNonz;
      }
    }
    myfile.close();

    //Now that the file has been read into the upper half of the matrix
    //form the full matrix
    FormFullA(nRows, &nNonz, &nzval, &colind, &rowptr);

    lflr::CSRMatrix spmat( nRows, nRows, nNonz, rowptr, colind, nzval);
    free(nzval); free(colind); free(rowptr);
    return spmat;
  }
  else
  {
    std::cout << "Error: File " << sparseMatrixFile << " is not found!" << std::endl;
    lflr::CSRMatrix spmat;
    return spmat;
  }
}

bool checkColIndices( lflr::CSRMatrix & spmat)
{
  bool sane = true;

  std::vector<int> row_offset = spmat.get_row_offset();
  std::vector<int> col_index = spmat.get_col_index();

  //Loop over rows, and for each row, get the set of column indices, check that they are in ascending order
  for (int r = 0; r < spmat.getRowSize(); r++)
    {
      //The column indices for this row go from col_index[ row_offset[r] ]...col_index[ row_offset[r+1] - 1]
      for(int col_idx_itr = row_offset[r]; col_idx_itr < row_offset[r+1] - 1; col_idx_itr++)
	{
	  if(col_index[col_idx_itr] >= col_index[col_idx_itr+1]) {
	    sane = false;
	    std::cout << "Check failed on row " << r <<" "<< col_idx_itr <<" "<<col_index[col_idx_itr]<<" "<<col_index[col_idx_itr+1]<<std::endl;
	  }
	}
    }

  return sane;
}

