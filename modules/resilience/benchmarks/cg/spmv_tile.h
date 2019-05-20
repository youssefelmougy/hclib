#include <vector>
#include <cassert>
#include <string>
#include <iterator>
#include <map>
#include <set>
#include <iostream>
#include <algorithm>
#include "readMatrixMarketFile.h"

/*
  Container to hold the vector tiles, that will be passed along 
  from one iteration to next through promises and futures.
  Inherits from std::vector<double>
 */
class DoubleVec : public hclib::resilience::obj,
                 public std::vector<double> {
public:
  virtual ~DoubleVec() {}

  bool equals(obj* obj2){
    DoubleVec *o2 = (DoubleVec*)obj2;
    return *this == *o2;
  }


};

class spmv_tile
{
  /* Class to create tiles of a CSR sparse matrix and dense vector
     to perform decomposed sparse matrix vector multiply. 
     Setup for structural mechanics problems, so 

     ******ASSUMES SQUARE MATRICES**************

   */

public:
  spmv_tile( );

  spmv_tile( lflr::CSRMatrix * matrix, int ntiles, int tile_id);

  void init_vector_in_tile( std::vector<double>& );

  void setup_vec_tiles_dependency();

  const int getStartRow() const { return start_row_; };
  const int getEndRow() const { return end_row_; };
  
  int start_row_for_tile(const int & tile);
  int getWhichTile(const int & index);

  const std::vector<int> & get_dependent_tiles_set() const {return vec_tiles_needed;};
  /*
  std::vector<double> & get_vector_of_tile() {return result_vec;};
  */
  void do_sanity_check(const std::map<int, DoubleVec* > & ) const;
  void do_tile_matvec(const std::map<int, DoubleVec* > &, DoubleVec * ) const;

private:

  //Basic ownership information
  lflr::CSRMatrix * matrix_ ; //pointer to the full sparse matrix
  std::vector<int> col_index_;
  std::vector<int> row_offset_;
  std::vector<double> values_; 
  int ntiles_ ;  //Total number of tiles of the sparse matrix 
  int tile_idx_; //which tile am I
  int start_row_;//what global row index do I start at
  int end_row_;  //what global row index do I end at
  /*
  std::vector<double> result_vec; //result (subset of full vec) I will publish at the end of each iteration
  */
  //Dependency information
  std::vector<int> vec_tiles_needed; //tile indices of vectors I need INCLUDING myself

  //For a specific row, maintain a list (vec<int>) of subset of vector elements I need
  //based on the column index of my non zero elements. This list is separate for each
  //vector tile I depend on, INCLUDING myself.

  //Now create such a map for each row I'm in charge of, store in a vector.

  //The result is a vector of maps

  //The elements of the outer vector span start_row_ ---> end_row_
  //Each element of this outer vector is an inner map.
  //The key for the inner map spans the vec_tiles_needed for this row.
  //The value for the inner map is the list of the subset indices needed from this vector for this row

  //Phew
  
  std::vector<std::map<int, std::vector<int> > > mapsOfVecTilesubsets;  

};

spmv_tile::spmv_tile( lflr::CSRMatrix * matrix, int ntiles, int tile_id) :
    matrix_(matrix), ntiles_(ntiles), tile_idx_(tile_id)
{
  col_index_ = matrix->get_col_index();
  row_offset_ = matrix->get_row_offset();
  values_ = matrix->get_values();

  if(tile_id < matrix->getRowSize()%ntiles ) { //this tile has +1 elems
    start_row_ = start_row_for_tile(tile_id);//tile_id * ( matrix->getRowSize()/ntiles + 1 );
    end_row_ = start_row_ + matrix->getRowSize()/ntiles;
  } else {
    start_row_ = start_row_for_tile(tile_id);//tile_id * ( matrix->getRowSize()/ntiles ) + matrix->getRowSize()%ntiles;
    end_row_ = start_row_ + matrix->getRowSize()/ntiles - 1;
  }
  
}

int spmv_tile::start_row_for_tile(const int & tile)
{
  if(tile < matrix_->getRowSize()%ntiles_ ) { //this tile has +1 elems
    return tile * ( matrix_->getRowSize()/ntiles_ + 1 );
  } else {
    return tile * ( matrix_->getRowSize()/ntiles_ ) + matrix_->getRowSize()%ntiles_;
  }
}

//Initialize the dense vector, for this tile
void spmv_tile::init_vector_in_tile(std::vector<double>& prev_vec )
{
  for(int i = start_row_; i<= end_row_ ; i++){
    /*
    result_vec.push_back(1.0);//((double)i);
    */
    prev_vec.push_back(1.0);//((double)i);
  }

}

//By processing the column indices of non zero elements
//of rows I own, identify the vector tiles I depend on
void spmv_tile::setup_vec_tiles_dependency()
{
  int at_tile, which_tile;

  //loop over rows I own
  for(int r = start_row_; r <= end_row_; r++) {
    std::map<int, std::vector<int> > mapOfSubindices; 

    at_tile = getWhichTile( col_index_[ row_offset_[r] ]  ); //First tile we depend on
    std::vector<int> subsets;

    int cnt = 0;
    //col indices for this row go from col_index[ row_offset[i] ] -> col_index[ row_offset[i+1] - 1 ]
    for(int j = row_offset_[r]  ; j < row_offset_[r+1]; j++) {
      int col = col_index_[j];
      
      which_tile = getWhichTile(col);

      if(which_tile == at_tile) {
        subsets.push_back( col - start_row_for_tile(which_tile) ); 
      } else if (which_tile > at_tile) {
        //We've gone past the current tile we're checking:
	// 1) Add this tile to the set of tiles this tile depends on (set.insert will insert only if not previously present)
	// 2) Commit the subsets list and move on
	//vec_tiles_needed.insert(at_tile); //1
	if(std::find( vec_tiles_needed.begin(),vec_tiles_needed.end(), at_tile) == vec_tiles_needed.end() ) {
	  vec_tiles_needed.push_back(at_tile); //1 Check and insert only if not inserted already
	}
        mapOfSubindices.insert( std::make_pair(at_tile, subsets) ); cnt += subsets.size();//2
        subsets.clear();
        at_tile = which_tile;
        subsets.push_back( col - start_row_for_tile(which_tile) ); 
      } else {
        std::cout << "Something went wrong on "<< tile_idx_ <<" row " << r << " column index " << col << std::endl;
      }
    }//end j for loop
    //Don't forget the last tile, which we would NOT have jumped past
    //vec_tiles_needed.insert(at_tile); //1
    if(std::find( vec_tiles_needed.begin(),vec_tiles_needed.end(), at_tile) == vec_tiles_needed.end() ) {
      vec_tiles_needed.push_back(at_tile); //1 Check and insert only if not inserted already
    }    
    mapOfSubindices.insert( std::make_pair(at_tile, subsets) ); //2
    cnt += subsets.size();     subsets.clear();

    assert( cnt == row_offset_[r+1] - row_offset_[r] );
    //if(tile_idx_ == 0) { std::cout << "Row: " << r << " " << cnt - row_offset_[r+1] + row_offset_[r]  << std::endl; }

    //Now that we've done processing this row, commit the map just accumulated to outer map
    mapsOfVecTilesubsets.push_back( mapOfSubindices );
  }
  //std::cout << "Num tiles I depend on "<< vec_tiles_needed.size() << " row maps are " << mapsOfVecTilesubsets.size() << std::endl;  
}

int spmv_tile::getWhichTile(const int & index)
{
  int M = matrix_->getRowSize();
  int K = ntiles_; 
  if( M%K == 0 ){ //equipartitioned tiles
    return index/( M/K );
  }  else {
    return ( index < (M%K)*(M/K + 1) ? index/( M/K + 1) : (index - M%K)/(M/K) );
  }
  
}

void spmv_tile::do_sanity_check(const std::map<int, DoubleVec* > & vec_tiles_I_seek) const
{

  assert(vec_tiles_I_seek.size() == vec_tiles_needed.size() );

  //Loop over rows I own
  //For each row:
  //  1) Get the inner map <int, std::vector<int>, key is vector_tile, value is subset of indices
  //  2) Loop over this map keys i.e. vector tiles depended on
  //  3) For each key, get the pointer to the vector corresponding to that tile
  //  4) Pluck the subset of that vector of that tile, pack into a tempry buffer
  //  5) Compare that with the column indices of this row (the two should be int==double equal)

  //loop over rows I own
  for(int r = start_row_; r <= end_row_; r++) {  
    std::vector<double> multiplicand; //buffer to be packed for multiplication of this row, empty to begin with
    std::map<int, std::vector<int> > inner_map = mapsOfVecTilesubsets[r-start_row_]; //1

    //iterate over this map. Keys are ids of vector tiles dependend on, by this row
    for( std::map<int, std::vector<int> >::iterator it = inner_map.begin(); it != inner_map.end(); ++it) { //2
      //it->first is the key i.e. the tile index of vector I depend on. Check if its present in tiles_I_seek
      int which_tile = it->first;
      //assert( vec_tiles_I_seek.find(which_tile) != vec_tiles_I_seek.end() );  
      std::vector<int> subindices = it->second; //subset of indices needed from which_tile 
      std::vector<double> seeking_vector = *( vec_tiles_I_seek.find(which_tile)->second );//3
      for(std::vector<int>::iterator tile_vec_itr = subindices.begin(); tile_vec_itr != subindices.end(); ++tile_vec_itr) {
	multiplicand.push_back(seeking_vector[*tile_vec_itr] ); //4, 5
      }
      
    }//end loop inner map

    //check that the size of multiplicand buffer equals nnZs for this row
    assert( multiplicand.size() == (row_offset_[r+1] - row_offset_[r]) );

    //Now check each element. The value in the multiplicand buffer must equal each column index
    std::vector<double>::iterator mult_itr = multiplicand.begin();
    for(int j = row_offset_[r]  ; j < row_offset_[r+1]; j++) {
      assert( (double)(col_index_[j]) == *mult_itr );
      mult_itr++;
    }

  }//end loop of r
  
}


void spmv_tile::do_tile_matvec(const std::map<int, DoubleVec* > & vec_tiles_I_seek, DoubleVec* result_vec) const
{
  assert(vec_tiles_I_seek.size() == vec_tiles_needed.size() );

  //Loop over rows I own
  //For each row:
  //  1) Get the inner map <int, std::vector<int>, key is vector_tile, value is subset of indices
  //  2) Loop over this map keys i.e. vector tiles depended on
  //  3) For each key, get the pointer to the vector corresponding to that tile
  //  4) Pluck the subset of that vector of that tile, multiply with corresponding column element of this row
  //  5) Update the running sum for this row, insert the value in result_vec at the end

  //loop over rows I own
  for(int r = start_row_; r <= end_row_; r++) {  
    int val_ptr = row_offset_[r]; //pointer that keeps track of nnZ elems of this row
    double res = 0.0; 

    std::map<int, std::vector<int> > inner_map = mapsOfVecTilesubsets[r-start_row_]; //1
    //iterate over this map. Keys are ids of vector tiles dependend on, by this row
    for( std::map<int, std::vector<int> >::iterator it = inner_map.begin(); it != inner_map.end(); ++it) { //2
      //it->first is the key i.e. the tile index of vector I depend on. Check if its present in tiles_I_seek
      int which_tile = it->first;
      //assert( vec_tiles_I_seek.find(which_tile) != vec_tiles_I_seek.end() );      
      std::vector<int> subindices = it->second; //subset of indices needed from which_tile
      std::vector<double> seeking_vector = *( vec_tiles_I_seek.find(which_tile)->second );//3
      for(std::vector<int>::iterator tile_vec_itr = subindices.begin(); tile_vec_itr != subindices.end(); ++tile_vec_itr) {
	res += values_[val_ptr] * seeking_vector[*tile_vec_itr]; //Actual multiplication of row elem with col elem from right vecor tile
	val_ptr++; 
      }
    }//end loop inner map
    
    //assert( val_ptr == row_offset_[r+1] ); 
    result_vec->push_back(res);
    /*
    result_vec[r-start_row_] = res;   
    */
  }//end loop of r
  
}

