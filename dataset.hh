/* ----------------------------------------------------------------/
This is open source software; you can redistribute it and/or modify
it under the terms of the BSD-3 License. If software is modified
to produce derivative works, such modified software should be
clearly marked, so as not to confuse it with the version available
from LANL. Full text of the BSD-3 License can be found in the
License file of the repository.
/---------------------------------------------------------------- */

#ifndef _DATASET_HH
#define _DATASET_HH

#include <iostream>
#include <cmath>
#include <cassert>

//-----------------------------------------------------------------------------
// similar to full matrix data type
//-----------------------------------------------------------------------------

template<typename TYPE>
class TDataMat {

private:
  TYPE * mat ;
  int nr ;
  int nc ;
  bool status ; // true  ==> has already been set
                // false ==> must be set

  TDataMat( TDataMat const & ) {} // disable copy CTOR

public:
  TDataMat( int _nr=0, int _nc=0 ) : nr(_nr), nc(_nc), status(false) {
    if ( nr>0 && nc>0 ) { setup( _nr, _nc ) ; }
  }
  ~TDataMat() {
    if ( status ) { delete[] mat ; }
  }
  // initial setup
  void setup( int _nr, int _nc, TYPE t=TYPE(0) ) {
    assert( _nr>0 && _nc>0 ) ;
    if ( status ) { delete[] mat ; }
    status = true ;
    nr  = _nr ;
    nc  = _nc ;
    mat = new TYPE[nr*nc] ; 
    assert( mat != 0 ) ;
    for(int i=0;i<nr*nc; ++i) { mat[i] = t ; }
  }
  // access methods
  inline TYPE const operator() ( int i, int j ) const {
    assert( 0<=i && i<nr ) ;
    assert( 0<=j && j<nc ) ;
    return mat[j*nr+i] ;
  }
  inline TYPE & operator() ( int i, int j ) {
    assert( 0<=i && i<nr ) ;
    assert( 0<=j && j<nc ) ;
    return mat[j*nr+i] ;
  }
  // size methods
  inline int size() const { return nr  ; }
  inline int size_loc( int i ) const { 
    assert( 0<=i && i<nr ) ;
    return nc  ; 
  }
} ;

template<typename TYPE>
class TDataVec {
  
public:
  TYPE * vec ;
  int n ;
  bool status ;
  
  TDataVec( TDataVec const & ) {} // disable copy CTOR

public:
  TDataVec( int _n=0 ) : n(_n), status(false) {
    if ( n>0 ) { setup( _n ) ; }
  }
  ~TDataVec() {
    if ( status ) { delete[] vec ; }
  }
  // initial setup 
  void setup(int _n, TYPE t=0) {
    assert( !status ) ;
    status = bool(true) ;
    n = _n ;
    vec = new TYPE[n] ; 
    assert( vec != 0 ) ;
    for(int i=0;i<n; ++i) { vec[i] = t ; }
  }
  // access methods
  inline TYPE const operator() ( int i ) const {
    assert( 0<=i && i<n ) ;
    return vec[i] ;
  }
  inline TYPE & operator() ( int i ) {
    assert( 0<=i && i<n ) ;
    return vec[i] ;
  }
  inline TYPE & back() { 
    assert( n>0 ) ;
    return vec[n-1] ; 
  }
  // size method
  inline int size() const { return n ; }
  // find method (implement a binary search)
  inline bool find( TYPE t ) {
    
    int ia(0), ib(n-1) ;
    int ic = (ia+ib)/2 ;
    
    bool retval = 
      (t==vec[ia]) || (t==vec[ib]) || (t==vec[ic]) ;

    while( !retval && ib-ia>1 ) {
      if      ( t<vec[ic] ) { ib = ic ; }
      else if ( t>vec[ic] ) { ia = ic ; }
      ic = (ia+ib)/2 ;
      retval = 
	(t==vec[ia]) || (t==vec[ib]) || (t==vec[ic]) ;
      
    }
    return retval ;
  }
  // find method (implement a binary search)
  inline int ipos( TYPE t ) {
    
    int ia(0), ib(n-1) ;
    int ic = (ia+ib)/2 ;
    
    bool retval = 
      (t==vec[ia]) || (t==vec[ib]) || (t==vec[ic]) ;

    while( !retval && ib-ia>1 ) {
      if      ( t<vec[ic] ) { ib = ic ; }
      else if ( t>vec[ic] ) { ia = ic ; }
      ic = (ia+ib)/2 ;
      retval = 
	(t==vec[ia]) || (t==vec[ib]) || (t==vec[ic]) ;
      
    }
    int ival = -1 ;
    if      (t==vec[ia]) { ival=ia ; } 
    else if (t==vec[ib]) { ival=ib ; } 
    else if (t==vec[ic]) { ival=ic ; }
    else                 { assert(!retval) ; } // if none, error exit
    return ival ;
  }
} ;

class dataset {

  typedef TDataVec<unsigned int> VecUInt ;
  typedef TDataVec<bool>         VecBool ;

private:
  int nrow ;
  VecUInt * iset ;
  VecBool * bset ;

public:
  dataset() {}
  ~dataset() {}
  void setup( int _nrow ) {
    nrow = _nrow ; 
    iset = new VecUInt[nrow] ;
    bset = new VecBool[nrow] ;
  }
  void setup( int i, int nnz ) {
    assert( 0<=i && i<nrow ) ;
    iset[i].setup(nnz) ;
    bset[i].setup(nnz) ;
  }
  int size() const { 
    return nrow ; 
  }
  int size_loc( int i ) const { 
    assert( 0<=i && i<nrow ) ;
    return iset[i].size() ; 
  }
  //only for accessing indices
  int operator() ( int i, int j ) const {
    assert( 0<=i && i<nrow ) ;
    assert( 0<=j && j<iset[i].size() ) ;
    return iset[i](j) ;
  }
  // "get"-methods
  bool get_index( int i, int j ) const {
    if ( !( 0<=j && j<iset[i].size() ) ) {
      MSG("wrong args in dataset::get_index ") ; 
      VAL( j ) ; PRT( iset[i].size() ) ;
    }
    assert( 0<=i && i<nrow ) ;
    assert( 0<=j && j<iset[i].size() ) ;
    return iset[i](j) ;
  }
  bool get_value( int i, int j ) const {
    if ( !( 0<=j && j<iset[i].size() ) ) {
      MSG("wrong args in dataset::get_value ") ; 
      VAL( j ) ; PRT( iset[i].size() ) ;
    }
    assert( 0<=i && i<nrow ) ;
    assert( 0<=j && j<iset[i].size() ) ;
    return bset[i](j) ;
  }
  // "set"-methods
  void set_index( int i, int j, int ival ) {
    assert( 0<=i && i<nrow ) ;
    assert( 0<=j && j<iset[i].size() ) ;
    iset[i](j) = ival ;
  }
  void set_value( int i, int j, int bval ) {
    assert( 0<=i && i<nrow ) ;
    assert( 0<=j && j<iset[i].size() ) ;
    bset[i](j) = bval ;
  }
} ;

#endif // end of _DATASET_HH
