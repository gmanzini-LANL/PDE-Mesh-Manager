/* ----------------------------------------------------------------/
This is open source software; you can redistribute it and/or modify
it under the terms of the BSD-3 License. If software is modified
to produce derivative works, such modified software should be
clearly marked, so as not to confuse it with the version available
from LANL. Full text of the BSD-3 License can be found in the
License file of the repository.
/---------------------------------------------------------------- */

#ifndef _MESH3D_GMV_OUTPUT_HH
#define _MESH3D_GMV_OUTPUT_HH

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cassert>
#include<string>
using namespace std ;

#include "mesh3D.hh"

class mesh3Dv_GMV_output {
protected:
  mesh_3Dv & mesh ;

private:
  int offset ;

  // auxiliary method
  string integer_to_string( int n, string str_fmt ) ;
  
protected:
  void write_vrtx_data( ostream & file_gmv ) ;
  void write_regn_data( ostream & file_gmv ) ;
  void write_materials( ostream & file_gmv ) ;

protected:
  void write_vrtx_data( ostream & file_gmv, int iR, vector<int> & vrtx_idx ) ;
  void write_regn_data( ostream & file_gmv, int iR, vector<int> & vrtx_idx ) ;
  void write_materials( ostream & file_gmv, int iR ) ;

public:
  mesh3Dv_GMV_output ( mesh_3Dv & _mesh, int _offset=1 ) : 
    mesh(_mesh), offset(_offset) {}
  ~mesh3Dv_GMV_output () {}
  
  int  get_offset()              { return offset ; }
  void set_offset( int _offset ) { offset=_offset ; }

  virtual void write_mesh( string _fname=string("mesh3D") )  ;
  virtual void write_mesh( string _fname, int iR )  ;
} ;

// auxiliary method
string mesh3Dv_GMV_output :: integer_to_string( int n, string str_fmt ) {
  const int n_buf = 10 ;
  char buf[n_buf] ;
  sprintf(buf,str_fmt.c_str(),n);
    return string(buf) ;
}

//-------------- ROUTINES FOR WRITING ALL CELLS
void mesh3Dv_GMV_output :: write_mesh( string fname ) {
  cout << "start  mesh3Dv_GMV_output::write_mesh"<<endl<<flush ;
  cout << "write \""<<fname<<"\""<<endl<<flush ;

  string gmv_fname = fname + string(".gmv") ;
  ofstream file_gmv(gmv_fname.c_str()) ;

  // header of gmv file in ascii format
  file_gmv << "gmvinput ascii" << endl ;

  // write vertex data file
  write_vrtx_data( file_gmv ) ;

  // write region data file
  write_regn_data( file_gmv ) ;

  // write materials (ONLY for GMV output)
  write_materials( file_gmv ) ;
  
  // closure of gmv file in ascii format
  file_gmv << "endgmv" << endl ;
  file_gmv.close() ;

  cout << "end of mesh3Dv_writer_REGN_FACE_format::write_mesh"<<endl<<flush ;
}
void mesh3Dv_GMV_output :: write_vrtx_data( ostream & file_gmv ) {
  // write vertex identifiers, coordinates & info
  file_gmv << "nodev " << mesh.n_vertex() << endl ;
  for ( int iV=0 ; iV<mesh.n_vertex() ; ++iV ) {
    file_gmv << setprecision(16)
	     << "     "  << mesh.coords_V(iV,0) << "   " 
	     << "     "  << mesh.coords_V(iV,1) << "   " 
	     << "     "  << mesh.coords_V(iV,2) << endl ;
  }
  file_gmv << endl ;
}
void mesh3Dv_GMV_output :: write_regn_data( ostream & file_gmv ) {
  // write vertex identifiers
  file_gmv << "cells " << mesh.n_region() << endl << endl ;
  for ( int iR=0 ; iR<mesh.n_region() ; ++iR ) {
    int nRF = mesh.n_regn_face(iR) ;
    file_gmv << "general " << nRF << endl ;
    for ( int ilF=0; ilF<nRF; ++ilF ) {
      int iF = mesh.regn_face( iR, ilF ) ;
      vector<int> face_vlist ;
      mesh.get_face_vrtx( iF, face_vlist ) ;
      int nFV = face_vlist.size() ;
      file_gmv << "  " << nFV ;
    }
    file_gmv << endl << " " ;
    for ( int ilF=0; ilF<nRF; ++ilF ) {
      int iF = mesh.regn_face( iR, ilF ) ;
      vector<int> face_vlist ;
      mesh.get_face_vrtx( iF, face_vlist ) ;
      int nFV = face_vlist.size() ;
      for ( int ilV=0; ilV<nFV; ++ilV ) {
	int iV = face_vlist[ilV] ;
	file_gmv << " " << iV+get_offset() ;
      }
      file_gmv << "  " ;
    }
    file_gmv << endl << endl ;
  }
}
void mesh3Dv_GMV_output :: write_materials( ostream & file_gmv ) {
  // (no offset)
  int n_matr =  0 ;
  int nR = mesh.n_region() ;
  vector<int> regn_matr(nR) ;
  for ( int iR=0; iR<nR; ++iR ) {
    int matr_R = mesh.get_fR(iR) ;
    if ( matr_R > -1 ) {        // should be more correct to say != UNSET
     regn_matr[iR] = matr_R ;
     ++n_matr ;
    }
  }
  // write material field for cells (0)
  file_gmv << "material " << n_matr << " 0 " << endl ;
  for ( int i_matr=0; i_matr<n_matr; ++i_matr  ) {
    file_gmv << string("layer-") + integer_to_string(i_matr+offset,"%d") << endl ;
  }
  int n = 0 ;
  for ( int iR=0; iR<nR; ++iR ) {
    if ( mesh.get_fR(iR) > -1 ) { // should be more correct to say != UNSET
       file_gmv << regn_matr[iR]+offset << " " ;
       if ( n%40==0 ) { file_gmv << endl ; }
    }
  }
}

//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------

//-------------- ROUTINES FOR WRITING A SINGLE CELL
// origin is centered at cell barycenter coordinates
void mesh3Dv_GMV_output :: write_mesh( string fname, int iR ) {
  assert( 0<=iR && iR<mesh.n_region() ) ;
  cout << "start  mesh3Dv_GMV_output::write_mesh"<<endl<<flush ;
  cout << "write \""<<fname<<"\""<<endl<<flush ;

  // vrtx renumbering
  vector<int> vrtx_idx(mesh.n_vertex()) ;
  for ( int iV=0; iV<mesh.n_vertex(); ++iV ) { vrtx_idx[iV] = -1 ; }

  // open GMV file
  string gmv_fname = fname + string(".gmv") ;
  ofstream file_gmv(gmv_fname.c_str()) ;

  // header of gmv file in ascii format
  file_gmv << "gmvinput ascii" << endl ;

  // write vertex data file
  write_vrtx_data( file_gmv, iR, vrtx_idx ) ;

  // write region data file
  write_regn_data( file_gmv, iR, vrtx_idx ) ;

  // write materials (ONLY for GMV output)
  write_materials( file_gmv, iR ) ;
  
  // closure of gmv file in ascii format
  file_gmv << "endgmv" << endl ;
  file_gmv.close() ;

  cout << "end of mesh3Dv_writer_REGN_FACE_format::write_mesh"<<endl<<flush ;
}
void mesh3Dv_GMV_output :: write_vrtx_data( ostream & file_gmv, int iR, vector<int> & vrtx_idx ) {
  // build region iR
  vector<int> regn_vlist ;
  mesh.get_regn_vrtx( iR, regn_vlist ) ;
  int nRV = regn_vlist.size() ;

  // write vertex identifiers, coordinates & info
  
  file_gmv << "nodev " << nRV << endl ;
  for ( int ilV=0 ; ilV<nRV ; ++ilV ) {
    int iV = regn_vlist[ilV] ;
    vrtx_idx[iV] = ilV ;
    file_gmv << setprecision(16) << "  " 
	     << "   "  << mesh.coords_V(iV,0)-mesh.coords_R(iR,0) 
	     << "   "  << mesh.coords_V(iV,1)-mesh.coords_R(iR,1) 
	     << "   "  << mesh.coords_V(iV,2)-mesh.coords_R(iR,2) 
	     << endl ;
  }
  file_gmv << endl ;
}
void mesh3Dv_GMV_output :: write_regn_data( ostream & file_gmv, int iR, vector<int> & vrtx_idx ) {
  // write vertex identifiers
  file_gmv << "cells 1" << endl << endl ;
  int nRF = mesh.n_regn_face( iR ) ;
  file_gmv << "general " << nRF << endl ;
  for ( int ilF=0; ilF<nRF; ++ilF ) {
    int iF = mesh.regn_face( iR, ilF ) ;
    vector<int> face_vlist ;
    mesh.get_face_vrtx( iF, face_vlist ) ;
    int nFV = face_vlist.size() ;
    file_gmv << "  " << nFV ;
  }
  file_gmv << endl << " " ;
  for ( int ilF=0; ilF<nRF; ++ilF ) {
    int iF = mesh.regn_face( iR, ilF ) ;
    vector<int> face_vlist ;
    mesh.get_face_vrtx( iF, face_vlist ) ;
    int nFV = face_vlist.size() ;
    for ( int ilV=0; ilV<nFV; ++ilV ) {
      int iV = face_vlist[ ilV ] ;
      file_gmv << " " << vrtx_idx[ iV ] +get_offset() ;
    }
    file_gmv << "  " ;
  }
  file_gmv << endl << endl ;
}
void mesh3Dv_GMV_output :: write_materials( ostream & file_gmv, int iR ) {
  // write material field for a single cell
  file_gmv << "material 1  0 " << endl ;
  file_gmv << string("layer-1") << endl << "1" << endl ;
}

#endif // end of _MESH3D_GMV_OUTPUT_CC
