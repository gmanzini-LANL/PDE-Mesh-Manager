/* ----------------------------------------------------------------/
This is open source software; you can redistribute it and/or modify
it under the terms of the BSD-3 License. If software is modified
to produce derivative works, such modified software should be
clearly marked, so as not to confuse it with the version available
from LANL. Full text of the BSD-3 License can be found in the
License file of the repository.
/---------------------------------------------------------------- */

#ifndef _MESH_3D_GREADER_HH
#define _MESH_3D_GREADER_HH

#include <iostream>
#include<fstream>
#include <vector>
#include <cassert>
#include<string>
using namespace std ;

//--------------------------------------------------------------------------------------------
// auxiliary functions useful for debugging
//--------------------------------------------------------------------------------------------
void print_logfile( vector<double> & xV, vector<double> & yV,  vector<double> & zV, vector<int> & fV, 
		    vector<int> & regn_vlist, vector<int> & fR, int offset=1 ) {
  
  // open log file
  ofstream LOGF("mesh.log") ;
  
  // print vertex data
  int nV = fV.size() ;
  LOGF << "number of vertices " << nV << endl ;
  for ( int iV=0; iV<nV; ++iV ) {
    LOGF << iV+offset << "  " << xV[iV] << "  " << yV[iV] << "  " << zV[iV] << "  " << fV[iV] << endl ;
  }
  LOGF << "---------------------------------" << endl ;
  
  // print regn data
  int nR = regn_vlist.back() ;
  LOGF << "number of regions " << nR << endl ;
  int kV = 0;
  for ( int iR=0; iR<nR; ++iR ) {
    int nRV = regn_vlist[kV++] ;
    LOGF << nRV << "\t< " ;
    for ( int iRV=0; iRV<nRV-1; ++iRV ) {
      LOGF << regn_vlist[kV++]+offset << ", " ;
    }
    LOGF << regn_vlist[kV++]+offset << " >  " ;
    LOGF << fR[iR] << " " << endl ;
  }

  // close log file
  LOGF.close() ;
}

void print_regn_face( vector<double> & xV, vector<double> & yV,  vector<double> & zV, vector<int> & fV, 
		      vector<int> & regn_flist, vector<int> & fR, int offset=1 ) {

  // open log file
  ofstream LOGF("mesh.log") ;

  // print vertex data
  int nV = fV.size() ;
  LOGF << "number of vertices " << nV << endl ;
  for ( int iV=0; iV<nV; ++iV ) {
    LOGF << iV+offset << "  " << xV[iV] << "  " << yV[iV] << "  " << zV[iV] << "  " << fV[iV] << endl ;
  }
  LOGF << "---------------------------------" << endl ;

  // print regn data
  int nR = regn_flist.back() ;
  LOGF << "number of regions " << nR << endl ;
  int kV = 0;
  for ( int iR=0; iR<nR; ++iR ) {
    int nRF = regn_flist[kV++] ;
    LOGF << " iR=" << iR+offset << " -->nRF=" << nRF << endl ;
    for ( int ilF=0; ilF<nRF; ++ilF ) {
      int nFV = regn_flist[kV++] ;
      LOGF << "  ilF=" << ilF << " -->";
      for ( int ilV=0; ilV<nFV; ++ilV ) {
	LOGF << "  " << regn_flist[kV++]+offset ;
      }
      LOGF << endl ;
    }
  }

  // close log file
  LOGF.close() ;
}

//--------------------------------------------------------------------------------------------
// base class for 3D mesh readers to be used in public derivations
//--------------------------------------------------------------------------------------------

class mesh3D_reader {

protected:
  static istream & eatline(istream & s) {
    while ( s.get() != '\n' && s.good() ) {}
    return s ;
  }
  
  static istream & eatchar(istream & s) { s.get() ; return s ; }
  
  static istream & eatcomments(istream & s) {
    char c = s.peek() ;
    while ( ( c == '!' || c == '%' || c == '#' || c == ';' || c == '$')
	    && s.good() ) { s >> eatline ; c = s.peek() ; }
    return s ;
  }
  
  void error_message( string & s ) {
    cerr << "fatal error:\n" 
	 << "mesh_reader --> read_mesh() cannot open file " << s << endl ;
    exit(0) ;
  }

  virtual void read_vrtx_data( ifstream & input_file, vector<double> & xV, vector<double> & yV, vector<double> & zV, vector<int> & fV )=0 ;
  virtual void read_regn_data( ifstream & input_file, vector<int> & regn_vlist, vector<int> & fR )=0 ;
  
public:
  virtual void read_the_mesh  ( vector<double> & xV, vector<double> & yV,  vector<double> & zV, vector<int> & fV, 
				vector<int> & regn_vlist, vector<int> & fR )=0 ;
} ;

//--------------------------------------------------------------------------------------------
// read input of MESH_3D in TETGEN format
//--------------------------------------------------------------------------------------------
class mesh3D_reader_TETGEN_format : public mesh3D_reader {
      
private:
  string file_name ;
  int offset ; // fortran offset=1, for example tetgen

  int nodelem, nodface ;
  
  // face/region
  void read_fixed_format_int_bnd_markers( ifstream & inp_file, int nRec, int nItm, vector<int> & vlist, vector<int> & flag_list ) ;
  void read_fixed_format_ext_bnd_markers( ifstream & inp_file, int nRec, int nItm, vector<int> & vlist, vector<int> & flag_list ) ;

protected:
  virtual void read_vrtx_data( ifstream & input_file, vector<double> & xV, vector<double> & yV, vector<double> & zV, vector<int> & fV ) ;
  virtual void read_regn_data( ifstream & input_file, vector<int> & regn_vlist, vector<int> & fR ) ;
    
public:
  mesh3D_reader_TETGEN_format( string _file_name, int _offset=1 ) :
    file_name(_file_name), offset(_offset), nodelem(4), nodface(3) {}
  ~mesh3D_reader_TETGEN_format() {}

  virtual void read_the_mesh( vector<double> & xV, vector<double> & yV, vector<double> & zV, vector<int> & fV,
                              vector<int> & regn_vlist, vector<int> & fR ) ;
} ; 
//--------------------------------------------------------------------------------------------
void mesh3D_reader_TETGEN_format :: read_the_mesh ( vector<double> & xV, vector<double> & yV, vector<double> & zV, vector<int> & fV, 
						    vector<int> & regn_vlist, vector<int> & fR ) {

  MSG("start  mesh3D_reader_TETGEN::read_the_mesh"<<endl<<flush) ;
  MSG("read \""<<file_name<<"\""<<endl<<flush) ;
  
  // read .NODE file of tetgen output
  string node_file_name = file_name + ".node" ;
  ifstream node_file( node_file_name.c_str() ) ;
  if ( node_file.good() ) {
    read_vrtx_data( node_file, xV, yV, zV, fV ) ;
  } else {
    error_message(node_file_name) ;
  }
  node_file.close() ;

  // read .ELE file of tetgen output
  string ele_file_name = file_name + ".ele" ;
  ifstream ele_file( ele_file_name.c_str() ) ;
  if ( ele_file.good() ) {
    read_regn_data( ele_file, regn_vlist, fR ) ;
  } else {
    error_message(ele_file_name) ;
  }
  ele_file.close() ;

  MSG("end-->mesh3D_reader_TETGEN::read_the_mesh"<<endl<<endl<<flush) ;
}
//--------------------------------------------------------------------------------------------
void mesh3D_reader_TETGEN_format :: read_vrtx_data ( ifstream & input_file, 
						     vector<double> & xV, vector<double> & yV, vector<double> & zV, vector<int> & fV ) {
  MSG("start  read_vrtx_data"<<endl<<flush) ;
  
  // declarations
  double x(0.), y(0.), z(0.) ;
  int id(0), nV(0), ndim(0), n_attr(0), nbmrk(0), f(0) ;
  
  // read file header
  input_file >> eatcomments >> nV >> ndim >> n_attr >> nbmrk >> eatline ;
  assert( ndim==3 ) ;

  // resize the vertex list of mesh
  MSG("---#vertices: ") ; PRT( nV ) ;
  // read mesh vertices
  if ( nbmrk==0 ) { 
    for ( int iV=0 ; iV<nV ; ++iV ) {
      input_file >> eatcomments >> id >> x >> y >> z >> eatline ;
      assert( id-offset==iV ) ;
      xV.push_back( x ) ;
      yV.push_back( y ) ;
      zV.push_back( z ) ;
      fV.push_back( UNSET ) ;
    }
  } else if ( nbmrk==1 ) { 
    for ( int iV=0 ; iV<nV ; ++iV ) {
      input_file >> eatcomments >> id >> x >> y >> z >> f >> eatline ;
      assert( id-offset==iV ) ;
      xV.push_back( x ) ;
      yV.push_back( y ) ;
      zV.push_back( z ) ;
      fV.push_back( f ) ;
    } 
  } else {
    assert(0) ;
  }
  MSG("end of read_vrtx_data"<<endl) ;
}
//--------------------------------------------------------------------------------------------
void mesh3D_reader_TETGEN_format :: 
read_regn_data( ifstream & ele_file, vector<int> & regn_vlist, vector<int> & fR ) {
  MSG("start  read_regn_data"<<endl<<flush) ;

  // declarations
  int nR=0, id=0, nbmrk=0 ;
  int kV(0) ;
  
  // read file's header record
  ele_file >> eatcomments >> nR >> nodelem >> nbmrk >> eatline ;

  // read the vertex region list of mesh
  MSG("---#regions: ") ; PRT(nR) ;
  if      ( nodelem>0 && nbmrk==0 ) { read_fixed_format_int_bnd_markers( ele_file, nR, nodelem, regn_vlist, fR ) ; }
  else if ( nodelem>0 && nbmrk==1 ) { read_fixed_format_ext_bnd_markers( ele_file, nR, nodelem, regn_vlist, fR ) ; }
  else { assert(0) ; }
  regn_vlist.push_back( nR ) ;

  MSG("end of read_regn_data"<<endl<<flush) ;
}
//--------------------------------------------------------------------------------------------
void mesh3D_reader_TETGEN_format :: 
read_fixed_format_int_bnd_markers( ifstream & inp_file, int nRec, int nVrt, vector<int> & vlist, vector<int> & flag_list ) {
  int id(0), kV(0) ;
  for ( int iR=0 ; iR<nRec ; ++iR ) {
    inp_file >> eatcomments >> id ;
    vlist.push_back( nVrt ) ;
    for ( int k=0; k<nVrt; ++k ) {
      inp_file >> kV ;
      vlist.push_back( kV-offset) ; // offset = 1 (FORTRAN style)
    }
    inp_file >> eatline ;
    flag_list.push_back( UNSET ) ;
  }
}
//--------------------------------------------------------------------------------------------
void mesh3D_reader_TETGEN_format :: 
read_fixed_format_ext_bnd_markers( ifstream & inp_file, int nRec, int nVrt, vector<int> & vlist, vector<int> & flag_list ) {
  int id(0), kV(0), ext_flag(UNSET) ;
  for ( int iR=0 ; iR<nRec ; ++iR ) {
    inp_file >> eatcomments >> id ;
    vlist.push_back( nVrt ) ;
    for ( int k=0; k<nVrt; ++k ) {
      inp_file >> kV ;
      vlist.push_back( kV-offset) ; // offset = 1 (FORTRAN style)
    }
    inp_file >> ext_flag >> eatline ;
    flag_list.push_back( ext_flag ) ;
  }
  vlist.push_back( nRec ) ;
}
//--------------------------------------------------------------------------------------------
class mesh3D_reader_REGN_FACE_format : public mesh3D_reader {
  
private:
  string file_name ;
  int offset ; // fortran offset=1, for example tetgen

  int nodelem, nodface ;

  // face/region
  void read_fixed_format_int_bnd_markers( ifstream & inp_file, int nRec, vector<int> & flist, vector<int> & flag_list ) ;
  void read_fixed_format_ext_bnd_markers( ifstream & inp_file, int nRec, vector<int> & flist, vector<int> & flag_list ) ;
  // face/region
  void read_fixed_format_int_bnd_markers_reversed_face( ifstream & inp_file, int nRec, vector<int> & flist, vector<int> & flag_list ) ;
  void read_fixed_format_ext_bnd_markers_reversed_face( ifstream & inp_file, int nRec, vector<int> & flist, vector<int> & flag_list ) ;

protected:
  virtual void read_vrtx_data( ifstream & input_file, vector<double> & xV, vector<double> & yV, vector<double> & zV, vector<int> & fV ) ;
  virtual void read_regn_data( ifstream & input_file, vector<int> & regn_flist, vector<int> & fR ) ;
  
public:
  mesh3D_reader_REGN_FACE_format( string _file_name, int _offset=1 ) : 
    file_name(_file_name), offset(_offset), nodelem(4), nodface(3) {}
  ~mesh3D_reader_REGN_FACE_format() {}
  
  virtual void read_the_mesh( vector<double> & xV, vector<double> & yV, vector<double> & zV, vector<int> & fV, 
			      vector<int> & regn_flist, vector<int> & fR ) ;

} ;
void mesh3D_reader_REGN_FACE_format :: read_the_mesh ( vector<double> & xV, vector<double> & yV, vector<double> & zV, vector<int> & fV, 
						       vector<int> & regn_flist, vector<int> & fR ) {

  MSG("start  mesh3D_reader_REGN_FACE_format::read_the_mesh"<<endl<<flush) ;
  MSG("read \""<<file_name<<"\""<<endl<<flush) ;

  // read .NODE file of tetgen output
  string node_file_name = file_name + ".node" ;
  ifstream node_file( node_file_name.c_str() ) ;
  if ( node_file.good() ) {
    read_vrtx_data( node_file, xV, yV, zV, fV ) ;
  } else {
    error_message(node_file_name) ;
  }
  node_file.close() ;

  // read .ELE file of tetgen output
  string ele_file_name = file_name + ".ele" ;
  ifstream ele_file( ele_file_name.c_str() ) ;
  if ( ele_file.good() ) {
    read_regn_data( ele_file, regn_flist, fR ) ;
  } else {
    error_message(ele_file_name) ;
  }
  ele_file.close() ;
  
  // screen output
  MSG("end-->mesh3D_reader_REGN_FACE_format::read_the_mesh"<<endl<<endl<<flush) ;
}
//--------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------
void mesh3D_reader_REGN_FACE_format :: read_vrtx_data ( ifstream & input_file, 
							vector<double> & xV, 
							vector<double> & yV,
							vector<double> & zV,
							vector<int>    & fV ) {
  MSG("start mesh3D_reader_REGN_FACE_format::read_vrtx_data"<<endl) ;
  
  // declarations
  double x(0.), y(0.), z(0.) ;
  int id(0), nV(0), ndim(0), n_attr(0), nbmrk(0), f(0) ;
  
  // read file header
  input_file >> eatcomments >> nV >> ndim >> n_attr >> nbmrk >> eatline ;
  assert( ndim==3 ) ;

  // resize the vertex list of mesh
  MSG("---#vertices: ") ; PRT( nV ) ;
  // read mesh vertices
  if ( nbmrk==0 ) { 
    for ( int iV=0 ; iV<nV ; ++iV ) {
      input_file >> eatcomments >> id >> x >> y >> z >> eatline ;
      assert( id-offset==iV ) ;
      xV.push_back( x ) ;
      yV.push_back( y ) ;
      zV.push_back( z ) ;
      fV.push_back( UNSET ) ;
    }
  } else if ( nbmrk==1 ) { 
    for ( int iV=0 ; iV<nV ; ++iV ) {
      input_file >> eatcomments >> id >> x >> y >> z >> f >> eatline ;
      assert( id-offset==iV ) ;
      xV.push_back( x ) ;
      yV.push_back( y ) ;
      zV.push_back( z ) ;
      fV.push_back( f ) ;
    } 
  } else {
    assert(0) ;
  }
  MSG("end-->mesh_reader_REGN_FACE_format::read_vrtx_data"<<endl) ;
}

// nR  nbmrk
// iR  nRF  [fR]
// ilF nFV  V0 V1 V2
void mesh3D_reader_REGN_FACE_format :: read_regn_data( ifstream    & ele_file, 
						       vector<int> & regn_flist,
						       vector<int> & fR ) {
  MSG("start read_regn_data"<<endl) ;

  // declarations
  int nR=0, id=0, nbmrk=0 ;
  int kV(0) ;
  
  // read file's header record
  ele_file >> eatcomments >> nR >> nbmrk >> eatline ;

#if 1
  // read the vertex region list of mesh
  MSG("---#regions: ") ; PRT(nR) ;
  if      ( nbmrk==0 ) { read_fixed_format_int_bnd_markers( ele_file, nR, regn_flist, fR ) ; }
  else if ( nbmrk==1 ) { read_fixed_format_ext_bnd_markers( ele_file, nR, regn_flist, fR ) ; }
  else { assert(0) ; }
  regn_flist.push_back( nR ) ;
#else
  // read the vertex region list of mesh, reverse the face orientation!!!
  // to be used for the meshes from KL
  MSG("---#regions: ") ; PRT(nR) ;
  if      ( nbmrk==0 ) { read_fixed_format_int_bnd_markers_reversed_face( ele_file, nR, regn_flist, fR ) ; }
  else if ( nbmrk==1 ) { read_fixed_format_ext_bnd_markers_reversed_face( ele_file, nR, regn_flist, fR ) ; }
  else { assert(0) ; }
  regn_flist.push_back( nR ) ;
#endif

  MSG("end-->read_regn_data"<<endl) ;
}
void mesh3D_reader_REGN_FACE_format :: 
read_fixed_format_int_bnd_markers( ifstream & inp_file, int nRec, vector<int> & flist, vector<int> & flag_list ) {
  int R_ig(0), F_il(0), kV(0) ;
  int nRF(0), nFV(0) ;
  for ( int iR=0 ; iR<nRec ; ++iR ) {
    inp_file >> eatcomments >> R_ig >> nRF ;
    flist.push_back( nRF ) ;
    for ( int ilF=0; ilF<nRF; ++ilF ) {
      inp_file >> eatcomments >> F_il >> nFV ;
      flist.push_back( nFV ) ;
      for ( int ilV=0; ilV<nFV; ++ilV ) {
	inp_file >> kV ;
	flist.push_back( kV-offset) ; // offset = 1 (FORTRAN style)
      }
    }
    inp_file >> eatline ;
    flag_list.push_back( UNSET ) ;
  }
}
void mesh3D_reader_REGN_FACE_format :: 
read_fixed_format_ext_bnd_markers( ifstream & inp_file, int nRec, vector<int> & flist, vector<int> & flag_list ) {
  int R_ig(0), F_il(0), kV(0) ;
  int nRF(0), nFV(0), ext_flag(0) ;
  for ( int iR=0 ; iR<nRec ; ++iR ) {
    inp_file >> eatcomments >> R_ig >> nRF >> ext_flag ;
    flist.push_back( nRF ) ;
    for ( int ilF=0; ilF<nRF; ++ilF ) {
      inp_file >> eatcomments >> F_il >> nFV ;
      flist.push_back( nFV ) ;
      for ( int ilV=0; ilV<nFV; ++ilV ) {
	inp_file >> kV ;
	flist.push_back( kV-offset) ; // offset = 1 (FORTRAN style)
      }
    }
    flag_list.push_back( ext_flag ) ;
  }
}
//--------------------------------------------------------------------------------------------
void mesh3D_reader_REGN_FACE_format :: 
read_fixed_format_int_bnd_markers_reversed_face( ifstream & inp_file, int nRec, vector<int> & flist, vector<int> & flag_list ) {
  int R_ig(0), F_il(0), kV(0) ;
  int nRF(0), nFV(0) ;
  vector<int> tmp ;
  for ( int iR=0 ; iR<nRec ; ++iR ) {
    inp_file >> eatcomments >> R_ig >> nRF ;
    flist.push_back( nRF ) ;
    for ( int ilF=0; ilF<nRF; ++ilF ) {
      inp_file >> eatcomments >> F_il >> nFV ;
      flist.push_back( nFV ) ;
      tmp.resize(nFV) ;
      for ( int ilV=0; ilV<nFV; ++ilV ) {
	inp_file >> kV ;
	tmp[ilV] = kV-offset ; // offset = 1 (FORTRAN style)
      }
      for ( int ilV=nFV-1; ilV>=0; --ilV ) {
	flist.push_back( tmp[ilV] ) ;
      }
    }
    inp_file >> eatline ;
    flag_list.push_back( UNSET ) ;
  }
}
void mesh3D_reader_REGN_FACE_format :: 
read_fixed_format_ext_bnd_markers_reversed_face( ifstream & inp_file, int nRec, vector<int> & flist, vector<int> & flag_list ) {
  int R_ig(0), F_il(0), kV(0) ;
  int nRF(0), nFV(0), ext_flag(0) ;
  vector<int> tmp ;
  for ( int iR=0 ; iR<nRec ; ++iR ) {
    inp_file >> eatcomments >> R_ig >> nRF >> ext_flag ;
    flist.push_back( nRF ) ;
    for ( int ilF=0; ilF<nRF; ++ilF ) {
      inp_file >> eatcomments >> F_il >> nFV ;
      flist.push_back( nFV ) ;
      tmp.resize(nFV) ;
      for ( int ilV=0; ilV<nFV; ++ilV ) {
	inp_file >> kV ;
	tmp[ilV] = kV-offset ; // offset = 1 (FORTRAN style)
      }
      for ( int ilV=nFV-1; ilV>=0; --ilV ) {
	flist.push_back( tmp[ilV] ) ;
      }
    }
    flag_list.push_back( ext_flag ) ;
  }
}

//--------------------------------------------------------------------------------------------
// read input of MESH_3D in MSH format
// actually, only a sub-set of the records of a *.msh file are used
//--------------------------------------------------------------------------------------------
class mesh3D_reader_MSH_format : public mesh3D_reader {
      
private:
  string file_name ;
  int offset ; // fortran offset=1, for example tetgen
  
  // mesh sizes (available in file's header)
  int nV, nR, nF, nE ;
  
  // read file's header
  void read_header( ifstream & input_file ) ;
  void print_vector( vector< vector<int> > & vec, int nvec ) ;
  
protected:
  virtual void read_vrtx_data( ifstream & input_file, vector<double> & xV, vector<double> & yV, vector<double> & zV, vector<int> & fV ) ;
  virtual void read_regn_data( ifstream & input_file, vector<int> & regn_vlist, vector<int> & fR ) ;
    
public:
  mesh3D_reader_MSH_format( string _file_name, int _offset=1 ) :
    file_name(_file_name), offset(_offset) {}
  ~mesh3D_reader_MSH_format() {}

  virtual void read_the_mesh( vector<double> & xV, vector<double> & yV, vector<double> & zV, vector<int> & fV,
                              vector<int> & regn_vlist, vector<int> & fR ) ;
} ; 
//--------------------------------------------------------------------------------------------
void mesh3D_reader_MSH_format :: 
read_the_mesh ( vector<double> & xV, vector<double> & yV, vector<double> & zV, vector<int> & fV,                    
		vector<int> & regn_flist, vector<int> & fR ) {
  
  MSG("start  mesh3D_reader_MSH::read_the_mesh"<<endl<<flush) ;
  MSG("read \""<<file_name<<"\""<<endl<<flush) ;
  
  PRT(offset) ;
  
  // open *.msh file
  string msh_file_name = file_name + ".msh" ;
  ifstream msh_file( msh_file_name.c_str() ) ;
  if ( msh_file.good() ) {
    read_header( msh_file ) ;
    read_vrtx_data ( msh_file, xV, yV, zV, fV ) ;
    read_regn_data ( msh_file, regn_flist, fR ) ;
  } else {
    error_message( msh_file_name ) ;
  }
  msh_file.close() ;
  
  MSG("end-->mesh3D_reader_MSH::read_the_mesh"<<endl<<endl<<flush) ;
}
//--------------------------------------------------------------------------------------------
void  mesh3D_reader_MSH_format ::
read_header( ifstream & input_file ) {
  string str_dummy, str_gridname ;
  int n_version ;
  
  input_file >> str_dummy >> eatline >> eatline >> eatline ;
  input_file >> str_dummy >> eatline >> n_version    >> eatline ;
  input_file >> str_dummy >> eatline >> str_gridname >> eatline ;
  input_file >> str_dummy >> eatline ;
  input_file >> str_dummy >> eatline >> nV >> eatline ;
  input_file >> str_dummy >> eatline >> nR >> eatline ;
  input_file >> str_dummy >> eatline >> nF >> eatline ;
  input_file >> str_dummy >> eatline >> nE >> eatline ;

  PRT(n_version) ;
  PRT(str_gridname) ;
  VAL(nV) ; VAL(nR) ; VAL(nF) ; PRT(nE) ;
}
//--------------------------------------------------------------------------------------------
void mesh3D_reader_MSH_format ::
read_vrtx_data ( ifstream & input_file, 
		 vector<double> & xV, vector<double> & yV, vector<double> & zV, vector<int> & fV ) {

  MSG("start  read_vrtx_data"<<endl<<flush) ;
  
  // declarations
  double x(0.), y(0.), z(0.) ;
  int    id(0) ;

  // control print
  MSG("---#vertices: ") ; PRT( nV ) ;

  // read mesh vertices
  input_file >> eatline ;
  for ( int iV=0 ; iV<nV ; ++iV ) {
    input_file >> x >> y >> z >> eatline ;
    xV.push_back( x ) ;
    yV.push_back( y ) ;
    zV.push_back( z ) ;
    fV.push_back( UNSET ) ;
  }
}
//--------------------------------------------------------------------------------------------
// re-add offset to compare with input data
void mesh3D_reader_MSH_format ::
print_vector( vector< vector<int> > & vec, int nvec ) {
  LINE(---) ;
  for ( int i=0; i<nvec; ++i ) {
    int n = vec[i].size() ;
    cout << i+offset << "-> " << n << " <- " ;
    for ( int il=0; il<n; ++il ) {
      cout << vec[i][il]+offset << "  " ;
    }
    cout << endl ;
  }
}
//--------------------------------------------------------------------------------------------
void mesh3D_reader_MSH_format ::
read_regn_data( ifstream & input_file, vector<int> & regn_flist, vector<int> & fR ) {

  string str_dummy("") ;

  // aux STL vectors
  vector< vector<int> > R_flist(nR) ; // *
  //vector<int> R_vlist[nR] ;
  vector< vector<int> > F_vlist(nF) ; // *
  //vector<int> F_elist[nF] ;
  vector< vector<int> > F_rlist(nF) ; // *
  //vector<int> E_vlist[nE] ;

  bool log_msh = false ;
  ofstream LOG("./msh.log") ;

  { // Volumes->faces
    input_file >> str_dummy >> eatline ;
    if ( log_msh ) { LOG << str_dummy << endl ; }
    LINE(---) ; 
    MSG("MSH reading: " << str_dummy << endl) ;
    for ( int iR=0 ; iR<nR ; ++iR ) {
      int nRF(-1) ;
      input_file >> nRF ;
      if ( log_msh ) { LOG << nRF << "-->" ; }
      for ( int ilF=0; ilF<nRF; ++ilF ) {
	int iF(-1) ;
	input_file >> iF ;
	if ( log_msh ) { LOG << "  " << iF ; }
	R_flist[iR].push_back( iF-offset) ; // offset = 1 (FORTRAN style)
      }
      input_file >> eatline ;
      if ( log_msh ) { LOG << endl ; }
    }
  }
  
  { // Volumes->sommets
    input_file >> str_dummy >> eatline ;
    if ( log_msh ) { LOG << str_dummy << endl ; }
    LINE(---) ; 
    MSG("MSH reading: " << str_dummy << endl) ;
    for ( int iR=0 ; iR<nR ; ++iR ) {
      int nRV(-1) ;
      input_file >> nRV ;
      for ( int ilV=0; ilV<nRV; ++ilV ) {
	int iV(-1) ;
	input_file >> iV ;
	if ( log_msh ) { LOG << "  " << iV ; }
	//R_vlist[iR].push_back( iV-offset) ; // offset = 1 (FORTRAN style)
      }
      input_file >> eatline ;
      if ( log_msh ) { LOG << endl ; }
    }
  }
  
  { // Faces->aretes
    input_file >> str_dummy >> eatline ;
    if ( log_msh ) { LOG << str_dummy << endl ; }
    LINE(---) ; 
    MSG("MSH reading: " << str_dummy << endl) ;
    for ( int iF=0; iF<nF; ++iF ) { 
      int nFE(-1) ;
      input_file >> nFE ;
      for ( int ilE=0; ilE<nFE; ++ilE ) {
	int iE(-1) ;
	input_file >> iE ;
	if ( log_msh ) { LOG << "  " << iE ; }
	//F_elist[iF].push_back( iE-offset) ; // offset = 1 (FORTRAN style)
      }
      input_file >> eatline ;
      if ( log_msh ) { LOG << endl ; }
    }
  }
  
  { // Faces->sommets
    input_file >> str_dummy >> eatline ;
    if ( log_msh ) { LOG << str_dummy << endl ; }
    LINE(---) ; 
    MSG("MSH reading: " << str_dummy << endl) ;
    for ( int iF=0 ; iF<nF ; ++iF ) {
      int nFV(-1) ;
      input_file >> nFV ;
      for ( int ilV=0; ilV<nFV; ++ilV ) {
	int iV(-1) ;
	input_file >> iV ;
	if ( log_msh ) { LOG << "  " << iV ; }
	F_vlist[iF].push_back( iV-offset) ; // offset = 1 (FORTRAN style)
      }
      input_file >> eatline ;
      if ( log_msh ) { LOG << endl ; }
    }
  }

  { // Faces->volumes
    input_file >> str_dummy >> eatline ;
    if ( log_msh ) { LOG << str_dummy << endl ; }
    LINE(---) ; 
    MSG("MSH reading: " << str_dummy << endl) ;
    for ( int iF=0 ; iF<nF ; ++iF ) {
      for ( int ilR=0; ilR<2; ++ilR ) {
	int iR(-1) ;
	input_file >> iR ;
	if ( log_msh ) { LOG << "  " << iR ; }
	F_rlist[iF].push_back( iR-offset) ; // offset = 1 (FORTRAN style)
      }
      input_file >> eatline ;
      if ( log_msh ) { LOG << endl ; }
    }
  }

  { // Aretes->sommets
    input_file >> str_dummy >> eatline ;
    if ( log_msh ) { LOG << str_dummy << endl ; }
    LINE(---) ; 
    MSG("MSH reading: " << str_dummy << endl) ;
    for ( int iE=0 ; iE<nE ; ++iE ) {
      for ( int ilV=0; ilV<2; ++ilV ) {
	int iV(-1) ;
	input_file >> iV ;
	if ( log_msh ) { LOG << "  " << iV ; }
	//E_vlist[iE].push_back( iV-offset) ; // offset = 1 (FORTRAN style)
      }
      input_file >> eatline ;
      if ( log_msh ) { LOG << endl ; }
    }
  }
  
  // print to compare with input data
  // (useful for debugging)
  if ( 0 ) {
    print_vector( R_flist, nR ) ;
    //print_vector( R_vlist, nR ) ;
    print_vector( F_vlist, nF ) ;
    //print_vector( F_elist, nF ) ;
    print_vector( F_rlist, nF ) ;
    //print_vector( E_vlist, nF ) ;
  }

  //print_vector( F_vlist, nF ) ;

  // reconstruct regn-face structure
  for ( int iR=0; iR<nR; ++iR ) {
    int nRF = R_flist[iR].size() ;
    regn_flist.push_back( nRF ) ;
    for ( int ilF=0; ilF<nRF; ++ilF ) {
      int iF = R_flist[iR][ilF] ;
      int  nFV = F_vlist[iF].size() ;
      regn_flist.push_back( nFV ) ;
      bool ok_regn_face = F_rlist[iF][0] == iR ;
      if ( ok_regn_face ) {
	for ( int ilV=0; ilV<nFV; ++ilV ) {
	  int iV = F_vlist[iF][ilV] ;
	  regn_flist.push_back( iV ) ;
	}
      } else {
	for ( int ilV=0; ilV<nFV; ++ilV ) {
	  int iV = F_vlist[iF][nFV-1-ilV] ;
	  regn_flist.push_back( iV ) ;
	}
      }
    }
    fR.push_back( UNSET ) ;
  }
  regn_flist.push_back(nR) ;
}

//--------------------------------------------------------------------------------------------
// base class for mesh readers for public derivations
//--------------------------------------------------------------------------------------------
class mesh2D_reader {
  
protected: 
  static istream & eatline(istream & s) {
    while ( s.get() != '\n' && s.good() ) {}
    return s ;
  } 
    
  static istream & eatchar(istream & s) { s.get() ; return s ; }
      
  static istream & eatcomments(istream & s) {
    char c = s.peek() ;
    while ( ( c == '!' || c == '%' || c == '#' || c == ';' || c == '$')
            && s.good() ) { s >> eatline ; c = s.peek() ; }
    return s ;
  }
  
  void error_message( string & s ) {
    cerr << "fatal error:\n"
         << "mesh_reader --> read_mesh() cannot open file " << s << endl ;
    exit(0) ;         
  }
  
  virtual void read_vrtx_data( ifstream & input_file, vector<double> & xV, vector<double> & yV, vector<int> & fV )=0 ;
  virtual void read_regn_data( ifstream & input_file, vector<int> & regn_vlist, vector<int> & fR )=0 ;
  
public:
  virtual void read_mesh( vector<double> & xV, vector<double> & yV, vector<int> & fV, vector<int> & regn_vlist, vector<int> & fR )=0 ;
} ;

//--------------------------------------------------------------------------------------------
// read input of MESH_2D in "General Format", which is a generalization of TRIANGLE format
//--------------------------------------------------------------------------------------------

class mesh2D_reader_GeneralFormat : public mesh2D_reader {
  
private:
  string file_name ;
  int offset ; // fortran offset

  int nodelem ;

  void read_regn_free_format_int_bnd_markers ( ifstream & ele_file, int nR, vector<int> & regn_vlist, vector<int> & fR ) ;
  void read_regn_free_format_ext_bnd_markers ( ifstream & ele_file, int nR, vector<int> & regn_vlist, vector<int> & fR ) ;
  void read_regn_fixed_format_int_bnd_markers( ifstream & ele_file, int nR, vector<int> & regn_vlist, vector<int> & fR ) ;
  void read_regn_fixed_format_ext_bnd_markers( ifstream & ele_file, int nR, vector<int> & regn_vlist, vector<int> & fR ) ;

protected:
  virtual void read_vrtx_data( ifstream & input_file, vector<double> & xV, vector<double> & yV, vector<int> & fV ) ;
  virtual void read_regn_data( ifstream & input_file, vector<int> & regn_vlist, vector<int> & fR ) ;
  
public:
  mesh2D_reader_GeneralFormat( string _file_name, int _offset=1 ) : file_name(_file_name), offset(_offset) {}
  ~mesh2D_reader_GeneralFormat() {}
  
  virtual void read_mesh( vector<double> & xV, vector<double> & yV, vector<int> & fV, vector<int> & regn_vlist, vector<int> & fR ) ;
} ;
void mesh2D_reader_GeneralFormat :: read_mesh ( vector<double> & xV, vector<double> & yV, vector<int> & fV, 
						vector<int> & regn_vlist, vector<int> & fR ) {
  // read NODE file of triangle output
  string node_file_name = file_name + ".node" ;
  ifstream node_file( node_file_name.c_str() ) ;
  if ( node_file.good() ) {
    read_vrtx_data( node_file, xV, yV, fV ) ;
  } else {
    error_message(node_file_name) ;
  }
  node_file.close() ;
  // read ELE file of triangle output
  string ele_file_name = file_name + ".ele" ;
  ifstream ele_file( ele_file_name.c_str() ) ;
  if ( ele_file.good() ) {
    read_regn_data( ele_file, regn_vlist, fR ) ;
  } else {
    error_message(ele_file_name) ;
  }
  ele_file.close() ;
}
void mesh2D_reader_GeneralFormat :: read_vrtx_data ( ifstream & input_file, 
						     vector<double> & xV, 
						     vector<double> & yV,
						     vector<int>    & fV ) {
  MSG("start mesh2D_reader_GeneralFormat::read_vrtx_data"<<endl) ;

  // declarations
  double x(0.), y(0.) ;
  int id(0), nV(0), ndim(0), n_attr(0), nbmrk(0), f(0) ;
  
  // read file header
  input_file >> eatcomments >> nV >> ndim >> n_attr >> nbmrk >> eatline ;
  
  // resize the vertex list of mesh
  MSG("---#vertices: ") ; PRT( nV ) ;
  // read mesh vertices
  if ( nbmrk==0 ) { 
    for ( int iV=0 ; iV<nV ; ++iV ) {
      input_file >> eatcomments >> id >> x >> y >> eatline ;
      assert( id-offset==iV ) ;
      xV.push_back( x ) ;
      yV.push_back( y ) ;
      fV.push_back( UNSET ) ;
    }
  } else if ( nbmrk==1 ) { 
    for ( int iV=0 ; iV<nV ; ++iV ) {
      input_file >> eatcomments >> id >> x >> y >> f >> eatline ;
      assert( id-offset==iV ) ;
      xV.push_back( x ) ;
      yV.push_back( y ) ;
      fV.push_back( f ) ;
    } 
  } else {
    assert(0) ;
  }
  MSG("end of mesh_reader_GeneralFormat::read_vrtx_data"<<endl) ;
}
void mesh2D_reader_GeneralFormat :: 
read_regn_free_format_int_bnd_markers( ifstream & ele_file, int nR, vector<int> & regn_vlist, vector<int> & fR ) {
  int id(0), nRV(0), kV(0) ;
  for ( int iR=0 ; iR<nR ; ++iR ) {
    ele_file >> eatcomments >> id >> nRV ;
    regn_vlist.push_back( nRV ) ;
    for ( int k=0; k<nRV; ++k ) {
      ele_file >> kV ;
      regn_vlist.push_back( kV-offset) ; // offset = 1 (FORTRAN style)
    }
    ele_file >> eatline ;
    fR.push_back( UNSET ) ;
  }
}
void mesh2D_reader_GeneralFormat ::
read_regn_free_format_ext_bnd_markers( ifstream & ele_file, int nR, vector<int> & regn_vlist, vector<int> & fR ) {
  int id(0), nRV(0), kV(0), regn_flag(UNSET) ;
  for ( int iR=0 ; iR<nR ; ++iR ) {
    ele_file >> eatcomments >> id >> nRV ;
    regn_vlist.push_back( nRV ) ;
    for ( int k=0; k<nRV; ++k ) {
      ele_file >> kV ;
      regn_vlist.push_back( kV-offset) ; // offset = 1 (FORTRAN style)
    }
    ele_file >> regn_flag >> eatline ;
    fR.push_back( regn_flag ) ;
  }
  regn_vlist.push_back( nR ) ;
}
void mesh2D_reader_GeneralFormat :: 
read_regn_fixed_format_int_bnd_markers( ifstream & ele_file, int nR, vector<int> & regn_vlist, vector<int> & fR ) {
  int id(0), kV(0) ;
  for ( int iR=0 ; iR<nR ; ++iR ) {
    ele_file >> eatcomments >> id ;
    regn_vlist.push_back( nodelem ) ;
    for ( int k=0; k<nodelem; ++k ) {
      ele_file >> kV ;
      regn_vlist.push_back( kV-offset) ; // offset = 1 (FORTRAN style)
    }
    ele_file >> eatline ;
    fR.push_back( UNSET ) ;
  }
}
void mesh2D_reader_GeneralFormat :: 
read_regn_fixed_format_ext_bnd_markers( ifstream & ele_file, int nR, vector<int> & regn_vlist, vector<int> & fR ) {
  int id(0), kV(0), regn_flag(UNSET) ;
  for ( int iR=0 ; iR<nR ; ++iR ) {
    ele_file >> eatcomments >> id ;
    regn_vlist.push_back( nodelem ) ;
    for ( int k=0; k<nodelem; ++k ) {
      ele_file >> kV ;
      regn_vlist.push_back( kV-offset) ; // offset = 1 (FORTRAN style)
    }
    ele_file >> regn_flag >> eatline ;
    fR.push_back( regn_flag ) ;
  }
  regn_vlist.push_back( nR ) ;
}
void mesh2D_reader_GeneralFormat :: read_regn_data( ifstream    & ele_file, 
						    vector<int> & regn_vlist,
						    vector<int> & fR ) {
  MSG("start read_regn_data"<<endl) ;

  // declarations
  int nR=0, id=0, nbmrk=0 ;
  int kV(0) ;
  
  // read file header
  ele_file >> eatcomments >> nR >> nodelem >> nbmrk >> eatline ;

  // read the vertex region list of mesh
  MSG("---#regions: ") ; PRT(nR) ;
  if      ( nodelem==0 && nbmrk==0 ) { read_regn_free_format_int_bnd_markers ( ele_file, nR, regn_vlist, fR ) ; }
  else if ( nodelem==0 && nbmrk==1 ) { read_regn_free_format_ext_bnd_markers ( ele_file, nR, regn_vlist, fR ) ; }
  else if ( nodelem >0 && nbmrk==0 ) { read_regn_fixed_format_int_bnd_markers( ele_file, nR, regn_vlist, fR ) ; }
  else if ( nodelem >0 && nbmrk==1 ) { read_regn_fixed_format_ext_bnd_markers( ele_file, nR, regn_vlist, fR ) ; }
  else { assert(0) ; }
  regn_vlist.push_back( nR ) ;

  MSG("end of read_regn_data"<<endl) ;
}

//--------------------------------------------------------------------------------------------
// To add a new mesh reader/builder:
// (i)  write the class for reader object
// (ii) write the corresponding method in ExtFileInput
//--------------------------------------------------------------------------------------------
// MESH READER/BUILDER
// objects of this class perform two actions:
// ------------------------------------------
// (i)  read  different mesh formats
// (ii) build the mesh into the input mesh reference
//--------------------------------------------------------------------------------------------
class ExtFileInput {
private:
  static const int C_offset = 0 ; // C/C++ offset
  static const int F_offset = 1 ; // fortran offset
  void set_regn_face_list( vector<int> & regn_flist, vector<int> & regn_vlist ) ;

public:
  ExtFileInput() {}
  ~ExtFileInput() {}
  
  // readers for different mesh formats, mesh builders
  void TETGEN_format      ( mesh_3Dv & mesh, string fname=string("mesh"), int offset=C_offset ) ;
  void REGN_FACE_format   ( mesh_3Dv & mesh, string fname=string("mesh"), int offset=C_offset ) ;
  void MSH_format         ( mesh_3Dv & mesh, string fname=string("mesh"), int offset=F_offset ) ;
  void multi_layer_format ( mesh_3Dv & mesh, string fname=string("mesh"), int n_layer=1, int offset=C_offset, bool tilted_planes=bool(false) ) ;
  
} ;
// build the face list as in regn_face format for the TETGEN tetrahedra
void ExtFileInput :: set_regn_face_list( vector<int> & regn_flist, vector<int> & regn_vlist ) {

  const int nRF = 4 ;
  const int nFV = 3 ;

  vector<int> vlist(nRF) ;
  int perm[nRF][nFV] = { {0,2,1}, {1,2,3}, {0,1,3}, {0,3,2} } ; 

  int nR  = regn_vlist.back() ;
  
  int k = 0;
  for ( int iR=0; iR<nR; ++iR ) { // loop on the regions

    int nRV = regn_vlist[k++] ;
    assert( nRV==nRF ) ;

    for ( int ilV=0; ilV<nRV; ++ilV ) {
      vlist[ilV] = regn_vlist[k++] ;
    }
    
    regn_flist.push_back(nRF) ;                    // number of faces of the region iR
    for ( int ilF=0; ilF<nRF; ++ilF ) {
      regn_flist.push_back(nFV) ;                  // number of vertices of the face ilF
      for ( int ip=0; ip<3; ++ip ) {
	regn_flist.push_back( vlist[ perm[ilF][ip] ] ) ;
      }
    }
  }
  regn_flist.push_back( nR ) ;
}
//--------------------------------------------------------------------------------------------
void ExtFileInput :: TETGEN_format( mesh_3Dv & mesh, string fname, int offset ) {
  MSG("start  ExtFileInput::TETGEN_format reading from disk file data"<<endl<<flush) ;
  vector<double> xV, yV, zV ;
  vector<int> fV, fR ;
  vector<int> regn_vlist, regn_flist ;
  vector<int> face_vlist, edge_vlist ;

  mesh3D_reader_TETGEN_format mesh_reader(fname,offset) ;
  mesh_reader.read_the_mesh( xV, yV, zV, fV, regn_vlist, fR ) ;
  //print_logfile( xV, yV, zV, fV, regn_vlist, fR, offset ) ;

  // build the list of region faces for tetrahedral meshes 
  // (to be used when data is input from TETGEN)
  set_regn_face_list( regn_flist, regn_vlist ) ;

  //  ---- start builder
  mesh3Dv_builder mesh_builder(mesh) ;
  mesh_builder . build_the_mesh( xV, yV, zV, fV, regn_flist, fR ) ; // !!!
  
  MSG("end of ExtFileInput::TETGEN_format reading from disk file data"<<endl<<flush) ;
}
//--------------------------------------------------------------------------------------------
void ExtFileInput :: REGN_FACE_format( mesh_3Dv & mesh, string fname, int offset ) {
  MSG("start  ExtFileInput::REGN_FACE_format reading from disk file data"<<endl<<flush) ;
  vector<double> xV, yV, zV ;
  vector<int> fV, fR ;
  vector<int> regn_flist ;

  mesh3D_reader_REGN_FACE_format mesh_reader(fname,offset) ;
  mesh_reader.read_the_mesh( xV, yV, zV, fV, regn_flist, fR ) ;
  //print_regn_face( xV, yV, zV, fV, regn_flist, fR, offset ) ;

  //  ---- start builder
  mesh3Dv_builder mesh_builder(mesh) ;
  mesh_builder . build_the_mesh( xV, yV, zV, fV, regn_flist, fR ) ; // !!!

  MSG("end of ExtFileInput::REGN_FACE_format reading from disk file data"<<endl<<flush) ;
}
//--------------------------------------------------------------------------------------------
void ExtFileInput :: MSH_format( mesh_3Dv & mesh, string fname, int offset ) {
  MSG("start  ExtFileInput::MSH_format reading from disk file data"<<endl<<flush) ;
  vector<double> xV, yV, zV ;
  vector<int> fV, fR ;
  vector<int> regn_flist ;

  mesh3D_reader_MSH_format mesh_reader(fname,offset) ;
  mesh_reader.read_the_mesh( xV, yV, zV, fV, regn_flist, fR ) ;
  //print_regn_face( xV, yV, zV, fV, regn_flist, fR, offset ) ;

  //  ---- start builder
  mesh3Dv_builder mesh_builder(mesh) ;
  mesh_builder . build_the_mesh( xV, yV, zV, fV, regn_flist, fR ) ; // !!!

  MSG("end of ExtFileInput::MSH_format reading from disk file data"<<endl<<flush) ;
}
//--------------------------------------------------------------------------------------------
void ExtFileInput :: multi_layer_format( mesh_3Dv & mesh, string fname, int n_layer, int offset, bool tilted_planes ) {
  MSG("begin  building layered mesh from 2D data"<<endl) ;

  vector<double> xV_2D, yV_2D  ;
  vector<int>    fV_2D, fR_2D  ;
  vector<int>    regn_vlist_2D ;

  mesh2D_reader_GeneralFormat mesh_reader(fname,offset) ;
  mesh_reader.read_mesh( xV_2D, yV_2D, fV_2D, regn_vlist_2D, fR_2D ) ;
  
  vector<double> xV, yV, zV ;
  vector<int>    regn_flist, fV, fR ;
  
  PRT(n_layer) ;
  PRT(offset) ;

  int    nL = n_layer ;
  double dz = 1./double(nL) ;

  // randomized planes
  vector<double> rnd_x(nL+1), rnd_y(nL+1) ;
  bool rnd_planar_coordinate_Z = tilted_planes ;
  if ( rnd_planar_coordinate_Z && false ) {
    for ( int jL=0; jL<=nL; ++jL ) {
      double fac = double(jL*(nL-jL)) / (double(nL*nL)/4.) ;

      double rnd_dz = 0.4*dz*(0.5-double(rand())/double(RAND_MAX)) ;
      assert(  -dz/2.<rnd_dz && rnd_dz<dz/2. ) ;

      rnd_x[jL] = jL==0||jL==nL ? 0. : rnd_dz ; //  0.45*fac*dz ;
      rnd_y[jL] = jL==0||jL==nL ? 0. : rnd_dz ; // -0.45*fac*dz ;
    }
  } else {
    for ( int jL=0; jL<=nL; ++jL ) {
      rnd_x[jL] = 0. ;
      rnd_y[jL] = 0. ;
    }
  }
  
  // VERTEX DATA
  int nV_2D = xV_2D.size() ;
  for ( int jL=0; jL<=nL; ++jL ) {
    for ( int iV=0; iV<nV_2D; ++iV ) {
      xV.push_back( xV_2D[iV] ) ;
      yV.push_back( yV_2D[iV] ) ;
      // ------------------------------------------------------------------------
      double z_vrtx = dz*double(jL) + rnd_x[jL]*xV_2D[iV] + rnd_y[jL]*yV_2D[iV] ;
      zV.push_back( z_vrtx ) ;
      //VAL(jL) ; VALA(iV,xV_2D) ; VALA(iV,yV_2D) ; VAL(dz*double(jL)) ; PRT(z_vrtx) ;
      // ------------------------------------------------------------------------
      fV.push_back( fV_2D[iV] ) ;
    } 
  }
  assert( xV.size() == nV_2D*(nL+1) ) ;

  // REGION DATA
  // loop on region layers
  int nR_2D = regn_vlist_2D.back() ;
  for ( int iL=0; iL<nL; ++iL ) {

    int k = 0 ;
    for ( int iR=0; iR<nR_2D; ++iR ) {      
      // get data of 2D region iR
      int nRV = regn_vlist_2D[k++] ;
      vector<int> tmp_vlist ;
      for ( int j=0; j<nRV; ++j ) { tmp_vlist.push_back( regn_vlist_2D[k++] ) ; } 
      
      // set fR for the region iR + iL*nR
      fR.push_back( iL ) ;
      
      // number of faces = nRV+2
      regn_flist.push_back( nRV+2 ) ;
      
      // number of vertices of bottom face
      regn_flist.push_back( nRV ) ;
      for ( int j=nRV-1; j>=0; --j ) {
	regn_flist.push_back( tmp_vlist[j]+nV_2D*iL   ) ;
      }

      // number of vertices of top face
      regn_flist.push_back( nRV ) ;
      for ( int j=0; j<nRV; ++j ) {
	regn_flist.push_back( tmp_vlist[j]+nV_2D*(iL+1) ) ;
      }
      
      // loop on vertical faces
      for ( int j=0; j<nRV; ++j ) {
	int j1 = (j+1)%nRV ;
	regn_flist.push_back( 4 ) ;
	regn_flist.push_back( tmp_vlist[j ]+nV_2D*iL     ) ;
	regn_flist.push_back( tmp_vlist[j1]+nV_2D*iL     ) ;
	regn_flist.push_back( tmp_vlist[j1]+nV_2D*(iL+1) ) ;
	regn_flist.push_back( tmp_vlist[j ]+nV_2D*(iL+1) ) ;
      }

    } // end loop on regions
  } // end loop on region's layers
  regn_flist.push_back( nR_2D*nL ) ;

  // output the mesh log file
  //print_regn_face( xV, yV, zV, fV, regn_flist, fR, offset ) ;

  //  ---- start builder
  mesh3Dv_builder mesh_builder(mesh) ;
  mesh_builder . build_the_mesh( xV, yV, zV, fV, regn_flist, fR ) ; // !!!
  
  MSG("end of building layered mesh from 2D data"<<endl) ;
}

#endif // end of _MESH_3D_GREADER_HH
