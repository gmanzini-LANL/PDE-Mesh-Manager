/* ----------------------------------------------------------------/
This is open source software; you can redistribute it and/or modify
it under the terms of the BSD-3 License. If software is modified
to produce derivative works, such modified software should be
clearly marked, so as not to confuse it with the version available
from LANL. Full text of the BSD-3 License can be found in the
License file of the repository.
/---------------------------------------------------------------- */

#ifndef _MESH3D_WRITER_HH
#define _MESH3D_WRITER_HH

#include <cassert>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
using namespace std ;

// include class declarations
#include "mesh3D.hh"

// add to MeshDefs
namespace aux_MeshDefs {
  static const int UNDEF    = -999 ;
  static const int BND_CELL = -998 ;
  static const int BND_FACE = -997 ;
  static const int BND_NODE = -996 ;
}
using namespace aux_MeshDefs ;

/*

  OUTPUT DATASETS IN A SINGLE FILE USING THE MSH FORMAT (F. HUBERT)
  (USED IN THE 3-D BENCHMARK OF FVCA-VI CONFERENCE)

  // offset is set to 1 (as default value)
  //
  // The file has a very strong structure determined by the initial header lines.
  //
  // The header may be in english (default) or french (as originally was) and the
  // switching between the languages is driven by lang_flag
  //
  // The (input/ouput) is unaffected by the language choice, but must respect the 
  // file structure)
  //
  // Due to this very different structure, this class has its own implementation and
  // can be used standalone. 
  //
  // The class is also used in mesh3Dv_writer::write_mesh_MSH(...) to write a mesh in
  // MSH format using an mesh3Dv_writer object.

 */
class mesh3Dv_writer_MSH {
protected:
  enum lang { FRENCH=0, ENGLISH=1 } ; 

protected:
  mesh_3Dv & mesh ;
  const int offset ;
  const int lang_flag ;

public:
  mesh3Dv_writer_MSH( mesh_3Dv & _mesh, int _offset=1 ) : 
    mesh(_mesh), offset(_offset), lang_flag(ENGLISH) {}
  ~mesh3Dv_writer_MSH() {}
  
  // write DATASETS
  void write_header  ( ostream & OUTF, string mesh_name_str ) ;
  void write_coords  ( ostream & OUTF ) ;
  void write_RegnFace( ostream & OUTF ) ;
  void write_RegnVrtx( ostream & OUTF ) ;
  void write_FaceEdge( ostream & OUTF ) ;
  void write_FaceVrtx( ostream & OUTF ) ;
  void write_FaceRegn( ostream & OUTF ) ;
  void write_EdgeVrtx( ostream & OUTF ) ;
  // ---
  void write_mesh_MSH( string fname=string("mesh3D") ) ;
} ;

// ----------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------
//
// general class to drive writer methods in different format
//
// ----------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------

class mesh3Dv_writer {
protected:
  static const int DIM = 3 ;
  mesh_3Dv & mesh ;
  int offset ;

  void write_node_coords( string fname=string("mesh3D") ) ;

  int n_vflag ; // number of external vertex flags (=0,1)
  int n_fflag ; // number of external face   flags (=0,1)
  int n_rflag ; // number of external region flags (=0,1)
  
public:
  mesh3Dv_writer( mesh_3Dv & _mesh, int _offset=0  ) : 
    mesh(_mesh), offset(_offset), n_vflag(0), n_fflag(0), n_rflag(0) {}
  ~mesh3Dv_writer() {}

  void write_mesh_RF ( string fname=string("mesh3D") ) ;
  void write_mesh_Fb ( string fname=string("mesh3D") ) ;
  void write_mesh_MSH( string fname=string("mesh3D") ) ;

  // for Sukumar output
  void write_node_coords_Suku_format( string fname ) ;
  void write_mesh_RF_Suku_format    ( string fname ) ;

  int  get_offset()              { return offset ; }
  void set_offset( int _offset ) { offset=_offset ; }
} ;

// ----------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------
//
// methods for writing a mesh file in MSH format
//
// ----------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------

void mesh3Dv_writer_MSH :: write_header( ostream & OUTF, string mesh_name_str ) {
  switch( lang_flag ) {
  case FRENCH:
    // french header
    OUTF << "Maillage cree par mesh3Dv_writer_MSH" << endl << endl << endl ;
    OUTF << "Version" << endl << "1" << endl ;
    OUTF << "Nom du maillage"        << endl << mesh_name_str << endl ;
    OUTF << "Infos sur le maillage"  << endl ;
    OUTF << "Nombre de sommets"      << endl << "  " << mesh.n_vertex() << endl ;
    OUTF << "Nombre de volumes"      << endl << "  " << mesh.n_region() << endl ;
    OUTF << "Nombre de faces"        << endl << "  " << mesh.n_face  () << endl ;
    OUTF << "Nombre d\'aretes"       << endl << "  " << mesh.n_edge  () << endl ;
    break;
  case ENGLISH:
    // english header
    OUTF << "Mesh created by mesh3Dv_writer_MSH" << endl << endl << endl ;
    OUTF << "Version"   << endl << "1" << endl ;
    OUTF << "Mesh name" << endl << mesh_name_str << endl ;
    OUTF << "Information on the mesh"   << endl ;
    OUTF << "Number of vertices"        << endl << "  " << mesh.n_vertex() << endl ;
    OUTF << "Number of control volumes" << endl << "  " << mesh.n_region() << endl ;
    OUTF << "Number of faces"           << endl << "  " << mesh.n_face  () << endl ;
    OUTF << "Nomber of edges"           << endl << "  " << mesh.n_edge  () << endl ;
    break;
  default:
    assert(false) ;
    break ;
  }
}
void mesh3Dv_writer_MSH :: write_coords( ostream & OUTF ) {
  // vertex coordinates
  switch( lang_flag ) {
  case FRENCH:  OUTF << "Sommets "  << mesh.n_vertex() << endl ; break ;
  case ENGLISH: OUTF << "Vertices " << mesh.n_vertex() << endl ; break ;
  }
  for ( int iV=0; iV<mesh.n_vertex(); ++iV ) {
    OUTF << setprecision(14) 
	 << scientific
	 << "   " << mesh.coords_V( iV, 0 ) 
	 << "   " << mesh.coords_V( iV, 1 ) 
	 << "   " << mesh.coords_V( iV, 2 ) << endl ;
  }
}
void mesh3Dv_writer_MSH :: write_RegnFace( ostream & OUTF ) {
  // RegnFace
  switch( lang_flag ) {
  case FRENCH:  OUTF << "Volumes->faces " << mesh.n_region() << endl ; break;
  case ENGLISH: OUTF << "Volumes->faces " << mesh.n_vertex() << endl ; break ;
  }
  for ( int iR=0; iR<mesh.n_region(); ++iR ) {
    vector<int> flist ;
    mesh.get_regn_face( iR, flist ) ;
    OUTF << "  " << flist.size() ;
    for ( int j=0; j<flist.size(); ++j ) {
      OUTF << " " << flist[j]+offset ;
    }
    OUTF << endl ;
  }
}
void mesh3Dv_writer_MSH :: write_RegnVrtx( ostream & OUTF ) {
  // RegnVrtx
  switch( lang_flag ) {
  case FRENCH:  OUTF << "Volumes->sommets "  << mesh.n_region() << endl ; break;
  case ENGLISH: OUTF << "Volumes->Vertices " << mesh.n_region() << endl ; break;
  }
  for ( int iR=0; iR<mesh.n_region(); ++iR ) {
    vector<int> vlist ;
    mesh.get_regn_vrtx( iR, vlist ) ;
    OUTF << "  " << vlist.size() ;
    for ( int j=0; j<vlist.size(); ++j ) {
      OUTF << " " << vlist[j]+offset ;
    }
    OUTF << endl ;
  }
}
void mesh3Dv_writer_MSH :: write_FaceEdge( ostream & OUTF ) {
  // FaceEdge
  switch( lang_flag ) {
  case FRENCH:  OUTF << "Faces->Aretes " << mesh.n_face() << endl ; break;
  case ENGLISH: OUTF << "Faces->Edges "  << mesh.n_face() << endl ; break;
  }
  for ( int iF=0; iF<mesh.n_face(); ++iF ) {
    vector<int> elist ;
    mesh.get_face_edge( iF, elist ) ;
    OUTF << "  " << elist.size() ;
    for ( int j=0; j<elist.size(); ++j ) {
      OUTF << " " << elist[j]+offset ;
    }
    OUTF << endl ;
  }
}
void mesh3Dv_writer_MSH :: write_FaceVrtx( ostream & OUTF ) {
  // FaceVrtx
  switch( lang_flag ) {
  case FRENCH:  OUTF << "Faces->Sommets "  << mesh.n_face() << endl ; break;
  case ENGLISH: OUTF << "Faces->Vertices " << mesh.n_face() << endl ; break;
  }
  for ( int iF=0; iF<mesh.n_face(); ++iF ) {
    vector<int> vlist ;
    mesh.get_face_vrtx( iF, vlist ) ;
    OUTF << "  " << vlist.size() ;
    for ( int j=0; j<vlist.size(); ++j ) {
      OUTF << " " << vlist[j]+offset ;
    }
    OUTF << endl ;
  }
}
void mesh3Dv_writer_MSH :: write_FaceRegn( ostream & OUTF ) {
  // FaceRegn
  switch( lang_flag ) {
  case FRENCH:  OUTF << "Faces->volumes "         << mesh.n_face() << endl ; break ;
  case ENGLISH: OUTF << "Faces->Control volumes " << mesh.n_face() << endl ; break ;
  }
  for ( int iF=0; iF<mesh.n_face(); ++iF ) {
    vector<int> rlist ;
    mesh.get_face_regn( iF, rlist ) ;
    if ( rlist.size()==1 ){
      OUTF << "  " << rlist[0]+offset << "  -1" << endl ;
    } else {
      OUTF << "  " << rlist[0]+offset 
	   << "  " << rlist[1]+offset << endl ; 
    } 
  }
}
// ------------------------------------------------------------------------------------------
void mesh3Dv_writer_MSH :: write_EdgeVrtx( ostream & OUTF ) {
  // EdgeVrtx
  switch( lang_flag ) {
  case FRENCH:  OUTF << "Aretes " << mesh.n_edge() << endl ; break;
  case ENGLISH: OUTF << "Edges "  << mesh.n_edge() << endl ; break;
  }
  for ( int iE=0; iE<mesh.n_edge(); ++iE ) {
    vector<int> vlist ;
    mesh.get_edge_vrtx( iE, vlist ) ;
    assert( vlist.size()==2 ) ;
    for ( int j=0; j<vlist.size(); ++j ) {
      OUTF << " " << vlist[j]+offset ;
    }
    OUTF << endl ;
  }
}
// --------------------------------------------------------------------------------------------
void mesh3Dv_writer_MSH :: write_mesh_MSH( string fname ) {
  fname += string(".msh") ;
  ofstream OUTF(fname.c_str()) ;
  write_header  ( OUTF, fname ) ;
  write_coords  ( OUTF ) ;
  write_RegnFace( OUTF ) ;
  write_RegnVrtx( OUTF ) ;
  write_FaceEdge( OUTF ) ;
  write_FaceVrtx( OUTF ) ;
  write_FaceRegn( OUTF ) ;
  write_EdgeVrtx( OUTF ) ;
  OUTF.close() ;
}

// ----------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------
//
// methods of class mesh3Dv_writer
//
// ----------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

void mesh3Dv_writer::write_node_coords( string fname ) {
  // get number of mesh vertices
  int nV = mesh.n_vertex() ;
  
  // open output file for node's coordinates
  string fnode = fname + string(".node") ;
  ofstream out_node(fnode.c_str()) ;
  
  // output the node file's header
  out_node << "# *.node file of 3D ATS mesh in REGN_FACE format " << endl ;
  out_node << "# offset = " << offset << endl ;
  out_node << "# fname  = " << fnode  << endl ;
  out_node << nV << "  3  0  " << n_vflag << endl ;

  // output the node's coordinates
  for ( int iV=0; iV<nV; ++iV ) {
    int vrtx_id = iV + offset ;
    out_node << setw(8)  << vrtx_id  << "     " 
	     << setw(20) << setprecision(16) << mesh.coords_V(iV,0) << "   "  
	     << setw(20) << setprecision(16) << mesh.coords_V(iV,1) << "   "  
	     << setw(20) << setprecision(16) << mesh.coords_V(iV,2) << "   "  ;
    if ( n_vflag==1 ) { out_node << " 1" ; }
    out_node << endl ;
  }
  out_node << "# output from mesh3Dv_writer.hh " << endl ;
  out_node.close() ;
}

void mesh3Dv_writer :: write_mesh_Fb( string fname ) {
  cout<<"start  mesh3Dv_writer_FACE_BASED_format::write_mesh"<<endl<<flush ;
  cout<<"write \""<<fname<<"\""<<endl<<flush ;
  
  int nV = mesh.n_vertex() ;
  int nF = mesh.n_face() ;
  int nR = mesh.n_region() ;

  int wV = int( log(nV) + 2 ) ; 
  int wF = int( log(nF) + 2 ) ;
  int wR = int( log(nR) + 2 ) ; 

  // write node coordinates
  write_node_coords( fname ) ;

  // open output file for face structure
  string fface = fname + string(".face") ;
  ofstream out_face(fface.c_str()) ;

  // output the faces file's header
  out_face << "# *.face file of 3D-mesh in FACE_based (Fb) format " << endl ;
  out_face << "# offset = " << offset << endl ;
  out_face << "# fname  = " << fface  << endl ;
  out_face << nF << "  " << n_fflag << endl ;

  for ( int iF=0; iF<nF; ++iF ) {
    
    // get node and cell lists
    int face_id = iF + offset ;
    vector<int> vlist, rlist ;
    mesh.get_face_vrtx( iF, vlist ) ;
    mesh.get_face_regn( iF, rlist ) ;

    // get size (and check)
    int nFV = vlist.size() ;
    int nFR = rlist.size() ;
    assert( nFV>3 ) ;
    assert( nFR==1 || nFR==2 ) ;
    
    // list of face's nodes
    out_face << setw(wV) << face_id << "   "
	     << setw(wV) << nFV     << "   " ;
    for ( int ilV=0; ilV<nFV; ++ilV ) {
      out_face << setw(wV) << vlist[ilV] + offset << " " ;
    }
    
    // connected cells
    out_face << "   " ;
    if ( nFR==2 ) {
      out_face << setw(wR) << rlist[0] + offset << "   " ;
      out_face << setw(wR) << rlist[1] + offset ;
    } else {
      out_face << setw(wR) << rlist[0] + offset << "   " ;
      out_face << setw(wR) << BND_FACE ;
    }

    // external flag
    if ( n_fflag==1 ) { out_face << " 1" ; }

    // end of face record
    out_face << endl ;
  }
  
  // close face file
  out_face.close() ;
  cout<<"end of mesh3Dv_writer_FACE_BASED_format::write_mesh"<<endl<<flush ;
}

void mesh3Dv_writer :: write_mesh_RF( string fname ) {
  cout <<"start  mesh3Dv_writer_REGN_FACE_format::write_mesh"<<endl<<flush ;
  cout <<"write \""<<fname<<"\""<<endl<<flush ;

  // write node coordinates
  write_node_coords( fname ) ;

  // set number of mesh cells
  int nR = mesh.n_region() ;

  // open output file for element's structure
  string fele  = fname + string(".ele") ;
  ofstream out_ele(fele.c_str()) ;

  // output the element file's header
  out_ele << "# *.ele file of 3D-mesh in REGN_FACE format " << endl ;
  out_ele << "# offset = " << offset << endl ;
  out_ele << "# fname  = " << fele << endl ;
  out_ele << nR << "  " << n_rflag << endl ;

  // output the element's structure
  for ( int iR=0 ; iR<mesh.n_region() ; ++iR ) {
    
    // get cell ID
    int cell_id = iR + offset ; 

    // get the list of face IDs
    vector<int> flist ;
    mesh.get_regn_face( iR, flist ) ;
    
    // get the list of face dirs (NOT USED)
    //vector<int> fdirs(flist.size()) ;
    //for ( int i=0; i<flist.size(); ++i ) { fdirs[i] = mesh.ok_regn_face(iR,i)? 1 : -1 ; }
    
    // set the number of the local faces
    int nRF = mesh.n_regn_face(iR) ;
    assert( nRF==flist.size() ) ;

    // write region header
    out_ele << cell_id << "  "  << nRF ;
    if ( n_rflag==1 ) { out_ele << mesh.get_fR(iR) ; }
    out_ele << endl ;

    // write the face's node
    for ( int ilF=0; ilF<nRF; ++ilF ) {

      // get the face ID
      int iF  = flist[ilF] ;

      // get face's node list
      vector<int> vlist ;
      mesh.get_face_vrtx( iF, vlist ) ;
      
      // get the face size
      int nFV = vlist.size() ;
      out_ele << "  " << ilF+offset << "  " << nFV << "  " ;
      
      // write the list of nodes in counterclockwise fashion
      if ( mesh.ok_regn_face(iR,ilF) ) { 
	for ( int ilV=0; ilV<nFV; ++ilV ) {
	  int vrtx_ID = vlist[ilV]+offset ;
	  out_ele << "  " << vrtx_ID ;
	}
      } else {
	for ( int ilV=0; ilV<nFV; ++ilV ) {
	  int vrtx_ID = vlist[nFV-1-ilV]+offset ;
	  out_ele << "  " << vrtx_ID ;
	}
      }
      // end of element record
      out_ele << endl ;
    }
  }
  
  out_ele << "# output from mesh3Dv_writer.hh " << endl ;
  out_ele.close() ;
  
  cout <<"end of mesh3Dv_writer_REGN_FACE_format::write_mesh"<<endl<<flush ;
}

void mesh3Dv_writer :: write_mesh_MSH( string fname ) {
  cout <<"start  mesh3Dv_writer_MSH_format::write_mesh"<<endl<<flush ;
  cout <<"write \""<<fname<<"\""<<endl<<flush ;
  
  const int fortran_offset = 1 ;
  mesh3Dv_writer_MSH mesh_writer(mesh,fortran_offset) ;
  mesh_writer.write_mesh_MSH( fname ) ;

  cout <<"end of mesh3Dv_writer_MSH_format::write_mesh"<<endl<<flush ;
}

// ---------------------------------------------------------------------

void mesh3Dv_writer::write_node_coords_Suku_format( string fname ) {
  // get number of mesh vertices
  int nV = mesh.n_vertex() ;
  
  // open output file for node's coordinates
  string fnode = fname + string(".node") ;
  ofstream out_node(fnode.c_str()) ;
  
  // output the node's coordinates
  for ( int iV=0; iV<nV; ++iV ) {
    int bnd_flag = mesh.is_boundary_vrtx( iV ) ? 1 : 0 ;
    int vrtx_id = iV + offset ;
    out_node << setw(20) << setprecision(16) << mesh.coords_V(iV,0) << "   "  
	     << setw(20) << setprecision(16) << mesh.coords_V(iV,1) << "   "  
	     << setw(20) << setprecision(16) << mesh.coords_V(iV,2) << "   " 
	     << setw(8)  << bnd_flag ;
    out_node << endl ;
  }
  out_node.close() ;
}

void mesh3Dv_writer :: write_mesh_RF_Suku_format( string fname ) {
  cout <<"start  mesh3Dv_writer_REGN_FACE_format::write_mesh"<<endl<<flush ;
  cout <<"write \""<<fname<<"\""<<endl<<flush ;

  // write node coordinates
  write_node_coords_Suku_format( fname ) ;

  // set number of mesh cells
  int nR = mesh.n_region() ;

  // open output file for element's structure
  string fele  = fname + string(".ele") ;
  ofstream out_ele(fele.c_str()) ;

  // output the element's structure
  for ( int iR=0 ; iR<mesh.n_region() ; ++iR ) {
    
    // get cell ID
    int cell_id = iR + offset ; 

    // get the list of face IDs
    vector<int> flist ;
    mesh.get_regn_face( iR, flist ) ;
    
    // set the number of the local faces
    int nRF = mesh.n_regn_face(iR) ;
    assert( nRF==flist.size() ) ;

    // write number of faces
    out_ele << nRF << "  " ;

    // write the face's node
    for ( int ilF=0; ilF<nRF; ++ilF ) {

      // get the face ID
      int iF  = flist[ilF] ;

      // get face's node list
      vector<int> vlist ;
      mesh.get_face_vrtx( iF, vlist ) ;
      
      // get the face size
      int nFV = vlist.size() ;
      out_ele << "  " << nFV ;
      
      // write the list of nodes in counterclockwise fashion
      if ( mesh.ok_regn_face(iR,ilF) ) { 
	for ( int ilV=0; ilV<nFV; ++ilV ) {
	  int vrtx_ID = vlist[ilV]+offset ;
	  out_ele << "  " << vrtx_ID ;
	}
      } else {
	for ( int ilV=0; ilV<nFV; ++ilV ) {
	  int vrtx_ID = vlist[nFV-1-ilV]+offset ;
	  out_ele << "  " << vrtx_ID ;
	}
      }
      // end of element record
      out_ele << "  " ;
    }
    out_ele << endl ;
  }
  out_ele.close() ;
  
  cout <<"end of mesh3Dv_writer_REGN_FACE_format::write_mesh"<<endl<<flush ;
}

#endif // end of --> _MESH3D_WRITER_CC
