/* ----------------------------------------------------------------/
This is open source software; you can redistribute it and/or modify
it under the terms of the BSD-3 License. If software is modified
to produce derivative works, such modified software should be
clearly marked, so as not to confuse it with the version available
from LANL. Full text of the BSD-3 License can be found in the
License file of the repository.
/---------------------------------------------------------------- */

#ifndef _MESH3D_MASTER_HH
#define _MESH3D_MASTER_HH

#include <algorithm>
#include <vector>
#include <iostream>
using namespace std ;

#include "mesh3D.hh"

class mesh3Dv_cmaster {

private:
  static const int DIM = 3 ;
  const double TOL ;

private:
  mesh_3Dv & mesh  ;
  bool mesh_status ;
  bool regn_status ;
  bool face_status ;
  bool edge_status ;
  bool vrtx_status ;
  bool ornt_status ;

public:
  mesh3Dv_cmaster( mesh_3Dv & _mesh ) : 
    TOL(1.e-14),
    mesh(_mesh), 
    mesh_status(true),
    regn_status(true), 
    face_status(true), 
    edge_status(true), 
    vrtx_status(true),
    ornt_status(true) {}

  // check the mesh
  void check_the_mesh() ;

  // check the regions
  bool check_the_regions() ;
  bool check_regn_regn(int i) ;
  bool check_regn_face(int i) ;
  bool check_regn_edge(int i) ;
  bool check_regn_vrtx(int i) ;

  // check the faces
  bool check_the_faces() ;
  bool check_face_regn(int i) ;
  bool check_face_face(int i) ;
  bool check_face_edge(int i) ;
  bool check_face_vrtx(int i) ;

  // check the edges
  bool check_the_edges() ;
  bool check_edge_regn(int i) ;
  bool check_edge_face(int i) ;
  bool check_edge_edge(int i) ;
  bool check_edge_vrtx(int i) ;

  // check the vertices
  bool check_the_vertices() ;
  bool check_vrtx_regn(int i) ;
  bool check_vrtx_face(int i) ;
  bool check_vrtx_edge(int i) ;
  bool check_vrtx_vrtx(int i) ;
  // ------------------
  //void check_this_item( int iE ) ;
  // ------------------
  bool check_local_regions () ;
  bool check_local_faces   () ;
  bool check_local_edges   () ; 
  bool check_local_vertices() ; 

  // check orientation and local numberings
  bool check_regn_face_orientation() ;
  bool check_regn_face_numbering  () ;
  bool check_face_edge_numbering  () ;
  bool check_vrtx_vrtx_numbering  () ;
  bool check_edge_regn_numbering  () ;

  // check some geometric quantities
  bool check_geometric_region_factor() ;
  bool check_region_face_orientation() ;
} ;

void mesh3Dv_cmaster::check_the_mesh() {
  cout<<"-------------------------------------------------------"<<endl<<flush ;
  cout<<"-------------------------------------------------------"<<endl<<flush ;
  cout<<"BEGIN CHECKING THE MESH"<<endl<<flush ;
  if ( regn_status ) { mesh_status &= check_the_regions () ; }
  if ( face_status ) { mesh_status &= check_the_faces   () ; }
  if ( edge_status ) { mesh_status &= check_the_edges   () ; }
  if ( vrtx_status ) { mesh_status &= check_the_vertices() ; }
  if ( mesh_status ) {
    cout<<"Mesh datasets are ok!"<<endl<<endl<<flush ;
  } else {
    cout<<"-->>>THERE ARE TROUBLES IN MESH DATASETS!!!"<<endl<<endl<<flush ;
  }
  if ( regn_status ) { mesh_status &= check_local_regions () ; }
  if ( face_status ) { mesh_status &= check_local_faces   () ; }
  if ( edge_status ) { mesh_status &= check_local_edges   () ; }
  if ( vrtx_status ) { mesh_status &= check_local_vertices() ; }
  if ( mesh_status ) {
    cout<<"Local mesh items are ok!"<<endl<<endl<<flush ;
  } else {
    cout<<"-->>>THERE ARE TROUBLES IN LOCAL MESH ITEMS!!!"<<endl<<endl<<flush ;
  } 
  if ( ornt_status ) { 
    mesh_status &= 
      check_regn_face_orientation() && 
      check_regn_face_numbering  () && 
      check_face_edge_numbering  () &&
      check_vrtx_vrtx_numbering  () ; 
  }
  cout<<"-------------------------------------------------------"<<endl<<flush ;
  cout<<"-------------------------------------------------------"<<endl<<flush ;
  if ( !check_geometric_region_factor() ) {
    cout<<"-->>>THERE ARE TROUBLES IN REGION GEOMETRIC VALUES!!!"<<endl<<endl<<flush ;
  }
  if ( !check_region_face_orientation() ) {
    cout<<"-->>>THERE ARE TROUBLES IN FACE GEOMETRIC VALUES!!!"<<endl<<endl<<flush ;
  }
}

#define check(F0,F1)							\
  bool mesh3Dv_cmaster :: check##F0##F1( int i ) {			\
    bool retval(true) ;							\
    vector<int> klist ;							\
    mesh.get##F0##F1( i, klist ) ;					\
    for ( int il=0; il<klist.size(); ++il ) {				\
      int k = klist[il] ;						\
      vector<int> plist ;						\
      mesh.get##F1##F0( k, plist ) ;					\
      sort( plist.begin(), plist.end() ) ;				\
      retval &= binary_search( plist.begin(), plist.end(), i ) ;	\
    }									\
    return retval ;							\
  }									

#define check_item(I,X)						\
  if ( !X ) {							\
    cout << endl << "-->" << (#I) << " = "    << I		\
	 << ", failed : " << (#X) << "  !!! " << endl ;		\
    local_status = false ;					\
  }

// check the regions
check(_regn,_regn) ;
check(_regn,_face) ;
check(_regn,_edge) ;
check(_regn,_vrtx) ;
// check the faces
check(_face,_regn) ;
check(_face,_face) ;
check(_face,_edge) ;
check(_face,_vrtx) ;
// check the edges
check(_edge,_regn) ;
check(_edge,_face) ;
check(_edge,_edge) ;
check(_edge,_vrtx) ;
// check the vertices
check(_vrtx,_regn) ;
check(_vrtx,_face) ;
check(_vrtx,_edge) ;
check(_vrtx,_vrtx) ;

bool mesh3Dv_cmaster :: check_geometric_region_factor() {
  cout<<"-->check geom regions  "<<endl<<flush ;
  bool local_status(true) ;
  int nR = mesh.n_region() ;
  double meas(0.), coords[DIM] ;
  for ( int s=0; s<DIM; ++s ) { coords[s]=0. ; }
  for ( int iR=0; iR<nR; ++iR ) {
    double meas_R = mesh.get_regn_measure(iR) ;
    meas += meas_R ;
    for ( int s=0; s<DIM; ++s ) {
      coords[s]  += meas_R * mesh.coords_R(iR,s) ;
    }
  }
  for ( int s=0; s<DIM; ++s ) { coords[s]/=meas ; }
  double coords_domain[DIM] ;
  for ( int s=0; s<DIM; ++s ) { 
    coords_domain[s] = ( mesh.min_coords(s)+mesh.max_coords(s) )/2. ; 
  }
  double measure_domain = 1. ;
  for ( int s=0; s<DIM; ++s ) { 
    measure_domain *= ( mesh.max_coords(s)-mesh.min_coords(s) ) ;
  }
  local_status = abs(meas-measure_domain)<TOL ;
  for ( int s=0; s<DIM; ++s ) { 
    local_status &= abs( coords[s]-coords_domain[s] )<TOL ;
  }

  double vol      = meas ;
  double vol_cube = measure_domain ;
  double xc       = coords[0] ;
  double yc       = coords[1] ;
  double zc       = coords[2] ;
  double xc_cube  = coords_domain[0] ;
  double yc_cube  = coords_domain[1] ;
  double zc_cube  = coords_domain[2] ;

  if ( local_status ) { 

    cout<<"--> OK (cube detected)"<<endl<<flush ;

  } else {
    
    cout << "vol = "      << vol << " "
	 << "vol_cube = " << vol_cube << endl ;
    
    cout << "xc = " << xc << " "
	 << "yc = " << yc << " "
	 << "zc = " << zc << endl ;
    
    cout << "xc_cube = " << xc_cube << " "
	 << "yc_cube = " << yc_cube << " "
	 << "zc_cube = " << zc_cube << endl ;

    cout << "TOL = " << TOL << endl ;
    cout << "vol-vol_cube = " << vol-vol_cube << endl ;
    cout << "xc-xc_cube   = " << xc-xc_cube   << endl ;
    cout << "yc-yc_cube   = " << yc-yc_cube   << endl ;
    cout << "zc-zc_cube   = " << zc-zc_cube   << endl ;

    cout<<endl<<">>>If your domain is a cube<<<<"<<endl<<flush ;
    cout<<endl<<">>>THE GEOMETRIC FACTORS OF MESH REGIONS ARE NOT CONSISTENT!<<<"<<endl<<endl<<flush ;
  }
  
  cout << "vol = "      << vol << " "
       << "vol_cube = " << vol_cube << endl ;
  
  cout << "xc = " << xc << " "
       << "yc = " << yc << " "
       << "zc = " << zc << endl ;
  
  cout << "xc_cube = " << xc_cube << " "
       << "yc_cube = " << yc_cube << " "
       << "zc_cube = " << zc_cube << endl ;
  
  return local_status ;
}

// face orientation
bool mesh3Dv_cmaster :: check_region_face_orientation() {
  cout<<"-->check geom faces    "<<endl<<flush ;
  bool local_status(true) ;
  int nR = mesh.n_region() ;
  //double vol(0.), xc(0.), yc(0.), zc(0.) ;
  for ( int iR=0; iR<nR; ++iR ) {
    int nRF = mesh.n_regn_face( iR ) ;
    for ( int ilF=0; ilF<nRF; ++ilF ) {
      int iF = mesh.regn_face(iR,ilF) ;
      double vFR[3] = { mesh.coords_F(iF,0)-mesh.coords_R(iR,0), 
			mesh.coords_F(iF,1)-mesh.coords_R(iR,1), 
			mesh.coords_F(iF,2)-mesh.coords_R(iR,2) } ;
      double nmF[3] = { mesh.get_nor(iF,0), mesh.get_nor(iF,1), mesh.get_nor(iF,2) } ;
      double vscal  = vFR[0]*nmF[0] + vFR[1]*nmF[1] + vFR[2]*nmF[2] ;
      double sgnRF  = mesh.ok_regn_face(iR,ilF) ? +1. : -1. ;
      local_status &= sgnRF*vscal>0 ;
      if ( sgnRF*vscal<0 ) {
	for ( int i=0; i<20; ++i ) { cout << "--" ; } 
	cout << endl ;
	cout<<"ERROR in check_region_face_orientation\n"<<endl<<flush ;
	cout << "iR = " << iR << " " << "ilF = " << ilF << " "
	     << "iF = " << iF << endl ;
	cout << "sgnRF = " << sgnRF << " "
	     << "vscal = " << vscal << endl ;
      } 
    }
  }
  if ( local_status ) { 
    cout<<"--> OK"<<endl<<endl<<flush ;
  } else {
    cout<<endl<<">>>ORIENTATION OF MESH FACES IN REGIONS IS NOT CONSISTENT!<<<"<<endl<<endl<<flush ;
  }
  return local_status ;
}

bool mesh3Dv_cmaster :: check_the_regions() {
  cout<<"-->check regions  "<<endl<<flush ;
  bool local_status(true) ;
  int nR = mesh.n_region() ;
  for ( int iR=0; iR<nR; ++iR ) {
    check_item( iR, check_regn_regn(iR) ) ;
    check_item( iR, check_regn_face(iR) ) ;
    check_item( iR, check_regn_edge(iR) ) ;
    check_item( iR, check_regn_vrtx(iR) ) ;
  }
  if ( local_status ) { 
    cout<<"--> OK"<<endl<<endl<<flush ;
  } else {
    cout<<endl<<">>>MESH REGIONS ARE NOT CONSISTENT!<<<"<<endl<<endl<<flush ;
  }
  return local_status ;
}

bool mesh3Dv_cmaster :: check_the_faces() {
  cout<<"-->check faces    "<<endl<<flush ;
  bool local_status(true) ;
  int nF = mesh.n_face() ;
  for ( int iF=0; iF<nF; ++iF ) {
    check_item( iF, check_face_regn(iF) ) ;
    check_item( iF, check_face_face(iF) ) ;
    check_item( iF, check_face_edge(iF) ) ;
    check_item( iF, check_face_vrtx(iF) ) ;
  }
  if ( local_status ) { 
    cout<<"--> OK"<<endl<<endl<<flush ;
  } else {
    cout<<endl<<">>>MESH FACES ARE NOT CONSISTENT!<<<"<<endl<<endl<<flush ;
  }
  return local_status ;
}

bool mesh3Dv_cmaster :: check_the_edges() {
  cout<<"-->check edges    "<<endl<<flush ;
  bool local_status(true) ;
  int nE = mesh.n_edge() ;
  for ( int iE=0; iE<nE; ++iE ) {
    check_item( iE, check_edge_regn(iE) ) ;
    check_item( iE, check_edge_face(iE) ) ;
    check_item( iE, check_edge_edge(iE) ) ;
    check_item( iE, check_edge_vrtx(iE) ) ;
  }
  if ( local_status ) { 
    cout<<"--> OK"<<endl<<endl<<flush ;
  } else {
    cout<<endl<<">>>MESH EDGES ARE NOT CONSISTENT!<<"<<endl<<flush ;
  }
  return local_status ;
}

bool mesh3Dv_cmaster :: check_the_vertices() {
  cout<<"-->check vertices "<<endl<<flush ;
  bool local_status(true) ;
  int nV = mesh.n_vertex() ;
  for ( int iV=0; iV<nV; ++iV ) {
    check_item( iV, check_vrtx_regn(iV) ) ;
    check_item( iV, check_vrtx_face(iV) ) ;
    check_item( iV, check_vrtx_edge(iV) ) ;
    check_item( iV, check_vrtx_vrtx(iV) ) ;
  }
  if ( local_status ) { 
    cout<<"--> OK"<<endl<<endl<<flush ;
  } else {
    cout<<endl<<">>>MESH VERTICES ARE NOT CONSISTENT!<<<"<<endl<<endl<<flush ;
  }
  return local_status ;
}

bool mesh3Dv_cmaster :: check_local_regions() {
  cout<<"-->check local regions  "<<endl<<flush ;
  bool local_status(true) ;
  vector<int> plist ;
  int nR = mesh.n_region() ;

  // check region consistency: 
  // get the list of connected/adjacent items of the region iR
  // for each item get the list of the connected regions and check
  // whether iR belongs to such lists
  for ( int iR=0; iR<nR; ++iR ) {
    // check regions: get_regn_regn, get_regn_regn
    vector<int> regn_rlist ;
    mesh.get_regn_regn(iR,regn_rlist) ;
    int nRR = regn_rlist.size() ;
    // check regions
    for ( int jlR=0; jlR<nRR; ++jlR ) {
      int jR = regn_rlist[jlR] ;
      plist.resize(0) ;
      mesh.get_regn_regn( jR, plist ) ; 
      sort( plist.begin(), plist.end() ) ;
      if ( !binary_search( plist.begin(), plist.end(), iR ) ) {
	local_status = bool(false) ;
	cout << "-->>region: " << iR << " region: " << jR << " failed!!!" << endl ; 
      }
    }
    // check faces: get_regn_face, get_face_regn
    vector<int> regn_flist ;
    mesh.get_regn_face(iR,regn_flist) ;
    int nRF = regn_flist.size() ;
    for ( int ilF=0; ilF<nRF; ++ilF ) {
      int iF = regn_flist[ilF] ;
      plist.resize(0) ;
      mesh.get_face_regn(iF,plist) ;
      sort( plist.begin(), plist.end() ) ;
      if ( !binary_search( plist.begin(), plist.end(), iR ) ) {
	local_status = bool(false) ;
	cout << "-->>region: " << iR << " face: " << iF << " failed!!!" << endl ; 
      }
    }
    // check edges: get_regn_edge, get_edge_regn
    vector<int> regn_elist ;
    mesh.get_regn_edge(iR,regn_elist) ;
    int nRE = regn_elist.size() ;
    for ( int ilE=0; ilE<nRE; ++ilE ) {
      int iE = regn_elist[ilE] ;
      plist.resize(0) ;
      mesh.get_edge_regn( iE, plist ) ;
      sort( plist.begin(), plist.end() ) ;
      if ( !binary_search( plist.begin(), plist.end(), iR ) ) {
	local_status = bool(false) ;
	cout << "-->>region: " << iR << " edge: " << iE << " failed!!!" << endl ; 
      }
    }
    // check vertices: get_regn_vrtx, get_vrtx_regn
    vector<int> regn_vlist ;
    mesh.get_regn_vrtx( iR, regn_vlist ) ;
    int nRV = regn_vlist.size() ;
    for ( int ilV=0; ilV<nRV; ++ilV ) {
      int iV = regn_vlist[ilV] ;
      plist.resize(0) ;
      mesh.get_vrtx_regn(iV,plist);
      sort( plist.begin(), plist.end() ) ;
      if ( !binary_search( plist.begin(), plist.end(), iR ) ) {
	local_status = bool(false) ;
	cout << "-->>region: " << iR << " vrtx: " << iV << " failed!!!" << endl ; 
      }
    }
  }
  if ( local_status ) { 
    cout<<"--> OK"<<endl<<endl<<flush ;
  } else {
    cout<<endl<<">>>MESH FACES ARE NOT LOCALLY CONSISTENT!<<<"<<endl<<endl<<flush ;
  }
  return local_status ;
}

bool mesh3Dv_cmaster :: check_local_faces() {
  cout<<"-->check local faces    "<<endl<<flush ;
  bool local_status(true) ;
  vector<int> plist ;
  int nF = mesh.n_face() ;
  for ( int iF=0; iF<nF; ++iF ) {
    // check regions
    vector<int> face_rlist ;
    mesh.get_face_regn( iF, face_rlist ) ;
    int nFR = face_rlist.size() ;
    for ( int ilR=0; ilR<nFR; ++ilR ) {
      int iR = face_rlist[ilR] ;
      plist.resize(0) ;
      mesh.get_regn_face(iR,plist) ;
      sort( plist.begin(), plist.end() ) ;
      if ( !binary_search( plist.begin(), plist.end(), iF ) ) {
	local_status = bool(false) ;
	cout << "-->>face: " << iF << " region: " << iR << " failed!!!" << endl ; 
      }
    }
    // check faces
    vector<int> face_flist ;
    mesh.get_face_face(iF,face_flist) ;
    int nFF = face_flist.size() ;
    for ( int jlF=0; jlF<nFF; ++jlF ) {
      int jF = face_flist[jlF] ;
      plist.resize(0) ;
      mesh.get_face_face( jF, plist ) ;
      sort( plist.begin(), plist.end() ) ;
      if ( !binary_search( plist.begin(), plist.end(), iF ) ) {
	local_status = bool(false) ;
	cout << "-->>face: " << iF << " face: " << jF << " failed!!!" << endl ; 
      }
    }
    // check edges
    vector<int> face_elist ;
    mesh.get_face_edge( iF, face_elist ) ;
    int nFE = face_elist.size() ;
    for ( int ilE=0; ilE<nFE; ++ilE ) {
      int iE = face_elist[ilE] ;
      plist.resize(0) ;
      mesh.get_edge_face( iE, plist ) ;
      sort( plist.begin(), plist.end() ) ;
      if ( !binary_search( plist.begin(), plist.end(), iF ) ) {
	local_status = bool(false) ;
	cout << "-->>face: " << iF << " edge: " << iE << " failed!!!" << endl ; 
      }
    }
    // check vertices
    vector<int> face_vlist ;
    mesh.get_face_vrtx( iF, face_vlist ) ;
    int nFV = face_vlist.size() ;
    for ( int ilV=0; ilV<nFV; ++ilV ) {
      int iV = face_vlist[ilV] ;
      plist.resize(0) ;
      mesh.get_vrtx_face(iV,plist) ;
      sort( plist.begin(), plist.end() ) ;
      if ( !binary_search( plist.begin(), plist.end(), iF ) ) {
	local_status = bool(false) ;
	cout << "-->>face: " << iF << " vrtx: " << iV << " failed!!!" << endl ; 
      }
    }
  }
  if ( local_status ) { 
    cout<<"--> OK"<<endl<<endl<<flush ;
  } else {
    cout<<endl<<">>>MESH FACES ARE NOT LOCALLY CONSISTENT!<<<"<<endl<<endl<<flush ;
  }
  return local_status ;
}

bool mesh3Dv_cmaster :: check_local_edges() {
  cout<<"-->check local edges    "<<endl<<flush ;
  bool local_status(true) ;
  vector<int> plist ;
  int nE = mesh.n_edge() ;
  for ( int iE=0; iE<nE; ++iE ) {
    // check regions
    vector<int> edge_rlist ;
    mesh.get_edge_regn(iE,edge_rlist) ;
    int nER = edge_rlist.size() ;
    for ( int ilR=0; ilR<nER; ++ilR ) {
      int iR = edge_rlist[ilR] ;
      plist.resize(0) ;
      mesh.get_regn_edge(iR,plist) ;
      sort( plist.begin(), plist.end() ) ;
      if ( !binary_search( plist.begin(), plist.end(), iE ) ) {
	local_status = bool(false) ;
	cout << "-->>edge: " << iE << " region: " << iR << " failed!!!" << endl ; 
      }
    }
    // check faces
    vector<int> edge_flist ;
    mesh.get_edge_face(iE,edge_flist) ;
    int nEF = edge_flist.size() ;
    for ( int ilF=0; ilF<nEF; ++ilF ) {
      int iF = edge_flist[ilF] ;
      plist.resize(0) ;
      mesh.get_face_edge(iF,plist) ;
      sort( plist.begin(), plist.end() ) ;
      if ( !binary_search( plist.begin(), plist.end(), iE ) ) {
	local_status = bool(false) ;
	cout << "-->>edge: " << iE << " face: " << iF << " failed!!!" << endl ; 
      }
    }
    // check edges
    vector<int> edge_elist ;
    mesh.get_edge_edge(iE,edge_elist) ;
    int nEE = edge_elist.size() ;
    for ( int jlE=0; jlE<nEE; ++jlE ) {
      int jE = edge_elist[jlE] ;
      plist.resize(0) ;
      mesh.get_edge_edge(jE,plist) ;
      sort( plist.begin(), plist.end() ) ;
      if ( !binary_search( plist.begin(), plist.end(), iE ) ) {
	local_status = bool(false) ;
	cout << "-->>edge: " << iE << " edge: " << jE << " failed!!!" << endl ; 
      }
    }
    // check vertices
    vector<int> edge_vlist ;
    mesh.get_edge_vrtx(iE,edge_vlist) ;
    int nEV = edge_vlist.size() ;
    for ( int ilV=0; ilV<nEV; ++ilV ) {
      int iV = edge_vlist[ilV] ;
      plist.resize(0) ;
      mesh.get_vrtx_edge(iV,plist) ;
      sort( plist.begin(), plist.end() ) ;
      if ( !binary_search( plist.begin(), plist.end(), iE ) ) {
	local_status = bool(false) ;
	cout << "-->>edge: " << iE << " vrtx: " << iV << " failed!!!" << endl ; 
      }
    }
  }
  if ( local_status ) { 
    cout<<"--> OK"<<endl<<endl<<flush ;
  } else {
    cout<<endl<<">>>MESH EDGES ARE NOT LOCALLY CONSISTENT!<<<"<<endl<<endl<<flush ;
  }
  return local_status ;
}

bool mesh3Dv_cmaster :: check_local_vertices() {
  cout<<"-->check local vertices "<<endl<<flush ;
  bool local_status(true) ;
  vector<int> plist ;
  int nV = mesh.n_vertex() ;
  for ( int iV=0; iV<nV; ++iV ) {
    // check regions
    vector<int> vrtx_rlist ;
    mesh.get_vrtx_regn(iV,vrtx_rlist) ;
    int nVR = vrtx_rlist.size() ;
    for ( int ilR=0; ilR<nVR; ++ilR ) {
      int iR = vrtx_rlist[ilR] ;
      plist.resize(0) ;
      mesh.get_regn_vrtx(iR,plist) ;
      sort( plist.begin(), plist.end() ) ;
      if ( !binary_search( plist.begin(), plist.end(), iV ) ) {
	local_status = bool(false) ;
	cout << "-->>vertex: " << iV << " region: " << iR << " failed!!!" << endl ; 
      }
    }
    // check faces
    vector<int> vrtx_flist ;
    mesh.get_vrtx_face(iV,vrtx_flist) ;
    int nVF = vrtx_flist.size() ;
    for ( int ilF=0; ilF<nVF; ++ilF ) {
      int iF = vrtx_flist[ilF] ;
      plist.resize(0) ;
      mesh.get_face_vrtx(iF,plist) ;
      sort( plist.begin(), plist.end() ) ;
      if ( !binary_search( plist.begin(), plist.end(), iV ) ) {
	local_status = bool(false) ;
	cout << "-->>vertex: " << iV << " face: " << iF << " failed!!!" << endl ; 
      }
    }
    // check edges
    vector<int> vrtx_elist ;
    mesh.get_vrtx_edge(iV,vrtx_elist) ;
    int nVE = vrtx_elist.size() ;
    for ( int ilE=0; ilE<nVE; ++ilE ) {
      int iE = vrtx_elist[ilE] ;
      plist.resize(0) ;
      mesh.get_edge_vrtx(iE,plist) ;
      sort( plist.begin(), plist.end() ) ;
      if ( !binary_search( plist.begin(), plist.end(), iV ) ) {
	local_status = bool(false) ;
	cout << "-->>vertex: " << iV << " edge: " << iE << " failed!!!" << endl ; 
      }
    }
    // check vertices
    vector<int> vrtx_vlist ;
    mesh.get_vrtx_vrtx(iV,vrtx_vlist) ;
    int nVV = vrtx_vlist.size() ;
    for ( int jlV=0; jlV<nVV; ++jlV ) {
      int jV = vrtx_vlist[jlV] ;
      plist.resize(0) ;
      mesh.get_vrtx_vrtx(jV,plist) ;
      sort( plist.begin(), plist.end() ) ;
      if ( !binary_search( plist.begin(), plist.end(), iV ) ) {
	local_status = bool(false) ;
	cout << "-->>vertex: " << iV << " vrtx: " << jV << " failed!!!" << endl ; 
      }
    }
  }
  if ( local_status ) { 
    cout<<"--> OK"<<endl<<endl<<flush ;
  } else {
    cout<<endl<<">>>MESH VERTICES ARE NOT LOCALLY CONSISTENT!<<<"<<endl<<endl<<flush ;
  }
  return local_status ;
}

// loop  om all the mesh regions
// for every region, loop on the region's faces
// for every face:
// -if INTERNAL, search the same face from the region of the opposite side and
//               check that the face has (global) opposite orientation
// -if BOUNDARY, orientation must be posistive
bool mesh3Dv_cmaster :: check_regn_face_orientation() {
  cout<<"-->check regn_face orientation "<<endl<<flush ;
  bool local_status(true) ;
  int nR = mesh.n_region() ;
  for ( int iR=0; iR<nR; ++iR ) {                    // loop on mesh's regions
    int nRF = mesh.n_regn_face(iR) ;
    for ( int ilF=0; ilF<nRF; ++ilF ) {       // get region's faces
      int iF = mesh.regn_face(iR,ilF) ;
      if ( mesh.is_internal_face(iF) ) {      // F is an internal face:
	int iRp = mesh.ok_regn_face(iR,ilF) ? // check face orientation
	  mesh.face_regn(iF,1) :              // iF points out of iR
	  mesh.face_regn(iF,0) ;              // iF points into iR
	int nRpF = mesh.n_regn_face(iRp) ;
	for ( int ilpF=0; ilpF<nRpF; ++ilpF ) {  // orientation...
	  int ipF = mesh.regn_face(iRp,ilpF) ;
	  if ( ipF==iF ) {                       // ...got it! Check consistency!
	    local_status = mesh.ok_regn_face(iR,ilF) == !mesh.ok_regn_face(iRp,ilpF) ;
	    if ( !local_status ) {
	      cout << "-->>region: " << iR << " face (int): " << iF << " failed!!!" << endl ; 
	    }
	  }
	}
      } else {                                       // F is a boundary face:
	local_status = mesh.ok_regn_face(iR,ilF) ;   // it MUST point outward, check it! 
	if ( !local_status ) {
	  cout << "-->>region: " << iR << " face (bnd): " << iF << " failed!!!" << endl ; 
	}
      }
    }
  }
  if ( local_status ) { 
    cout<<"--> OK"<<endl<<endl<<flush ;
  } else {
    cout<<endl<<">>>FACES ARE NOT CONSISTENTLY ORIENTED IN REGIONS!<<<"<<endl<<endl<<flush ;
  }
  return local_status ;
}

// check that:
//   for every internal region the sequence of the faces is consistent with the 
//   sequence of the adjacent regions
bool mesh3Dv_cmaster :: check_regn_face_numbering() { 
  // the local numbering of region faces and neighbour regions should be the same
  // if the region is internal
  cout<<"-->check regn_face_regn numbering   "<<endl<<flush ;
  bool local_status(true) ;
  int nR = mesh.n_region() ;
  for ( int iR=0; iR<nR; ++iR ) {                    // loop on mesh's regions
    if ( mesh.is_internal_regn(iR) ) {
      vector<int> regn_rlist ;
      mesh.get_regn_regn(iR,regn_rlist) ;
      int nRF = mesh.n_regn_face(iR) ;
      for ( int ilF=0; ilF<nRF; ++ilF ) {              // get region's faces
	int iF = mesh.regn_face(iR,ilF) ;
	if ( mesh.is_internal_face(iF) ) {             // F is an internal face:
	  int iRp = mesh.ok_regn_face( iR, ilF ) ?     // check face orientation
	    mesh.face_regn(iF,1) :                     // iF points out of iR
	    mesh.face_regn(iF,0) ;                     // iF points into iR
	  local_status = iRp==regn_rlist[ilF] ;
	  if ( !local_status ) {
	    cout << "-->>region: " << iR << " face (int): " << iF << " failed!!!" << endl ; 
	  }
	} else {                                       // an internal rfegion cannot have a a boundary face!
	  cout << "-->>region: (int)" << iR << " face (bnd): " << iF << " failed!!!" << endl ; 
	}
      }
    }
  }
  if ( local_status ) { 
    cout<<"--> OK"<<endl<<endl<<flush ;
  } else {
    cout<< "\n>>>REGION FACES ARE NOT CONSISTENTLY NUMBERED WITH NEIGHBOUR REGIONS!<<<"<<endl <<endl<<flush ;
  }
  return local_status ;
}

// check that:
//   for every face the sequence of the edges is consistent with the 
//   sequence of the face vertices
bool mesh3Dv_cmaster :: check_face_edge_numbering() {
  cout<<"-->check face_edge numbering   "<<endl<<flush ;
  bool local_status(true) ;
  int nF = mesh.n_face() ;
  for ( int iF=0; iF<nF; ++iF ) {                    // loop on mesh's faces
    vector<int> face_vlist ;
    mesh.get_face_vrtx(iF,face_vlist) ;
    int nFE = mesh.n_face_edge(iF) ;
    for ( int ilE=0; ilE<nFE; ++ilE ) {        // get face's edges
      int iE  = mesh.face_edge(iF,ilE) ;       // get edge ID
      int iV0 = mesh.edge_vrtx(iE,0) ;         // get 1st edge vertex 
      int iV1 = mesh.edge_vrtx(iE,1) ;         // get 2nd edge vertex
      if ( !mesh.ok_face_edge(iF,ilE) ) {      // check orientation:
	iV0 = mesh.edge_vrtx(iE,1) ;           // swap the vertices if the edge
	iV1 = mesh.edge_vrtx(iE,0) ;           // has opposite orientation
      }
      local_status = (iV0==face_vlist[ilE]) && (iV1==face_vlist[(ilE+1)%nFE]) ;
      if ( !local_status ) {
	cout << "-->>face: " << iF << " edge: " << nFE << " failed!!!" << endl ; 
      }
    }
  }
  if ( local_status ) { 
    cout<<"--> OK"<<endl<<endl<<flush ;
  } else {
    cout<<endl<<">>>FACE EDGES AND VERTICES ARE NOT CONSISTENTLY ORIENTED IN FACES!<<<"<<endl<<endl<<flush ;
  }
  return local_status ;
}

bool mesh3Dv_cmaster :: check_vrtx_vrtx_numbering() {
  cout<<"-->check vrtx_vrtx numbering   "<<endl<<flush ;
  bool local_status(true) ;
  int nV = mesh.n_vertex() ;
  for ( int iV=0; iV<nV; ++iV ) {                    // loop on mesh's faces
    
    vector<int> vrtx_elist ;
    mesh.get_vrtx_edge(iV,vrtx_elist) ;
    int nVE = vrtx_elist.size() ;
    
    vector<int> vrtx_vlist ;
    mesh.get_vrtx_vrtx(iV,vrtx_vlist) ;
    int nVV = vrtx_vlist.size() ;

    for ( int jlV=0; jlV<nVV; ++jlV ) {
      int jV = vrtx_vlist[jlV] ;
      int iE = vrtx_elist[jlV] ;
      local_status = jV==mesh.edge_vrtx(iE,0) || jV==mesh.edge_vrtx(iE,1) ;
      // more sophisticated check using ok_vrtx_edge()
      //local_status = mesh.ok_vrtx_edge(iV,jlV) && jV==mesh.edge_vrtx(iE,1) ;
      if ( !local_status ) {
	cout << "-->>vertex: " << iV << " edge: " << iE << " failed!!!" << endl ; 
      }
    }
  }
  if ( local_status ) { 
    cout<<"--> OK"<<endl<<endl<<flush ;
  } else {
    cout<<endl<<">>>CONNECTED VERTICES AND EDGES ARE NOT CONSISTENTLY ORIENTED FOR VERTICES!<<<"<<endl<<endl<<flush ;
  }
  return local_status ;
}

#endif // end of _MESH3D_MASTER_HH
