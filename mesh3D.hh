/* ----------------------------------------------------------------/
This is open source software; you can redistribute it and/or modify
it under the terms of the BSD-3 License. If software is modified
to produce derivative works, such modified software should be
clearly marked, so as not to confuse it with the version available
from LANL. Full text of the BSD-3 License can be found in the
License file of the repository.
/---------------------------------------------------------------- */

#ifndef _MESH3DV_HH
#define _MESH3DV_HH

#include "dataset.hh"

const int UNSET = -999 ;

class mesh_3Dv {

  // friendships and type declarations 
  friend class vrtx_3Dv ;
  friend class edge_3Dv ;
  friend class face_3Dv ;
  friend class regn_3Dv ;

  friend class mesh3Dv_builder ;
  friend class mesh3Dv_printer ;
  friend class mesh3Dv_cmaster ;

  typedef TDataMat<double> dataMat ;
  typedef TDataMat<int>    dataInt ;
  typedef TDataVec<int>    dataVec ;

private:
  static const int DIM=3 ;

private:
  int nV, nE, nF, nR ;
  string mesh_name ;

  // coordinates
  dataMat V_coords ;

  // bounding box & mesh size
  bool   bb_status ;
  double bb_min[DIM], bb_max[DIM] ;
  
  bool   ms_status ;
  double ms_hmax, ms_hvol, ms_havg, ms_hmin ;
  
  // primary datasets
  dataInt EdgeVrtx ;
  dataset FaceEdge ;
  dataset RegnFace ;

  // tranposed datasets
  dataInt FaceRegn ;
  dataset EdgeFace ;
  dataset VrtxEdge ;

  // lists of boundary items
  dataVec bnd_vrtx ;
  dataVec bnd_edge ;
  dataVec bnd_face ;
  dataVec bnd_regn ;

  // external flags
  dataVec fV ;
  dataVec fE ;
  dataVec fF ;
  dataVec fR ;

  // aux geometrical quantities
  dataMat R_coords ;
  dataMat F_coords ;
  dataMat R_volume ;
  dataMat F_area   ;
  dataMat F_nor    ;

  // private method
  void shrink_list( vector<int> & tmp_list ) ;

public:
  mesh_3Dv() : bb_status(false), ms_status(false), mesh_name("mesh-3D") {}
  ~mesh_3Dv() {}

  // access methods
  int n_region() { return nR ; }
  int n_face  () { return nF ; }
  int n_edge  () { return nE ; }
  int n_vertex() { return nV ; }

  int n_bregion() { return bnd_regn.size() ; }
  int n_bface  () { return bnd_face.size() ; }
  int n_bedge  () { return bnd_edge.size() ; }
  int n_bvertex() { return bnd_vrtx.size() ; }

  int get_bnd_regn( int ilR ) ;
  int get_bnd_face( int ilF ) ;
  int get_bnd_edge( int ilE ) ;
  int get_bnd_vrtx( int ilV ) ;

  // TOPOLOGICAL METHODS:
  int regn_face( int iR, int ilF ) ;
  int face_edge( int iF, int ilE ) ;
  int edge_vrtx( int iE, int ilV ) ;

  int face_regn( int iF, int ilR ) ;
  int edge_face( int iE, int ilF ) ;
  int vrtx_edge( int iV, int ilE ) ;

  int vrtx_vrtx( int iV, int ilV ) ; // ISO vrtx_edge

  int n_regn_face( int iR ) ;
  int n_face_edge( int iF ) ;
  int n_edge_vrtx( int iE ) ;

  int n_face_regn( int iF ) ;
  int n_edge_face( int iE ) ;
  int n_vrtx_edge( int iV ) ;

  int n_vrtx_vrtx( int iV ) ;

  bool ok_regn_face( int iR, int ilF ) ;
  bool ok_face_edge( int iF, int ilE ) ;
  bool ok_edge_face( int iE, int ilF ) ;
  bool ok_vrtx_edge( int iV, int ilE ) ;

  // ...for regions
  void get_regn_regn( int iR, vector<int> & rlist ) ;
  void get_regn_face( int iR, vector<int> & flist ) ;
  void get_regn_edge( int iR, vector<int> & elist ) ;
  void get_regn_vrtx( int iR, vector<int> & vlist ) ;

  // ...for faces
  void get_face_regn( int iF, vector<int> & rlist ) ;
  void get_face_face( int iF, vector<int> & flist ) ;
  void get_face_edge( int iF, vector<int> & elist ) ;
  void get_face_vrtx( int iF, vector<int> & vlist ) ;

  // ...for edges
  void get_edge_regn( int iE, vector<int> & rlist ) ;
  void get_edge_face( int iE, vector<int> & flist ) ;
  void get_edge_edge( int iE, vector<int> & elist ) ;
  void get_edge_vrtx( int iE, vector<int> & vlist ) ;

  // ...for vertices
  void get_vrtx_regn( int iV, vector<int> & rlist ) ;
  void get_vrtx_face( int iV, vector<int> & flist ) ;
  void get_vrtx_edge( int iV, vector<int> & elist ) ;
  void get_vrtx_vrtx( int iV, vector<int> & vlist ) ;
  
  // Logical Methods fro detecting boundary items
  bool is_boundary_vrtx( int iV ) { return bnd_vrtx.find(iV) ; }
  bool is_boundary_edge( int iE ) { return bnd_edge.find(iE) ; }
  bool is_boundary_face( int iF ) { return bnd_face.find(iF) ; }
  bool is_boundary_regn( int iR ) { return bnd_regn.find(iR) ; }

  bool is_internal_vrtx( int iV ) { return !bnd_vrtx.find(iV) ; }
  bool is_internal_edge( int iE ) { return !bnd_edge.find(iE) ; }
  bool is_internal_face( int iF ) { return !bnd_face.find(iF) ; }
  bool is_internal_regn( int iR ) { return !bnd_regn.find(iR) ; }

  bool get_bnd_pos_vrtx( int iV ) { return bnd_vrtx.ipos(iV) ; }
  bool get_bnd_pos_edge( int iE ) { return bnd_edge.ipos(iE) ; }
  bool get_bnd_pos_face( int iF ) { return bnd_face.ipos(iF) ; }
  bool get_bnd_pos_regn( int iR ) { return bnd_regn.ipos(iR) ; }

  // Geometrical Methods (new)
  double get_nor( int iF, int s ) ;
  double get_tng( int iE, int s ) ;

  // Geometrical measures (new)
  double get_regn_measure( int iR ) ;
  double get_face_measure( int iF ) ;
  double get_edge_measure( int iF ) ;

  // eval bounding box
  inline double eval_bbox( string retstr ) ;
  inline void   eval_bbox() ;

  inline double min_coords( int s ) { 
    assert( 0<=s && s<DIM ) ;
    if ( !bb_status ) { eval_bbox() ; }
    return bb_min[s] ;
  }
  inline double max_coords( int s ) { 
    assert( 0<=s && s<DIM ) ;
    if ( !bb_status ) { eval_bbox() ; }
    return bb_max[s] ;
  }

  inline void bbox( double & xmin, double & ymin, double & zmin, 
		    double & xmax, double & ymax, double & zmax ) ;

  // used by problems
  inline double xmin()  { return bb_status ? bb_min[0] : eval_bbox("xmin") ; }
  inline double ymin()  { return bb_status ? bb_min[1] : eval_bbox("ymin") ; }
  inline double zmin()  { return bb_status ? bb_min[2] : eval_bbox("zmin") ; }
  inline double xmax()  { return bb_status ? bb_max[0] : eval_bbox("xmax") ; }
  inline double ymax()  { return bb_status ? bb_max[1] : eval_bbox("ymax") ; }
  inline double zmax()  { return bb_status ? bb_max[2] : eval_bbox("zmax") ; }

  // Geometrical Methods (with some problems)
  double coords_V( int iV, int k ) ;
  double coords_E( int iE, int k ) ;
  double coords_F( int iF, int k ) ;
  double coords_R( int iR, int k ) ;
  
  // arithmetic center of face iF
  double ari_coords_F( int iF, int k ) ; 

  // eval hmax
  void   eval_h() ;
  double h_max() ;
  double h_min() ;
  double h_vol() ;
  double h_avg() ;

  // mesh name
  void   set_mesh_name( string _mesh_name ) { mesh_name = _mesh_name ; } 
  string get_mesh_name( bool tex_flag=false ) { 
    string retval("") ;
    for ( int i=0; i<mesh_name.length(); ++i ) {
       if ( tex_flag && mesh_name[i]==char('_') ) { retval += "\\" ; }
       retval += mesh_name[i] ;
    }
    return retval ; 
  } 

  // set external flags (useful for gmv option "explode")
#if 1
  void set_fV( int iV, int new_fV ) { fV(iV) = new_fV ; }
  void set_fE( int iE, int new_fE ) { fE(iE) = new_fE ; }
  void set_fF( int iF, int new_fF ) { fF(iF) = new_fF ; }
  void set_fR( int iR, int new_fR ) { fR(iR) = new_fR ; }

  int  get_fV( int iV ) { return fV(iV) ; }
  int  get_fE( int iE ) { return fE(iE) ; }
  int  get_fF( int iF ) { return fF(iF) ; }
  int  get_fR( int iR ) { return fR(iR) ; }
#endif

// TEST 3
  void set_coords_V( int iV, double _x, double _y, double _z ) {
    V_coords(iV,0) = _x ;
    V_coords(iV,1) = _y ;
    V_coords(iV,2) = _z ;
  }
  void set_coords_F( int iF, double _x, double _y, double _z ) {
    F_coords(iF,0) = _x ;
    F_coords(iF,1) = _y ;
    F_coords(iF,2) = _z ;
  }
  void set_coords_R( int iR, double _x, double _y, double _z ) {
    R_coords(iR,0) = _x ;
    R_coords(iR,1) = _y ;
    R_coords(iR,2) = _z ;
  }
} ;
// --------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------
void mesh_3Dv :: shrink_list( vector<int> & tmp_list ) {
  sort( tmp_list.begin(), tmp_list.end() ) ;
  int k = 0 ;
  for ( int i=1; i<tmp_list.size(); ++i ) {
    if ( tmp_list[i]!=tmp_list[k] ) {
      tmp_list[++k]=tmp_list[i] ;
    }
  }
  tmp_list.resize(k+1) ;
}
// --------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------
  // TOPOLOGICAL METHODS:
int mesh_3Dv :: regn_face( int iR, int ilF ) {
  assert( 0<=iR  && iR<nR ) ;
  assert( 0<=ilF && ilF<RegnFace.size_loc(iR) ) ;
  return RegnFace(iR,ilF) ;
}
int mesh_3Dv :: face_edge( int iF, int ilE ) {
  assert( 0<=iF  && iF<nF ) ;
  assert( 0<=ilE && ilE<FaceEdge.size_loc(iF) ) ;
  return FaceEdge(iF,ilE) ;
}
int  mesh_3Dv :: edge_vrtx( int iE, int ilV ) {
  assert( 0<=iE  && iE<nE  ) ;
  assert( ilV==0 || ilV==1 ) ;
  return EdgeVrtx(iE,ilV) ;
}
int mesh_3Dv :: face_regn( int iF, int ilR ) {
  assert( 0<=iF  && iF<nF  ) ;
  assert( ilR==0 || ilR==1 ) ;
  return FaceRegn(iF,ilR) ;
}
int mesh_3Dv :: edge_face( int iE, int ilF ) {
  assert( 0<=iE  && iE<nE ) ;
  assert( 0<=ilF && ilF<EdgeFace.size_loc(iE) ) ;
  return EdgeFace(iE,ilF) ;
}
int  mesh_3Dv :: vrtx_edge( int iV, int ilE ) {
  assert( 0<=iV  && iV<nV ) ;
  assert( 0<=ilE && ilE<VrtxEdge.size_loc(iV) ) ;
  return VrtxEdge(iV,ilE) ;
}

int mesh_3Dv :: vrtx_vrtx( int iV, int ilE ) {
  assert( 0<=iV && iV<nV ) ;
  int iE  = VrtxEdge(iV,ilE) ; 
  int iV0 = EdgeVrtx(iE,0) ;
  int iV1 = EdgeVrtx(iE,1) ;
  return iV0==iV ? iV1 : iV0 ;
}

int mesh_3Dv :: n_regn_face( int iR ) {
  assert( 0<=iR  && iR<nR ) ;
  return RegnFace.size_loc(iR) ;
}
int mesh_3Dv :: n_face_edge( int iF ) {
  assert( 0<=iF  && iF<nF ) ;
  return FaceEdge.size_loc(iF) ;
}
int mesh_3Dv :: n_edge_vrtx( int iE ) {
  assert( 0<=iE  && iE<nE ) ;
  return EdgeVrtx.size_loc(iE) ;
}
// vrtx_vrtx ISO vrtx_edge
int mesh_3Dv :: n_vrtx_vrtx( int iV ) {
  assert( 0<=iV  && iV<nV ) ;
  return VrtxEdge.size_loc(iV) ;
}

int mesh_3Dv :: n_face_regn( int iF ) {
  assert( 0<=iF  && iF<nF ) ;
  return FaceRegn.size_loc(iF) ;
}
int mesh_3Dv :: n_edge_face( int iE ) {
  assert( 0<=iE  && iE<nE ) ;
  return EdgeFace.size_loc(iE) ;
}
int mesh_3Dv :: n_vrtx_edge( int iV ) {
  assert( 0<=iV  && iV<nV ) ;
  return VrtxEdge.size_loc(iV) ;
}

bool mesh_3Dv :: ok_regn_face( int iR, int ilF ) {
  assert( 0<=iR  && iR<nR ) ;
  assert( 0<=ilF && ilF<RegnFace.size_loc(iR) ) ;
  return RegnFace.get_value(iR,ilF) ;
}
bool mesh_3Dv :: ok_face_edge( int iF, int ilE ) {
  assert( 0<=iF  && iF<nF ) ;
  assert( 0<=ilE && ilE<FaceEdge.size_loc(iF) ) ;
  return FaceEdge.get_value(iF,ilE); 
}
bool mesh_3Dv :: ok_edge_face( int iE, int ilF ) {
  assert( 0<=iE  && iE<nE ) ;
  assert( 0<=ilF && ilF<EdgeFace.size_loc(iE) ) ;
  return EdgeFace.get_value(iE,ilF) ;
}
bool mesh_3Dv :: ok_vrtx_edge( int iV, int ilE ) {
  assert( 0<=iV  && iV<nV ) ;
  assert( 0<=ilE && ilE<VrtxEdge.size_loc(iV) ) ;
  return VrtxEdge.get_value(iV,ilE) ;
}
// --------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------
// access method of regions
void mesh_3Dv :: get_regn_regn( int iR, vector<int> & rlist ) {
  assert( 0<=iR && iR<nR ) ;
  for ( int ilF=0; ilF<RegnFace.size_loc(iR); ++ilF ) { 
    int iF  = RegnFace(iR,ilF) ;
    int iR0 = FaceRegn(iF,0) ;
    int iR1 = FaceRegn(iF,1) ;
    if      ( iR==iR1               ) { rlist.push_back( iR0 )   ; }
    else if ( iR==iR0 && iR1!=UNSET ) { rlist.push_back( iR1 )   ; }
    else if ( iR==iR0 && iR1==UNSET ) {} // boundary face
    else { // it should never happen
      MSG("-->>get_regn_regn: consistency error! "<<endl<<flush) ; 
      assert(false) ; 
    }
  }
}
void mesh_3Dv :: get_regn_face( int iR, vector<int> & flist ) {
  assert( 0<=iR && iR<nR ) ;
  flist.resize( RegnFace.size_loc(iR) ) ;
  for ( int ilF=0; ilF<RegnFace.size_loc(iR); ++ilF ) { flist[ilF] = RegnFace(iR,ilF) ; }
}
void mesh_3Dv :: get_regn_edge( int iR, vector<int> & elist ) {
  assert( 0<=iR  && iR<nR ) ;
  vector<int> tmp_vec ;
  for ( int ilF=0; ilF<RegnFace.size_loc(iR); ++ilF ) { 
    int iF = RegnFace(iR,ilF) ;
    for ( int ilE=0; ilE<FaceEdge.size_loc(iF); ++ilE ) {
      tmp_vec.push_back( FaceEdge(iF,ilE) ) ;
    } 
  }
  shrink_list( tmp_vec ) ;
  elist.resize( tmp_vec.size() ) ;
  for ( int i=0; i<elist.size(); ++i ) { elist[i]=tmp_vec[i] ; }
}
void mesh_3Dv :: get_regn_vrtx( int iR, vector<int> & vlist ) {
  assert( 0<=iR  && iR<nR ) ;
  vector<int> tmp_vec ;
  for ( int ilF=0; ilF<RegnFace.size_loc(iR); ++ilF ) { 
    int iF = RegnFace(iR,ilF) ;
    for ( int ilE=0; ilE<FaceEdge.size_loc(iF); ++ilE ) {
      int iE = FaceEdge(iF,ilE) ;
      tmp_vec.push_back( EdgeVrtx(iE,0) ) ;
      tmp_vec.push_back( EdgeVrtx(iE,1) ) ;
    } 
  }
  shrink_list( tmp_vec ) ;
  vlist.resize( tmp_vec.size() ) ;
  for ( int i=0; i<vlist.size(); ++i ) { vlist[i]=tmp_vec[i] ; }
}
// --------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------
// access method of faces
void mesh_3Dv :: get_face_regn( int iF, vector<int> & rlist ) {
  assert( 0<=iF && iF<nF ) ;
  rlist.resize(0) ;
  rlist.push_back(FaceRegn(iF,0)) ;
  if ( FaceRegn(iF,1)!=UNSET ) {
    rlist.push_back(FaceRegn(iF,1)) ;
  }
}
void mesh_3Dv :: get_face_face( int iF, vector<int> & flist ) {
  assert( 0<=iF  && iF<nF ) ;
  vector<int> tmp_vec ;
  for ( int ilE=0; ilE<FaceEdge.size_loc(iF); ++ilE ) {
    int iE = FaceEdge(iF,ilE) ;
    for ( int ilF=0; ilF<EdgeFace.size_loc(iE); ++ilF ) {
      int jF = EdgeFace(iE,ilF) ;
      if ( iF != jF ) { tmp_vec.push_back( jF ) ; }
    }
  } 
  shrink_list( tmp_vec ) ;
  flist.resize( tmp_vec.size() ) ;
  for ( int i=0; i<flist.size(); ++i ) { flist[i]=tmp_vec[i] ; }
}
void mesh_3Dv :: get_face_edge( int iF, vector<int> & elist ) {
  assert( 0<=iF && iF<nF ) ;
  elist.resize( FaceEdge.size_loc(iF) ) ;
  for ( int ilE=0; ilE<FaceEdge.size_loc(iF); ++ilE ) { elist[ilE] = FaceEdge(iF,ilE) ; }
}
void mesh_3Dv :: get_face_vrtx( int iF, vector<int> & vlist ) {
  assert( 0<=iF && iF<nF ) ;
  vector<int> tmp_vec ;
  for ( int ilE=0; ilE<FaceEdge.size_loc(iF); ++ilE ) {
    int iE = FaceEdge(iF,ilE) ;
    if ( FaceEdge.get_value(iF,ilE) ) {
      tmp_vec.push_back( EdgeVrtx(iE,0) ) ;
    } else {
      tmp_vec.push_back( EdgeVrtx(iE,1) ) ;
    }
  } 
  vlist.resize( tmp_vec.size() ) ;
  for ( int i=0; i<vlist.size(); ++i ) { vlist[i]=tmp_vec[i] ; }
}
// --------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------
// access method of edges
void mesh_3Dv :: get_edge_regn( int iE, vector<int> & rlist ) {
  assert( 0<=iE && iE<nE ) ;
  vector<int> tmp_vec ;
  for ( int ilF=0; ilF<EdgeFace.size_loc(iE); ++ilF ) { 
    int iF = EdgeFace(iE,ilF) ;
    tmp_vec.push_back(FaceRegn(iF,0)) ;
    if ( FaceRegn(iF,1)!=UNSET ) {
      tmp_vec.push_back(FaceRegn(iF,1)) ;
    }
  }
  shrink_list( tmp_vec ) ;
  rlist.resize( tmp_vec.size() ) ;
  for ( int i=0; i<rlist.size(); ++i ) { rlist[i]=tmp_vec[i] ; }
}
void mesh_3Dv :: get_edge_face( int iE, vector<int> & flist ) {
  assert( 0<=iE && iE<nE ) ;
  flist.resize( EdgeFace.size_loc(iE) ) ;
  for ( int ilF=0; ilF<EdgeFace.size_loc(iE); ++ilF ) { flist[ilF] = EdgeFace(iE,ilF) ; }
}
void mesh_3Dv :: get_edge_edge( int iE, vector<int> & elist ) {
  assert( 0<=iE && iE<nE ) ;
  vector<int> tmp_vec ;
  for ( int ilV=0; ilV<EdgeVrtx.size_loc(iE); ++ilV ) { 
    int iV = EdgeVrtx(iE,ilV) ;
    for ( int ilE=0; ilE<VrtxEdge.size_loc(iV); ++ilE ) {
      int jE = VrtxEdge(iV,ilE) ;
      if ( iE!=jE ) { tmp_vec.push_back( jE ) ; }
    }
  }
  shrink_list( tmp_vec ) ;
  elist.resize( tmp_vec.size() ) ;
  for ( int i=0; i<elist.size(); ++i ) { elist[i]=tmp_vec[i] ; }
}
void mesh_3Dv :: get_edge_vrtx( int iE, vector<int> & vlist ) {
  assert( 0<=iE && iE<nE ) ;
  vlist.resize( 2 ) ;
  vlist[0] = EdgeVrtx(iE,0) ;
  vlist[1] = EdgeVrtx(iE,1) ;
}
// --------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------
// access method of vertices
void mesh_3Dv :: get_vrtx_regn( int iV, vector<int> & rlist ) {
  assert( 0<=iV && iV<nV ) ;
  vector<int> tmp_vec ;  
  for ( int ilE=0; ilE<VrtxEdge.size_loc(iV); ++ilE ) { 
    int iE = VrtxEdge(iV,ilE) ; 
    for ( int ilF=0; ilF<EdgeFace.size_loc(iE); ++ilF ) { 
      int iF = EdgeFace(iE,ilF) ;
      tmp_vec.push_back(FaceRegn(iF,0)) ;
      if ( FaceRegn(iF,1)!=UNSET ) {
	tmp_vec.push_back(FaceRegn(iF,1)) ;
      }
    }
  }
  shrink_list( tmp_vec ) ;
  rlist.resize( tmp_vec.size() ) ;
  for ( int i=0; i<rlist.size(); ++i ) { rlist[i]=tmp_vec[i] ; }
}
void mesh_3Dv :: get_vrtx_face( int iV, vector<int> & flist ) {
  assert( 0<=iV && iV<nV ) ;
  vector<int> tmp_vec ;  
  for ( int ilE=0; ilE<VrtxEdge.size_loc(iV); ++ilE ) { 
    int iE = VrtxEdge(iV,ilE) ; 
    for ( int ilF=0; ilF<EdgeFace.size_loc(iE); ++ilF ) { 
      tmp_vec.push_back(EdgeFace(iE,ilF)) ; 
    }
  }
  shrink_list( tmp_vec ) ;
  flist.resize( tmp_vec.size() ) ;
  for ( int i=0; i<flist.size(); ++i ) { flist[i]=tmp_vec[i] ; }
}
void mesh_3Dv :: get_vrtx_edge( int iV, vector<int> & elist ) {
  assert( 0<=iV && iV<nV ) ;
  elist.resize( VrtxEdge.size_loc(iV) ) ;
  for ( int ilE=0; ilE<VrtxEdge.size_loc(iV); ++ilE ) { elist[ilE] = VrtxEdge(iV,ilE) ; }
}
void mesh_3Dv :: get_vrtx_vrtx( int iV, vector<int> & vlist ) {
  assert( 0<=iV && iV<nV ) ;
  vlist.resize( VrtxEdge.size_loc(iV) ) ;
  for ( int ilE=0; ilE<VrtxEdge.size_loc(iV); ++ilE ) { 
    int iE  = VrtxEdge(iV,ilE) ; 
    int iV0 = EdgeVrtx(iE,0) ;
    int iV1 = EdgeVrtx(iE,1) ;
    vlist[ilE] = iV0==iV ? iV1 : iV0 ;
  }
}
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
// (new)
double mesh_3Dv :: get_regn_measure( int iR ) {
  assert( 0<=iR && iR<nR ) ;
  assert( R_volume.size()==nR && R_volume.size_loc(iR)==1 ) ;
  return R_volume(iR,0) ;
}
double mesh_3Dv :: get_face_measure( int iF ) {
  assert( 0<=iF && iF<nF ) ;
  assert( F_area.size()==nF && F_area.size_loc(iF)==1 ) ;
  return F_area(iF,0) ;
}
double mesh_3Dv :: get_edge_measure( int iE ) {
  assert( 0<=iE && iE<nE ) ;
  int iV0 = edge_vrtx(iE,0) ;
  int iV1 = edge_vrtx(iE,1) ;
  return sqrt( pow(coords_V(iV1,0)-coords_V(iV0,0),2) + 
	       pow(coords_V(iV1,1)-coords_V(iV0,1),2) + 
	       pow(coords_V(iV1,2)-coords_V(iV0,2),2) ) ; 
}
double mesh_3Dv :: get_nor( int iF, int s ) {
  assert( 0<=iF && iF<nF ) ;
  assert( 0<=s  && s<DIM ) ;
  assert( F_nor.size()==nF && F_nor.size_loc(iF)==3 ) ;
  return F_nor(iF,s) ;
}
double mesh_3Dv :: get_tng( int iE, int s ) {
  assert( 0<=iE && iE<nE ) ;
  assert( 0<=s  && s<DIM ) ;
  int    iV0   = edge_vrtx(iE,0) ;
  int    iV1   = edge_vrtx(iE,1) ;
  double len_E = get_edge_measure(iE) ;
  return ( coords_V(iV1,s)-coords_V(iV0,s) )/len_E ;
}
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
#if 0
double mesh_3Dv :: get_regn_volume( int iR ) {
  assert( 0<=iR && iR<nR ) ;
  assert( R_volume.size()==nR && R_volume.size_loc(iR)==1 ) ;
  return R_volume(iR,0) ;
}
double mesh_3Dv :: get_face_area( int iF ) {
  assert( 0<=iF && iF<nF ) ;
  assert( F_area.size()==nF && F_area.size_loc(iF)==1 ) ;
  return F_area(iF,0) ;
}
double mesh_3Dv :: get_edge_length( int iE ) {
  assert( 0<=iE && iE<nE ) ;
  int iV0 = edge_vrtx(iE,0) ;
  int iV1 = edge_vrtx(iE,1) ;
  return sqrt( pow(coords_V(iV1,0)-coords_V(iV0,0),2) + 
	       pow(coords_V(iV1,1)-coords_V(iV0,1),2) + 
	       pow(coords_V(iV1,2)-coords_V(iV0,2),2) ) ; 
}
#endif
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
int mesh_3Dv :: get_bnd_regn( int ilR ) { 
  assert( 0<=ilR && ilR<bnd_regn.size() ) ;
  return bnd_regn(ilR) ;
}
int mesh_3Dv :: get_bnd_face( int ilF ) { 
  assert( 0<=ilF && ilF<bnd_face.size() ) ;
  return bnd_face(ilF) ;
}
int mesh_3Dv :: get_bnd_edge( int ilE ) { 
  assert( 0<=ilE && ilE<bnd_edge.size() ) ;
  return bnd_edge(ilE) ;
}
int mesh_3Dv :: get_bnd_vrtx( int ilV ) { 
  assert( 0<=ilV && ilV<bnd_vrtx.size() ) ;
  return bnd_vrtx(ilV) ;
}
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
double mesh_3Dv :: coords_V( int iV, int k ) {
  assert( 0<=iV && iV<nV ) ;
  assert( 0<=k  && k<3   ) ;
  return V_coords(iV,k) ;
}
double mesh_3Dv :: coords_E( int iE, int k ) {
  assert( 0<=iE && iE<nE ) ;
  assert( 0<=k  && k<3   ) ;
  double retval = 0. ; 
  for ( int ilV=0; ilV<EdgeVrtx.size_loc(iE); ++ilV ) { 
    int iV = EdgeVrtx(iE,ilV) ;
    retval += coords_V(iV,k) ;
  }
  return retval/double(EdgeVrtx.size_loc(iE)) ;
}
// this implementation is exact only for constants
double mesh_3Dv :: ari_coords_F( int iF, int k ) {
  assert( 0<=iF && iF<nF ) ;
  assert( 0<=k  && k<3   ) ;
  double retval = 0. ; 
  for ( int ilE=0; ilE<FaceEdge.size_loc(iF); ++ilE ) { 
    int iE = FaceEdge(iF,ilE) ;
    retval += coords_E(iE,k) ;
  }
  return retval/double(FaceEdge.size_loc(iF)) ;
}
double mesh_3Dv :: coords_F( int iF, int k ) {
  assert( 0<=iF && iF<nF ) ;
  assert( 0<=k  && k<3   ) ;
  return F_coords( iF, k ) ;
  //return ari_coords_F( iF, k ) ; // DEBUG
}
double mesh_3Dv :: coords_R( int iR, int k ) {
  assert( 0<=iR && iR<nR ) ;
  assert( 0<=k  && k<3   ) ;
  return R_coords( iR, k ) ;
}
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
void mesh_3Dv :: eval_bbox() {
  bb_status = true ;
  for ( int s=0; s<DIM; ++s ) { 
    bb_min[s] = +1e+99 ;
    bb_max[s] = -1e+99 ;
  }
  for ( int ilV=0; ilV<n_bvertex(); ++ilV ) {
    int iV = bnd_vrtx(ilV) ;
    for ( int s=0; s<DIM; ++s ) {
      bb_min[s] = min( bb_min[s], V_coords(iV,s) ) ;
      bb_max[s] = max( bb_max[s], V_coords(iV,s) ) ;
    }
  }
}
double mesh_3Dv :: eval_bbox( string retstr ) {
  bb_status = true ;
  for ( int s=0; s<DIM; ++s ) { 
    bb_min[s] = +1e+99 ;
    bb_max[s] = -1e+99 ;
  }
  for ( int ilV=0; ilV<n_bvertex(); ++ilV ) {
    int iV = bnd_vrtx(ilV) ;
    for ( int s=0; s<DIM; ++s ) {
    bb_min[s] = min( bb_min[s], V_coords(iV,s) ) ;
    bb_max[s] = max( bb_max[s], V_coords(iV,s) ) ;
    }
  }
  double retval ;
  if      ( retstr=="xmin" ) { retval = bb_min[0] ; }
  else if ( retstr=="ymin" ) { retval = bb_min[1] ; }
  else if ( retstr=="zmin" ) { retval = bb_min[2] ; }
  else if ( retstr=="xmax" ) { retval = bb_max[0] ; }
  else if ( retstr=="ymax" ) { retval = bb_max[1] ; }
  else if ( retstr=="zmax" ) { retval = bb_max[2] ; }
  else                       { retval = 0.        ; }
  return retval ;
}
void mesh_3Dv :: bbox( double & xmin, double & ymin, double & zmin, 
		       double & xmax, double & ymax, double & zmax ) {
  if ( !bb_status ) { eval_bbox() ; }
  xmin = bb_min[0] ;
  ymin = bb_min[1] ;
  zmin = bb_min[2] ;
  xmax = bb_max[0] ;
  ymax = bb_max[1] ;
  zmax = bb_max[2] ;
}
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
void mesh_3Dv :: eval_h() {
  if ( !ms_status ) { 
    ms_status = true ;
    // hmin
    ms_hmin = 1.e+99 ;
    for ( int iE=0; iE<nE; ++iE ) {
      int    iV0 = EdgeVrtx(iE,0) ;
      int    iV1 = EdgeVrtx(iE,1) ;
      double len = 
	sqrt( pow( V_coords(iV0,0)-V_coords(iV1,0),2 ) + 
	      pow( V_coords(iV0,1)-V_coords(iV1,1),2 ) + 
	      pow( V_coords(iV0,2)-V_coords(iV1,2),2 ) ) ;
      ms_hmax = max( len, ms_hmax ) ;
      ms_hmin = min( len, ms_hmin ) ;
    }

    // evaluare max distance of a polygon (from a Jerome suggestion)
    ms_hmax = 0. ;
    for ( int iR=0; iR<nR; ++iR ) {
      vector<int> R_vlist ;
      get_regn_vrtx( iR, R_vlist ) ;
      int nRV = R_vlist.size() ;
      for ( int ilV=0; ilV<nRV; ++ilV ) {
	int    iV  = R_vlist[ilV] ;
	double xiV = coords_V( iV, 0 ) ;
	double yiV = coords_V( iV, 1 ) ;
	double ziV = coords_V( iV, 2 ) ;
	for ( int jlV=ilV+1; jlV<nRV; ++jlV ) {
	  int    jV  = R_vlist[jlV] ;
	  double xjV = coords_V( jV, 0 ) ;
	  double yjV = coords_V( jV, 1 ) ;
	  double zjV = coords_V( jV, 2 ) ;
	  double dist = sqrt( pow( xiV-xjV, 2 ) + pow( yiV-yjV, 2 ) + pow( ziV-zjV, 2 ) ) ;
	  ms_hmax = max( ms_hmax,  dist ) ;
	}
      }
    }
  
    // hvol
    ms_hvol = 0. ;
    for ( int iR=0; iR<nR; ++iR ) {
      double hvol = pow( R_volume(iR,0), 1./3. ) ;
      ms_hvol = max( hvol, ms_hvol ) ;
    }
    // havg
    double sum = 0. ;
    for ( int iR=0; iR<nR; ++iR ) {
      sum += R_volume(iR,0) ;
    }
    ms_havg = pow( sum/nR, 1./double(DIM) ) ;    
  }
}
double mesh_3Dv :: h_max() {
  if ( !ms_status ) { eval_h() ; }
  return ms_hmax ;
}
double mesh_3Dv :: h_min() {
  if ( !ms_status ) { eval_h() ; }
  return ms_hmin ;
}
double mesh_3Dv :: h_vol() {
  if ( !ms_status ) { eval_h() ; }
  return ms_hvol ;
}
double mesh_3Dv :: h_avg() {
  if ( !ms_status ) { eval_h() ; }
  return ms_havg ;
}
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------
#if 0
#define mesh_item(ITEM,T,fZ,fZZ)					\
  class ITEM {								\
  private:								\
      int i##T ;							\
  protected:								\
      mesh_3Dv & mesh ;							\
  public:								\
      ITEM( int _i##T, mesh_3Dv & _mesh ) :				\
	i##T(_i##T), mesh(_mesh) {}					\
      ~ITEM() {}							\
      int get_id()       { return i##T ; }				\
      int get_ext_flag() { return mesh.f##T(i##T) ; }			\
      double x()         { return mesh.coords_##T(i##T,0) ; }		\
      double y()         { return mesh.coords_##T(i##T,1) ; }		\
      double z()         { return mesh.coords_##T(i##T,2) ; }		\
      bool is_internal() { return mesh.is_internal_##T(i##T) ; }	\
      int  n_##fZ()      { return mesh.n_##fZZ( i##T ) ; }		\
      int  ok_##fZ()     { return mesh.ok_##fZZ( i##T ) ; }		\
      int  fZ(int ilT)   { return mesh.fZZ( i##T, ilT ) ; }		\
  } ;

mesh_item(vrtx_3Dv,V,vrtx,vrtx_vrtx) ;
mesh_item(edge_3Dv,E,vrtx,edge_vrtx) ;
mesh_item(face_3Dv,F,edge,face_edge) ;
mesh_item(regn_3Dv,R,face,regn_face) ;
#endif
// -------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------

//#include "regn3D.hh"
//#include "face3D.hh"
//#include "edge3D.hh"
//#include "vrtx3D.hh"

#endif // end of _MESH3DV_HH
