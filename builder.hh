/* ----------------------------------------------------------------/
This is open source software; you can redistribute it and/or modify
it under the terms of the BSD-3 License. If software is modified
to produce derivative works, such modified software should be
clearly marked, so as not to confuse it with the version available
from LANL. Full text of the BSD-3 License can be found in the
License file of the repository.
/---------------------------------------------------------------- */

#define USE_TETGEN 0

#if USE_TETGEN
#include "quad_regn.hh"
#endif

class mesh3Dv_builder {

  static const int DIM = 3 ;
  static const bool sorting_regn_face = true ;

  // output for debugging purposes
  void print_vec_regn_face() ;

private:
  struct aux_struct {
  public:
    int  iGlb ;
    int  iloc ;
    bool bval ;
    aux_struct( int _iGlb, int _iloc, bool _bval=false ) : 
      iGlb(_iGlb), iloc(_iloc), bval(_bval) {}
    ~aux_struct() {}
    bool operator< ( const aux_struct & S1 ) const { return iGlb< S1.iGlb ; }
  } ;
  
  typedef struct aux_struct regn_face_struct ;
  typedef struct aux_struct face_edge_struct ;

private:
  mesh_3Dv & mesh ;

private:
  // build aux data structure
  void build_face_vlist( vector< vector<regn_face_struct> > & vec_regn_face,
			 vector<int>                        & face_vlist,
			 vector<int>                        & regn_flist ) ;
  
  void build_face_elist( vector< vector<face_edge_struct> > & vec_face_edge,
			 vector<int>                        & edge_vlist,
			 vector<int>                        & face_vlist ) ;

  // mesh coordinates
  void set_coords( vector<double> & xV, vector<double> & yV, vector<double> & zV ) ;
  void set_flags ( vector<int> & fV, vector<int> & fE, vector<int> & fF, vector<int> & fR ) ;

  // build primary data structures
  void build_RegnFace( vector< vector<regn_face_struct> > vec_regn_face, int _nR ) ;
  void build_FaceEdge( vector< vector<face_edge_struct> > vec_face_edge, int _nF ) ;
  void build_EdgeVrtx( vector<int> & edge_vlist ) ;

  // build transposed data structures
  void build_VrtxEdge() ;
  void build_EdgeFace() ;
  void build_FaceRegn() ;

  // build boundary lists
  void build_boundary_lists() ;
  void shrink_list( vector<int> & tmp_list ) ;

  // build geometric quantities
  void setup_geom_factors   () ;
  void set_regn_geom_factors() ;
  void set_face_geom_factors() ;
  void new_region_intg_values( int iR, vector<int> & vlist, double & vol_R, double & xR, double & yR, double & zR ) ;
  
  // check and reset face orientation (at face level)
  // (class check_master has a checking for orientation from inside regions)
  void fix_wrong_face_orientation() ;

public:
  mesh3Dv_builder( mesh_3Dv & _mesh ) : mesh(_mesh) {}
  ~mesh3Dv_builder() {}

  void change_bbox( double new_xmin, double new_ymin, double new_zmin,
		    double new_xmax, double new_ymax, double new_zmax ) ;

  void build_the_mesh( vector<double> & xV, vector<double> & yV, vector<double> & zV,
		       vector<int> & fV, 
		       vector<int> & regn_flist, vector<int> & fR ) {

    // -- set coordinates in mesh
    set_coords( xV, yV, zV ) ;
    
    // -- core for building primary dataset
    const int nR = fR.size() ;
    vector< vector<regn_face_struct> > vec_regn_face(nR) ;
    vector<int> face_vlist ;
    build_face_vlist( vec_regn_face, face_vlist, regn_flist ) ;
    build_RegnFace( vec_regn_face, nR ) ;
    
    const int nF = face_vlist.back() ;
    vector< vector<face_edge_struct> > vec_face_edge(nF) ;
    vector<int> edge_vlist ;
    build_face_elist( vec_face_edge, edge_vlist, face_vlist ) ;
    build_FaceEdge( vec_face_edge, nF ) ;
    build_EdgeVrtx( edge_vlist ) ;

    // -- build missing external flags, put flags into mesh
    int nE = edge_vlist.back() ;
    vector<int> fE(nE), fF(nF) ;
    for ( int iE=0; iE<nE; ++iE ) { fE[iE]=UNSET ; }
    for ( int iF=0; iF<nF; ++iF ) { fF[iF]=UNSET ; }
    set_flags( fV, fE, fF, fR ) ;

    // -- tranpose datasets
    build_VrtxEdge() ;
    build_EdgeFace() ;
    build_FaceRegn() ;

    // -- build boundary lists
    build_boundary_lists() ;
    
    // -- compute/set last geometric quantities
    setup_geom_factors() ;

    // -- reset wrong face orientation
    fix_wrong_face_orientation() ;
  }
} ; 

// setup coordinates
void mesh3Dv_builder ::
set_coords( vector<double> & xV, vector<double> & yV, vector<double> & zV ) {
  assert( xV.size()==yV.size() && yV.size()==zV.size() ) ;
  int nV = xV.size() ;
  mesh.nV = nV ;
  mesh.V_coords.setup(nV,3) ;
  for ( int iV=0; iV<nV; ++iV ) {
    mesh.V_coords(iV,0) = xV[iV] ;
    mesh.V_coords(iV,1) = yV[iV] ;
    mesh.V_coords(iV,2) = zV[iV] ;
  }
}
// setup external flags
void mesh3Dv_builder :: set_flags ( vector<int> & fV, vector<int> & fE, vector<int> & fF, vector<int> & fR ) {
  // --
  int nV = fV.size() ;
  mesh.fV.setup(nV) ;
  for ( int iV=0; iV<nV; ++iV ) { mesh.fV(iV) = fV[iV] ; }
  // --
  int nE = fE.size() ;
  mesh.fE.setup(nE) ;
  for ( int iE=0; iE<nE; ++iE ) { mesh.fE(iE) = fE[iE] ; }
  // --
  int nF = fF.size() ;
  mesh.fF.setup(nF) ;
  for ( int iF=0; iF<nF; ++iF ) { mesh.fF(iF) = fF[iF] ; }
  // --
  int nR = fR.size() ;
  mesh.fR.setup(nR) ;
  for ( int iR=0; iR<nR; ++iR ) { mesh.fR(iR) = fR[iR] ; }
  // --
}
// build up primary datasets
void mesh3Dv_builder :: build_RegnFace( vector< vector<regn_face_struct> > vec_regn_face, int _nR ) {
  MSG("begin build_RegnFace"<<endl<<flush) ;
  int nR = _nR ;
  mesh.nR = nR ;
  mesh.RegnFace.setup( nR ) ;
  for ( int iR=0; iR<nR; ++iR ) {
    if ( sorting_regn_face ) { 
      sort( vec_regn_face[iR].begin(), vec_regn_face[iR].end() ) ;
    }
    int nRF = vec_regn_face[iR].size() ;
    mesh.RegnFace.setup(iR,nRF) ;
    for ( int il=0; il<nRF; ++il ) {
      int  iF   = vec_regn_face[iR][il].iGlb ;
      int  ilF  = vec_regn_face[iR][il].iloc ;
      bool bval = vec_regn_face[iR][il].bval ;
      ilF = sorting_regn_face ? il : ilF ; 
      mesh.RegnFace.set_index(iR,ilF,iF) ;
      mesh.RegnFace.set_value(iR,ilF,bval) ;
    }
  }
  MSG("end-->build_RegnFace"<<endl<<flush) ;
}
void mesh3Dv_builder :: build_FaceEdge( vector< vector<face_edge_struct> > vec_face_edge, int _nF ) {
  MSG("begin build_FaceEdge"<<endl<<flush) ;
  int nF = _nF ;
  mesh.nF = nF ;
  mesh.FaceEdge.setup( nF ) ;
  for ( int iF=0; iF<nF; ++iF ) {
    int nFE = vec_face_edge[iF].size() ;
    mesh.FaceEdge.setup(iF,nFE) ;
    for ( int il=0; il<nFE; ++il ) {
      int  iE   = vec_face_edge[iF][il].iGlb ;
      int  ilE  = vec_face_edge[iF][il].iloc ;
      bool bval = vec_face_edge[iF][il].bval ;
      mesh.FaceEdge.set_index(iF,ilE,iE) ;
      mesh.FaceEdge.set_value(iF,ilE,bval) ;
    }
  }
  MSG("end-->build_FaceEdge"<<endl<<flush) ;
}
void mesh3Dv_builder :: build_EdgeVrtx( vector<int> & edge_vlist ) {
  MSG("begin build_EdgeVrtx"<<endl<<flush) ;
  int nE = edge_vlist.back() ;
  mesh.nE = nE ;
  mesh.EdgeVrtx.setup( nE, 2 ) ;
  int k = 0 ;
  for ( int iE=0; iE<nE; ++iE ) {
    int nEV = edge_vlist[k++] ;
    int kV0 = edge_vlist[k++] ;
    int kV1 = edge_vlist[k++] ;
    assert( nEV==2 ) ;
    mesh.EdgeVrtx(iE,0) = kV0 ;
    mesh.EdgeVrtx(iE,1) = kV1 ;
  }
  MSG("end-->build_EdgeVrtx"<<endl<<flush) ;
}
// build up transposed datasets
void mesh3Dv_builder :: build_VrtxEdge() {
  MSG("begin build_VrtxEdge"<<endl<<flush) ;
  const int nV = mesh.nV ;
  const int nE = mesh.nE ;
  PRT(nV) ;
  vector< vector<aux_struct> > vec_aux(nV) ;
  for ( int iE=0; iE<nE; ++iE ) {
    int nEV = 2 ;
    for ( int ilV=0; ilV<nEV; ++ilV ) {
      int iV = mesh.EdgeVrtx(iE,ilV) ;
      assert( 0<=iV && iV<nV ) ;
      vec_aux[iV].push_back( aux_struct(iE,ilV,bool(ilV==0)) ) ;
    }
  }
  mesh.VrtxEdge.setup( nV ) ;
  for ( int iV=0; iV<nV; ++iV ) {
    int nVE = vec_aux[iV].size() ;
    mesh.VrtxEdge.setup( iV, nVE ) ;
    for ( int ilE=0; ilE<nVE; ++ilE ) {
      int  iE   = vec_aux[iV][ilE].iGlb ;
      bool bval = vec_aux[iV][ilE].bval ;
      assert( 0<=iE && iE<nE ) ;
      mesh.VrtxEdge.set_index( iV, ilE, iE   ) ;
      mesh.VrtxEdge.set_value( iV, ilE, bval );
    }
  }
  MSG("end-->build_VrtxEdge"<<endl<<flush) ;
}
void mesh3Dv_builder :: build_EdgeFace() {
  MSG("begin build_EdgeFace"<<endl<<flush) ;
  const int nE = mesh.nE ;
  const int nF = mesh.nF ;

  vector< vector<aux_struct> > vec_aux(nE) ;

  for ( int iF=0; iF<nF; ++iF ) {
    int nFE = mesh.FaceEdge.size_loc(iF) ;
    for ( int ilE=0; ilE<nFE; ++ilE ) {
      int iE = mesh.FaceEdge(iF,ilE) ;
      assert( 0<=iE && iE<nE ) ;
      bool bval = mesh.FaceEdge.get_value(iF,ilE) ;
      vec_aux[iE].push_back( aux_struct(iF,ilE,bval) ) ;
    }
  }
  mesh.EdgeFace.setup( nE ) ;
  for ( int iE=0; iE<nE; ++iE ) {
    int nEF = vec_aux[iE].size() ;
    mesh.EdgeFace.setup( iE, nEF ) ;
    for ( int ilF=0; ilF<nEF; ++ilF ) {
      int  iF   = vec_aux[iE][ilF].iGlb ;
      bool bval = vec_aux[iE][ilF].bval ;
      mesh.EdgeFace.set_index( iE, ilF, iF   ) ;
      mesh.EdgeFace.set_value( iE, ilF, bval );
    }
  }
  MSG("end-->build_EdgeFace"<<endl<<flush) ;
}
void mesh3Dv_builder :: build_FaceRegn() {
  MSG("begin build_FaceRegn"<<endl<<flush) ;
  const int nF = mesh.nF ;
  const int nR = mesh.nR ;
  vector< vector<aux_struct> > vec_aux(nF) ;
  for ( int iR=0; iR<nR; ++iR ) {
    int nRF = mesh.RegnFace.size_loc(iR) ;
    for ( int ilF=0; ilF<nRF; ++ilF ) {
      int  iF   = mesh.RegnFace(iR,ilF) ;
      bool bval = mesh.RegnFace.get_value(iR,ilF) ;
      vec_aux[iF].push_back( aux_struct(iR,ilF,bval) ) ;
      if ( vec_aux[iF].size()!=1 && vec_aux[iF].size()!=2 ) {
	MSG("Error in construction of FaceRegn dataset. Face: "<<endl) ;
	PRT(iF) ;
	MSG("belongs to more than 1 or 2 cells: ") ;
	int nFR = vec_aux[iF].size() ;
	for ( int ilR=0; ilR<nFR; ++ilR ) {
	  int  iR   = vec_aux[iF][ilR].iGlb ;
	  //bool bval = vec_aux[iF][ilR].bval ;
	  cout << "   "  << iR ;
	}
	cout << endl ;
	LINE(---) ;
      }
    }
  }
  mesh.FaceRegn.setup( nF, 2, UNSET ) ;
  for ( int iF=0; iF<nF; ++iF ) {
    int nFR = vec_aux[iF].size() ;
    assert( nFR==1 || nFR==2 ) ; 
    for ( int ilR=0; ilR<nFR; ++ilR ) {
      int  iR   = vec_aux[iF][ilR].iGlb ;
      //bool bval = vec_aux[iF][ilR].bval ;
      mesh.FaceRegn( iF, ilR ) = iR ;
    }
  }
  MSG("end-->build_FaceRegn"<<endl<<flush) ;
}
void mesh3Dv_builder :: shrink_list( vector<int> & tmp_list ) {
  sort( tmp_list.begin(), tmp_list.end() ) ;
  int k = 0 ;
  for ( int i=1; i<tmp_list.size(); ++i ) {
    if ( tmp_list[i]!=tmp_list[k] ) {
      tmp_list[++k]=tmp_list[i] ;
    }
  }
  tmp_list.resize(k+1) ;
}
void mesh3Dv_builder :: build_boundary_lists() {
  MSG("begin build_boundary_lists"<<endl<<flush) ;
  // introduce tmp vectors
  vector<int> tmp_vrtx ;
  vector<int> tmp_face ;
  vector<int >tmp_edge ;
  vector<int> tmp_regn ;
  
  // gather all boundary items
  for ( int iF=0; iF<mesh.nF; ++iF ) {
    if ( mesh.FaceRegn(iF,1)==UNSET ) {
      tmp_face.push_back( iF ) ;
      tmp_regn.push_back( mesh.FaceRegn(iF,0) ) ;
      for ( int ilE=0; ilE<mesh.FaceEdge.size_loc(iF); ++ilE ) {
	int iE = mesh.FaceEdge(iF,ilE) ;
	tmp_edge.push_back( iE ) ;
	tmp_vrtx.push_back( mesh.EdgeVrtx(iE,0) ) ;
	tmp_vrtx.push_back( mesh.EdgeVrtx(iE,1) ) ;
      }
    }
  }

  // shrink
  shrink_list( tmp_vrtx ) ;
  shrink_list( tmp_edge ) ;
  shrink_list( tmp_face ) ;
  shrink_list( tmp_regn ) ;

  // setup and copy
  mesh.bnd_vrtx.setup( tmp_vrtx.size() ) ;
  mesh.bnd_edge.setup( tmp_edge.size() ) ;
  mesh.bnd_face.setup( tmp_face.size() ) ;
  mesh.bnd_regn.setup( tmp_regn.size() ) ;

  for ( int i=0; i<tmp_vrtx.size(); ++i ) { mesh.bnd_vrtx(i) = tmp_vrtx[i] ; }
  for ( int i=0; i<tmp_edge.size(); ++i ) { mesh.bnd_edge(i) = tmp_edge[i] ; }
  for ( int i=0; i<tmp_face.size(); ++i ) { mesh.bnd_face(i) = tmp_face[i] ; }
  for ( int i=0; i<tmp_regn.size(); ++i ) { mesh.bnd_regn(i) = tmp_regn[i] ; }

  MSG("end-->build_boundary_lists"<<endl<<flush) ;
}

#include "build_face.hh"
#include "build_edge.hh"
#include "build_geom.hh"

// DEBUG
#if 0
void mesh3Dv_builder :: print_vec_regn_face() {
  int nR = vec_regn_face.size() ;
  for ( int iR=0; iR<nR; ++iR ) {
    regn_face_struct & rfs = vec_regn_face[iR] ;
    int iF   = rfs.iGlb ;
    int ilF  = rfs.iloc ;
    int bval = rfs.bval ;
  } 
}
#endif

// this method  is required when the input mesh does not necessarily satisfy the
// criterion that face orientation points for region "0" to region "1"
// (for example: the *.msh format provided by F. Hubert).
void mesh3Dv_builder :: fix_wrong_face_orientation() {
  // set "verbosity level"
  int noisy = 1 ;

  // this routine should be run after mesh construction is terminated
  int nF = mesh.n_face() ;
  for ( int iF=0; iF<nF; ++iF ) {

    bool need_fixing = false ;

    if ( mesh.is_internal_face(iF) ) {
    
      int iR0 = mesh.FaceRegn( iF, 0 ) ;
      int iR1 = mesh.FaceRegn( iF, 1 ) ;
    
      double nor_F[3], vFR_0[3], vFR_1[3] ;
      for ( int s=0; s<3; ++s ) {
	vFR_0[s] = mesh.coords_F( iF, s ) - mesh.coords_R( iR0, s ) ;
	vFR_1[s] = mesh.coords_F( iF, s ) - mesh.coords_R( iR1, s ) ;
	nor_F[s] = mesh.get_nor ( iF,  s ) ;
      }

      double sgn_0  = vFR_0[0]*nor_F[0] + vFR_0[1]*nor_F[1] + vFR_0[2]*nor_F[2]>0 ? 1. : -1 ;
      double sgn_1  = vFR_1[0]*nor_F[0] + vFR_1[1]*nor_F[1] + vFR_1[2]*nor_F[2]>0 ? 1. : -1 ;

      if ( sgn_0 * sgn_1 > 0. ) {
	MSGF("A MAJOR ERROR OCCURRED IN MESH CONSTRUCTION:"<<endl) ;
	MSGF("FACE ORIENTATION IS WRONG AND CANNOT BE RECOVERED"<<endl) ;

	PRT(iF) ;
	PRT(iR0) ;
	PRT(iR1) ;

	for ( int s=0; s<3; ++s ) {
	  VALM(iF,s,mesh.coords_F) ; VALM(iR0,s,mesh.coords_R) ; PRTM(iR1,s,mesh.coords_R) ; 
	}

	PRT( vFR_0[0]*nor_F[0] + vFR_0[1]*nor_F[1] + vFR_0[2]*nor_F[2] ) ;
	PRT( vFR_1[0]*nor_F[0] + vFR_1[1]*nor_F[1] + vFR_1[2]*nor_F[2] ) ;

	assert(sgn_0*sgn_1<0.) ;
      }

      if ( sgn_0<0. && sgn_1>0 ) {
	need_fixing = true ;
	if ( noisy>0 ) { 
	  MSG("change orientation of internal face:") ;
	  VAL(iF) ; VAL(iR0) ; PRT(iR1) ;
	}
      } else {
	// nothing to be done, the face is properly oriented
      }
      
    } else if ( mesh.is_boundary_face(iF) ) {
      
      int iR0 = mesh.FaceRegn( iF, 0 ) ;
      
      double nor_F[3], vFR_0[3] ;
      for ( int s=0; s<3; ++s ) {
	vFR_0[s] = mesh.coords_F( iF, s ) - mesh.coords_R( iR0, s ) ;
	nor_F[s] = mesh.get_nor ( iF,  s ) ;
      }
      
      double sgn_0 = vFR_0[0]*nor_F[0] + vFR_0[1]*nor_F[1] + vFR_0[2]*nor_F[2]>0 ? 1. : -1 ;
      
      if ( sgn_0<0 ) {
	need_fixing = true ;
	if ( noisy>0 ) { 
	  MSG("change orientation of boundary face: ") ; 
	  VAL(iF) ; VAL(iR0) ; PRT(sgn_0) ;
	}
      }
    }

    // reverse face orientation
    if ( need_fixing ) {
 
      vector<int> elist ;
      mesh.get_face_edge( iF, elist ) ;
      int nFE = elist.size() ;

      vector<bool> blist(nFE) ;
      for ( int ilE=0; ilE<nFE; ++ilE ) {
	blist[ilE] = mesh.FaceEdge.get_value( iF, ilE ) ; 
      }

      // reverse FaceEdge
      for ( int ilE=0; ilE<nFE; ++ilE ) {
	int  iE = elist[nFE-1-ilE] ;
	bool bE = blist[nFE-1-ilE] ;
	mesh.FaceEdge.set_index( iF, ilE,  iE ) ;
	mesh.FaceEdge.set_value( iF, ilE, !bE ) ;
      }

      // reverse the value of EdgeFace
      for ( int ilE=0; ilE<nFE; ++ilE ) {
	int  iE = elist[nFE-1-ilE] ;
	int nEF = mesh.n_edge_face(iE) ;
	for ( int ilF=0; ilF<nEF; ++ilF ) {
	  if ( mesh.EdgeFace(iE,ilF)==iF ) {
	    bool bval = mesh.EdgeFace.get_value( iE, ilF ) ;
	    mesh.EdgeFace.set_value( iE, ilF, !bval ) ;
	  }
	}
      }

      for ( int s=0; s<3; ++s ) {
	mesh.F_nor( iF, s ) *= -1. ;
      }

    } // end of --> if ( need_fixing 
  } // end of --> for ( int iF=0; ...
}


// fine DEBUG
