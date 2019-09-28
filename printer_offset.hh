/* ----------------------------------------------------------------/
This is open source software; you can redistribute it and/or modify
it under the terms of the BSD-3 License. If software is modified
to produce derivative works, such modified software should be
clearly marked, so as not to confuse it with the version available
from LANL. Full text of the BSD-3 License can be found in the
License file of the repository.
/---------------------------------------------------------------- */

class mesh3Dv_printer {

private:
  mesh_3Dv & mesh ;
  int offset ;

public:
  mesh3Dv_printer( mesh_3Dv & _mesh, int _offset=0 ) : mesh(_mesh), offset(_offset) {}
  ~mesh3Dv_printer() {}

  // print DATASETS
  void print_EdgeVrtx( ostream & LOGF ) ;
  void print_FaceEdge( ostream & LOGF ) ;
  void print_RegnFace( ostream & LOGF ) ;
  // ---
  void print_FaceRegn( ostream & LOGF ) ;
  void print_EdgeFace( ostream & LOGF ) ;
  void print_VrtxEdge( ostream & LOGF ) ;
  // ---
  void print_boundary_lists( ostream & LOGF ) ;
  // --- 
  void print_all_datasets() ;

  // print REGIONS
  void print_all_regions() ;
  void print_region( int iR, ostream & LOGF ) ;

  // print FACES
  void print_all_faces() ;
  void print_face( int iF, ostream & LOGF ) ;

  // print EDGES
  void print_all_edges() ;
  void print_edge( int iE, ostream & LOGF ) ;

  // print VERTICES
  void print_all_vertices() ;
  void print_vertex( int iV, ostream & LOGF ) ;
} ;

// --------------------------------------------------------------------------------------------
void mesh3Dv_printer :: print_EdgeVrtx( ostream & LOGF ) {
  // EdgeVrtx
  LOGF << "EdgeVrtx " << endl ;
  LOGF << "#of vertices = "  << mesh.EdgeVrtx.size() << endl ;
  for ( int iE=0; iE<mesh.EdgeVrtx.size(); ++iE ) {
    LOGF << "iE = " << iE+offset << ", { " ;
    if ( mesh.EdgeVrtx.size_loc(iE)>0 ) {
      for ( int ilV=0; ilV<mesh.EdgeVrtx.size_loc(iE)-1; ++ilV ) {
	LOGF << mesh.EdgeVrtx(iE,ilV)+offset << ", " ;
      }
      LOGF << mesh.EdgeVrtx(iE,mesh.EdgeVrtx.size_loc(iE)-1)+offset ;
    }
    LOGF << " } " << endl ;
  }
  LOGF << "-------------------" << endl ;
}
void mesh3Dv_printer :: print_FaceEdge( ostream & LOGF ) {
  // FaceEdge
  LOGF << "FaceEdge " << endl ;
  LOGF << "#of faces = "  << mesh.FaceEdge.size() << endl ;
  for ( int iF=0; iF<mesh.FaceEdge.size(); ++iF ) {
    LOGF << "iF = " << iF+offset << ", { "  ;
    if ( mesh.FaceEdge.size_loc(iF)>0 ) {
      for ( int j=0; j<mesh.FaceEdge.size_loc(iF)-1; ++j ) {
	LOGF << "(" << mesh.FaceEdge.get_value(iF,j) << ")" << mesh.FaceEdge(iF,j)+offset << ", " ;
      }
      int j = mesh.FaceEdge.size_loc(iF)-1 ;
      LOGF << "(" << mesh.FaceEdge.get_value(iF,j) << ")" 
	   << mesh.FaceEdge(iF,j)+offset ;
    }
    LOGF << " } " << endl ;
  }
  LOGF << "-------------------" << endl ;
}
void mesh3Dv_printer :: print_RegnFace( ostream & LOGF ) {
  // RegnFace
  LOGF << "RegnFace " << endl ;
  LOGF << "#of regns = "  << mesh.RegnFace.size() << endl ;
  for ( int iR=0; iR<mesh.RegnFace.size(); ++iR ) {
    LOGF << "iR = " << iR+offset << ", { "  ;
    if ( mesh.RegnFace.size_loc(iR)>0 ) {
      for ( int j=0; j<mesh.RegnFace.size_loc(iR)-1; ++j ) {
	LOGF << "(" << mesh.RegnFace.get_value(iR,j) << ")" << mesh.RegnFace(iR,j)+offset << ", " ;
      }
      int j = mesh.RegnFace.size_loc(iR)-1 ;
      LOGF << "(" << mesh.RegnFace.get_value(iR,j) << ")" 
	   << mesh.RegnFace(iR,j)+offset ;
    }
    LOGF << " } " << endl ;
  }
  LOGF << "-------------------" << endl ;
}
// ------------------------------------------------------------------------------------------
void mesh3Dv_printer :: print_VrtxEdge( ostream & LOGF ) {
  // VrtxEdge
  LOGF << "VrtxEdge " << endl ;
  LOGF << "#of vertices = "  << mesh.VrtxEdge.size() << flush << endl ;
  for ( int iV=0; iV<mesh.VrtxEdge.size(); ++iV ) {
    LOGF << "iV = " << iV+offset << ", { " ;
    if ( mesh.VrtxEdge.size_loc(iV)>0 ) {
      for ( int j=0; j<mesh.VrtxEdge.size_loc(iV)-1; ++j ) {
	LOGF << "(" << mesh.VrtxEdge.get_value(iV,j) << ")" << mesh.VrtxEdge(iV,j)+offset << ", " ;
      }
      int j = mesh.VrtxEdge.size_loc(iV)-1 ;
      LOGF << "(" << mesh.VrtxEdge.get_value(iV,j) << ")" 
	   << mesh.VrtxEdge(iV,j)+offset ; 
    }
    LOGF << " } " << endl ;
  }
  LOGF << "-------------------" << endl ;
}
void mesh3Dv_printer :: print_EdgeFace( ostream & LOGF ) {
  // EdgeFace
  LOGF << "EdgeFace " << endl ;
  LOGF << "#of edges = "  << mesh.EdgeFace.size() << endl ;
  for ( int iE=0; iE<mesh.EdgeFace.size(); ++iE ) {
    LOGF << "iE = " << iE+offset << ", { "  ;
    if ( mesh.EdgeFace.size_loc(iE)>0 ) {
      for ( int j=0; j<mesh.EdgeFace.size_loc(iE)-1; ++j ) {
	LOGF << "(" << mesh.EdgeFace.get_value(iE,j) << ")" << mesh.EdgeFace(iE,j)+offset << ", " ;
      }
      int j = mesh.EdgeFace.size_loc(iE)-1 ;
      LOGF << "(" << mesh.EdgeFace.get_value(iE,j) << ")" 
	   << mesh.EdgeFace(iE,j)+offset ;
    }
    LOGF << " } " << endl ;
  }
  LOGF << "-------------------" << endl ;
}
void mesh3Dv_printer :: print_FaceRegn( ostream & LOGF ) {
  // EdgeFace
  LOGF << "FaceRegn " << endl ;
  LOGF << "#of faces = "  << mesh.FaceRegn.size() << endl ;
  for ( int iF=0; iF<mesh.FaceRegn.size_loc(iF); ++iF ) {
    LOGF << "iF = " << iF+offset << ", { " 
	 << mesh.FaceRegn(iF,0)+offset << ", " << mesh.FaceRegn(iF,1)+offset ;
  }
  LOGF << " } " << endl ;
  LOGF << "-------------------" << endl ;
}

void mesh3Dv_printer :: print_boundary_lists( ostream & LOGF ) {
  // boundary lists
  LOGF << "boundary lists " << endl ;
  LOGF << "#of boundary vertices = "  << mesh.bnd_vrtx.size() << endl ;
  for ( int i=0; i<mesh.bnd_vrtx.size(); ++i ) {
    LOGF << "bnd_vrtx[" << i << "] = " << mesh.bnd_vrtx(i)+offset << endl ; 
  }
  LOGF << "-------------------" << endl ;
  LOGF << "#of boundary edges = "  << mesh.bnd_edge.size() << endl ;
  for ( int i=0; i<mesh.bnd_edge.size(); ++i ) {
    LOGF << "bnd_edge[" << i << "] = " << mesh.bnd_edge(i)+offset << endl ; 
  }
  LOGF << "-------------------" << endl ;
  LOGF << "#of boundary faces = "  << mesh.bnd_face.size() << endl ;
  for ( int i=0; i<mesh.bnd_face.size(); ++i ) {
    LOGF << "bnd_face[" << i << "] = " << mesh.bnd_face(i)+offset << endl ; 
  }
  LOGF << "-------------------" << endl ;
  LOGF << "#of boundary regions = "  << mesh.bnd_regn.size() << endl ;
  for ( int i=0; i<mesh.bnd_regn.size(); ++i ) {
    LOGF << "bnd_regn[" << i << "] = " << mesh.bnd_regn(i)+offset << endl ; 
  }
  LOGF << "-------------------" << endl ;
}
// --------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------
void mesh3Dv_printer :: print_all_datasets() {
  ofstream LOGF("dset.log") ;
  print_EdgeVrtx( LOGF ) ;
  print_FaceEdge( LOGF ) ;
  print_RegnFace( LOGF ) ;
  print_EdgeFace( LOGF ) ;
  print_VrtxEdge( LOGF ) ;
  print_FaceRegn( LOGF ) ;
  LOGF << "-------------------" << endl ;
  print_boundary_lists( LOGF ) ;
  LOGF.close() ;
}
// --------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------
# define print_list_offset(VEC)                                     \
  do {                                                       \
    LOGF << " { " ;					     \
    if ( VEC.size()>0 ) {				     \
      for ( int i=0 ; i<VEC.size()-1 ; ++i ) {		     \
	LOGF << VEC[i]+offset << ", " ;			     \
      }							     \
      LOGF << VEC.back()+offset ;				     \
    }							     \
    LOGF << " }" << endl ;				     \
  } while(0)
// --------------------------------------------------------------------------------------------
// REGIONS
// --------------------------------------------------------------------------------------------
void mesh3Dv_printer :: print_all_regions() {
  ofstream LOGF("region.log") ;
  LOGF << "#of regions: nR=" << mesh.n_region() << endl ;
  for ( int iR=0; iR<mesh.n_region(); ++iR ) {
    print_region( iR, LOGF ) ;
  }
  LOGF.close() ;
}
void mesh3Dv_printer :: print_region( int iR, ostream & LOGF ) {
  LOGF << "--------------------------------------------------------------" << endl ;
  LOGF << "region: iR =" << iR+offset << endl ;
  // faces
  int nRF = mesh.RegnFace.size_loc(iR);
  LOGF << "#faces: nRF=" << nRF << ", { " ;
  for ( int ilF=0; ilF<nRF-1; ++ilF ) {
    LOGF << "(" 
	 << mesh.ok_regn_face(iR,ilF) << ")" 
	 << mesh.regn_face   (iR,ilF)+offset << ", " ;
  }
  LOGF << "(" 
       << mesh.ok_regn_face(iR,nRF-1) << ")" 
       << mesh.regn_face   (iR,nRF-1)+offset << " }" << endl ;
  LOGF << flush ;
  // topological info
  vector<int> rlist, flist, elist, vlist ;
  mesh.get_regn_regn( iR, rlist ) ;
  mesh.get_regn_face( iR, flist ) ;
  mesh.get_regn_edge( iR, elist ) ;
  mesh.get_regn_vrtx( iR, vlist ) ;
  LOGF << "---" << endl ;
  LOGF << "#regions:  nRR=" << rlist.size() ; print_list_offset(rlist) ;
  LOGF << "#faces:    nRF=" << flist.size() ; print_list_offset(flist) ;
  LOGF << "#edges:    nRE=" << elist.size() ; print_list_offset(elist) ;
  LOGF << "#vertices: nRV=" << vlist.size() ; print_list_offset(vlist) ;
  // geometrical info
  LOGF << "---" << endl ;
  LOGF << "#center: ( " 
       << mesh.coords_R(iR,0) << ", " << mesh.coords_R(iR,1) << ", " 
       << mesh.coords_R(iR,2) << " )" << endl << flush ; 
  LOGF << "#volume:     " << mesh.get_regn_measure(iR) << endl << flush ; 
}
// --------------------------------------------------------------------------------------------
// FACES
// --------------------------------------------------------------------------------------------
void mesh3Dv_printer :: print_all_faces() {
  ofstream LOGF("face.log") ;
  LOGF << "#of faces: nF=" << mesh.n_face() << endl ;
  for ( int iF=0; iF<mesh.n_face(); ++iF ) {
    print_face( iF, LOGF ) ;
  }
  LOGF.close() ;
}
void mesh3Dv_printer :: print_face( int iF, ostream & LOGF ) {
  LOGF << "--------------------------------------------------------------" << endl ;
  LOGF << "face: iF =" << iF+offset << endl ;
  // faces
  int nFE = mesh.FaceEdge.size_loc(iF);
  LOGF << "#edges: nFE=" << nFE << ", { " ;
  for ( int ilE=0; ilE<nFE-1; ++ilE ) {
    LOGF << "(" 
	 << mesh.ok_face_edge(iF,ilE) << ")" 
	 << mesh.face_edge   (iF,ilE)+offset << ", " ;
  }
  LOGF << "(" 
       << mesh.ok_face_edge(iF,nFE-1) << ")" 
       << mesh.face_edge   (iF,nFE-1)+offset << " }" << endl ;
  LOGF << flush ;
  // topological info
  vector<int> rlist, flist, elist, vlist ;
  mesh.get_face_regn( iF, rlist ) ;
  mesh.get_face_face( iF, flist ) ;
  mesh.get_face_edge( iF, elist ) ;
  mesh.get_face_vrtx( iF, vlist ) ;
  LOGF << "---" << endl ;
  LOGF << "#regions:  nFR=" << rlist.size() ; print_list_offset(rlist) ;
  LOGF << "#faces:    nFF=" << flist.size() ; print_list_offset(flist) ;
  LOGF << "#edges:    nFE=" << elist.size() ; print_list_offset(elist) ;
  LOGF << "#vertices: nFV=" << vlist.size() ; print_list_offset(vlist) ;
  // geometrical info
  LOGF << "---" << endl ;
  LOGF << "#center: ( " 
       << mesh.coords_F(iF,0) << ", " 
       << mesh.coords_F(iF,1) << ", " 
       << mesh.coords_F(iF,2) << " )" << endl ; 
  LOGF << "#normal: ( " 
       << mesh.get_nor(iF,0) << ", " 
       << mesh.get_nor(iF,1) << ", " 
       << mesh.get_nor(iF,2) << " )" << endl ; 
  LOGF << "#area:     " << mesh.get_face_measure(iF) << endl ; 
}
// --------------------------------------------------------------------------------------------
// EDGES
// --------------------------------------------------------------------------------------------
void mesh3Dv_printer :: print_all_edges() {
  ofstream LOGF("edge.log") ;
  LOGF << "#of edges: nE=" << mesh.n_edge() << endl ;
  for ( int iE=0; iE<mesh.n_edge(); ++iE ) {
    print_edge( iE, LOGF ) ;
  }
  LOGF.close() ;
}
void mesh3Dv_printer :: print_edge( int iE, ostream & LOGF ) {
  LOGF << "--------------------------------------------------------------" << endl ;
  LOGF << "edge: iE =" << iE+offset << endl ;
  // edges
  int nEV = 2 ;
  LOGF << "#vertices: nVE=" << nEV << ", { " ;
  LOGF << mesh.edge_vrtx(iE,0)+offset << ", " 
       << mesh.edge_vrtx(iE,1)+offset << " }" << endl << flush ;
  // topological info
  vector<int> rlist, flist, elist, vlist ;
  mesh.get_edge_regn( iE, rlist ) ;
  mesh.get_edge_face( iE, flist ) ;
  mesh.get_edge_edge( iE, elist ) ;
  mesh.get_edge_vrtx( iE, vlist ) ;
  LOGF << "---" << endl ;
  LOGF << "#regions:  nER=" << rlist.size() ; print_list_offset(rlist) ;
  LOGF << "#faces:    nEF=" << flist.size() ; print_list_offset(flist) ;
  LOGF << "#edges:    nEE=" << elist.size() ; print_list_offset(elist) ;
  LOGF << "#vertices: nEV=" << vlist.size() ; print_list_offset(vlist) ;
  // geometrical info
  LOGF << "---" << endl ;
  LOGF << "#center: ( " 
       << mesh.coords_E(iE,0) << ", " 
       << mesh.coords_E(iE,1) << ", " 
       << mesh.coords_E(iE,2) << " )" << endl ; 
  LOGF << "#normal: ( " 
       << mesh.get_tng(iE,0) << ", " 
       << mesh.get_tng(iE,1) << ", " 
       << mesh.get_tng(iE,2) << " )" << endl ; 
  LOGF << "#length:     " << mesh.get_edge_measure(iE) << endl ; 
}
// --------------------------------------------------------------------------------------------
// VERTICES
// --------------------------------------------------------------------------------------------
void mesh3Dv_printer :: print_all_vertices() {
  ofstream LOGF("vertex.log") ;
  LOGF << "#of vertices: nV=" << mesh.n_vertex() << endl ;
  for ( int iV=0; iV<mesh.n_vertex(); ++iV ) {
    print_vertex( iV, LOGF ) ;
  }
  LOGF.close() ;
}
void mesh3Dv_printer :: print_vertex( int iV, ostream & LOGF ) {
  LOGF << "--------------------------------------------------------------" << endl ;
  LOGF << "vertex: iV =" << iV+offset << endl ;
  // coordinates
  LOGF << "#coords:= ( " 
       << mesh.V_coords(iV,0) << ", " 
       << mesh.V_coords(iV,1) << ", " 
       << mesh.V_coords(iV,2) << ") " << endl << flush ; 
  // topological info
  vector<int> rlist, flist, elist, vlist ;
  mesh.get_vrtx_regn( iV, rlist ) ;
  mesh.get_vrtx_face( iV, flist ) ;
  mesh.get_vrtx_edge( iV, elist ) ;
  mesh.get_vrtx_vrtx( iV, vlist ) ;
  LOGF << "---" << endl ;
  LOGF << "#regions:  nVR=" << rlist.size() ; print_list_offset(rlist) ;
  LOGF << "#faces:    nVF=" << flist.size() ; print_list_offset(flist) ;
  LOGF << "#edges:    nVE=" << elist.size() ; print_list_offset(elist) ;
  LOGF << "#vertices: nVV=" << vlist.size() ; print_list_offset(vlist) ;
}
