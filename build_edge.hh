/* ----------------------------------------------------------------/
This is open source software; you can redistribute it and/or modify
it under the terms of the BSD-3 License. If software is modified
to produce derivative works, such modified software should be
clearly marked, so as not to confuse it with the version available
from LANL. Full text of the BSD-3 License can be found in the
License file of the repository.
/---------------------------------------------------------------- */

struct edge_struct {
public:
  int V0, V1  ;
  int iF, ilF ; 
  edge_struct( int _V0, int _V1, int _iF, int _ilF ) : V0(_V0), V1(_V1), iF(_iF), ilF(_ilF) {
    if ( V0 > V1 ) {
      V0   = _V1 ; 
      V1   = _V0 ;
      ilF *= -1 ;
    }
  }
  ~edge_struct() {}
} ;

bool operator<( const edge_struct & E0, const edge_struct & E1 ) {
  return 
    ( E0.V0 <  E1.V0                                  ) ||
    ( E0.V0 == E1.V0 && E0.V1 <  E1.V1                ) ||
    ( E0.V0 == E1.V0 && E0.V1 == E1.V1 && E0.iF<E1.iF ) ;
}

bool operator==( const edge_struct & E0, const edge_struct & E1 ) {
  return E0.V0==E1.V0 && E0.V1==E1.V1 ; // may belong to different faces
}

bool operator!=( const edge_struct & E0, const edge_struct & E1 ) {
  return !( E0==E1 ) ;
}

void mesh3Dv_builder ::
build_face_elist( vector< vector<face_edge_struct> > & vec_face_edge,
		  vector<int>                        & edge_vlist,
		  vector<int>                        & face_vlist ) {
  
  vector<edge_struct> tmp_edge_struct ;
  int nF = face_vlist.back() ;
  
  int kV = 0;
  for ( int iF=0; iF<nF; ++iF ) {
    int nFV = face_vlist[kV++] ;
    for ( int ilV=0; ilV<nFV-1; ++ilV ) {
      int kV0 = face_vlist[kV] ;
      int kV1 = face_vlist[kV+1] ;
      tmp_edge_struct.push_back(edge_struct(kV0,kV1,iF,ilV+1)) ; // local pos. starts from 1 to store
      kV++ ;                                                     // the edge orientation in the face
    }
    int kV0 = face_vlist[kV] ;         // last  vertex of the face
    int kV1 = face_vlist[kV-(nFV-1)] ; // first vertex of the face: kV-(nFV-1)==0 !!!
    tmp_edge_struct.push_back(edge_struct(kV0,kV1,iF,nFV)) ;
    kV++ ;
  }
  
  sort( tmp_edge_struct.begin(), tmp_edge_struct.end() ) ;
  
  { // set up datasets
    // first instance
    int kE = 0 ; // local pointer
    int iE = 0 ; // number of edges 
    edge_vlist.push_back(2) ; // each edge has two vertices
    edge_vlist.push_back( tmp_edge_struct[kE].V0 ) ;
    edge_vlist.push_back( tmp_edge_struct[kE].V1 ) ;

    int iF  = tmp_edge_struct[iE].iF  ;
    int ilF = tmp_edge_struct[iE].ilF ;

    // be careful: ilF starts from 1 to keep sign (+/-) info!
    bool bval = ilF>0 ? true : false ;
    vec_face_edge[iF].push_back( face_edge_struct(iE,abs(ilF)-1,bval) ) ;
    
    for ( int il=1; il<tmp_edge_struct.size(); ++il ) { // il = running on all the faces
      if ( tmp_edge_struct[il]!=tmp_edge_struct[kE] ) {
	kE = il ; // update pointer to new face
	iE++ ;    // update counter
	edge_vlist.push_back(2) ; // each edge has two vertices
	edge_vlist.push_back( tmp_edge_struct[il].V0 ) ;
	edge_vlist.push_back( tmp_edge_struct[il].V1 ) ;
      }
      int iF  = tmp_edge_struct[il].iF   ;
      int ilF = tmp_edge_struct[il].ilF  ;
      bool bval = ilF>0 ? true : false ;
      vec_face_edge[iF].push_back( face_edge_struct(iE,abs(ilF)-1,bval) ) ;
    }
    edge_vlist.push_back( ++iE ) ;
  }
}
