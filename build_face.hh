/* ----------------------------------------------------------------/
This is open source software; you can redistribute it and/or modify
it under the terms of the BSD-3 License. If software is modified
to produce derivative works, such modified software should be
clearly marked, so as not to confuse it with the version available
from LANL. Full text of the BSD-3 License can be found in the
License file of the repository.
/--------------------------------------------------------------- */

struct face_struct {
public:
  vector<int> vlist ; // original vertex list
  vector<int> slist ; // sorted   vertex list
  int iR  ; // region of the face iF  
  int ilF ; // local position of iF inside iR
  face_struct( vector<int> & _vlist, int _iR, int _ilF ) : 
    vlist(_vlist), slist(_vlist), iR(_iR), ilF(_ilF) {
    sort( slist.begin(), slist.end() ) ; 
  }
  ~face_struct() {}
} ;

bool operator<( const face_struct & F0, const face_struct & F1 ) {
  bool retval = F0.slist.size()<F1.slist.size() ;
  if ( F0.slist.size()==F1.slist.size() ) { 
    int i = 0 ;
    while( F0.slist[i]==F1.slist[i] && i<F0.slist.size() ) { ++i ; } 
    if ( i<F0.slist.size() ) { 
      retval = F0.slist[i]<F1.slist[i] ;
    } else {
      retval = F0.iR<F1.iR ;
    }
  }
  return retval ;
}

bool operator==( const face_struct & F0, const face_struct & F1 ) {
  bool retval(false) ;
  if ( F0.slist.size()==F1.slist.size() ) { 
    retval = true ;
    for ( int i=0; i<F0.slist.size() && retval; ++i ) {
      retval &= F0.slist[i]==F1.slist[i] ;
    }
  }
  return retval ;
}

bool operator!=( const face_struct & F0, const face_struct & F1 ) {
  return !( F0==F1 ) ;
}

// INPUT:  regn_flist
// OUTPUT:
// (i)  build face_vlist from regn_flist
// (ii) collect face indices for mesh::RegnFace
void mesh3Dv_builder ::
build_face_vlist( vector< vector<regn_face_struct> > & vec_regn_face,   // output (construction)
		  vector<int>                        & face_vlist,      // output (new)
		  vector<int>                        & regn_flist ) {   // input
  
  vector<int> vlist ;
  vector<face_struct> vec_face_struct ;

  // fill vec_face_struct with instances of face_struct type
  // from all the region faces (any order is OK)
  int nR = regn_flist.back() ;
  int k = 0;
  for ( int iR=0; iR<nR; ++iR ) {         // loop on all the regions
    int nRF = regn_flist[k++] ;           // #of faces of region iR
    for ( int ilF=0; ilF<nRF; ++ilF ) {   // loop on the faces of region iR
      int nFV = regn_flist[k++] ;         // #of vertices of face ilF
      vlist.resize(nFV) ;                 // resize vlist
      for ( int ilV=0; ilV<nFV; ++ilV ) { // loop on the vertices of face ilF
	vlist[ilV] = regn_flist[k++] ;    // put face's vertices into vlist
      }
      vec_face_struct.push_back( face_struct(vlist,iR,ilF) ) ;
    }
  }

  // build an ORDERED set of faces represented by face_struct instances 
  // from all the region faces: internal faces are repetead twice
  sort( vec_face_struct.begin(), vec_face_struct.end() ) ;

  if (0) { // DEBUG!!!   print...
    int nF = vec_face_struct.size() ;
    for ( int k=0; k<nF; ++k ) { 
      //cout << k << vec_face_struct[k].iR  << vec_face_struct[k].ilF ;
      cout << setw(6) << k << " -->> " ;
      int nFV = vec_face_struct[k].vlist.size() ; // #of vertices of face ilF
      for ( int ilV=0; ilV<nFV; ++ilV ) { // loop on the vertices of face ilF
	cout << " " << setw(6) << vec_face_struct[k].slist[ilV] ;
      }
      cout << "<<-- " ; 
      int iR  = vec_face_struct[k].iR  ; 
      int ilF = vec_face_struct[k].ilF ;
      cout << " iR  = " << setw(6) << iR 
	   << " ilF = " << setw(3) << ilF ;
	//VAL(iR) ; VAL(ilF) ; 
      cout << " -->> " ;
      for ( int ilV=0; ilV<nFV; ++ilV ) { // loop on the vertices of face ilF
	cout << " " << setw(6) << vec_face_struct[k].vlist[ilV] ;
      }
      cout << endl ;
    }
  }

  { // make data structure face_vlist & set RegnFace
    // set the first face
    int iF = 0 ; // face index
    int kF = 0 ; // local pointer, be careful, kF==ilF is true only at this step!
    // --store the vertex list
    face_vlist.push_back( vec_face_struct[kF].vlist.size() ) ;
    for ( int ilV=0; ilV<vec_face_struct[kF].vlist.size(); ++ilV ) {
      face_vlist.push_back( vec_face_struct[kF].vlist[ilV] ) ;
    }
    { // --set RegnFace, face index is iF:=0
      int iR  = vec_face_struct[kF].iR  ;
      int ilF = vec_face_struct[kF].ilF ;
      vec_regn_face[iR].push_back( regn_face_struct(iF,ilF,bool(true)) ) ;
    }
    // set the remaining faces
    for ( int il=1; il<vec_face_struct.size(); ++il ) {
      if ( vec_face_struct[il]!=vec_face_struct[kF] ) {
	// found a new face
	iF++    ; // update the face counter
	kF = il ; // reset  the face pointer
	// --store the vertex list
	face_vlist.push_back( vec_face_struct[kF].vlist.size() ) ;
	for ( int ilV=0; ilV<vec_face_struct[kF].vlist.size(); ++ilV ) {
	  face_vlist.push_back( vec_face_struct[kF].vlist[ilV] ) ;
	}
      }
      // --set RegnFace, face index is iF, ilF is the other instance of the face iF
      // --take the orientation of the face iF in the first region
      int  iR   = vec_face_struct[il].iR  ;
      int  ilF  = vec_face_struct[il].ilF ;
      bool bval = kF==il                   ;
      vec_regn_face[iR].push_back( regn_face_struct(iF,ilF,bval) ) ;
    }
    // finally, set the number of mesh face
    face_vlist.push_back( iF+1 ) ;
  }
}

// ------- REMARKS ------- 
//
// In this implementation, the face orientation coincides with the
// orientation of the face in the first of the two regions that this
// face belongs to.
//
// This region-based face orientation is given by the sequence of
// vertices that is input to the builder for the face in the region.
//
// The first region is determined by sorting all structures of type
// face_struct. So, the first region is the one between two with the
// smallest identifier.
//
// Implicitly, we assume that the face orientation of this region is
// correct, while in the other region the face may be correct
// (opposite orientation, reversed vertex sequence) or wrong (same
// orientation, same vertex sequence).
//
// The other face does not influence the regn-face based construction.
//
// NOTE that the oldest implementation of mesh3D_writer that outputs
// meshes in REGN_FACE format was bugged for this reason, as it did
// not distinguish between the orientation in the two regions.
//
// Nonetheless, everything worked properly because the first region
// was always written in the correct way, and the second (wrong) was
// never used to build mesh faces.
