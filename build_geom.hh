/* ----------------------------------------------------------------/
This is open source software; you can redistribute it and/or modify
it under the terms of the BSD-3 License. If software is modified
to produce derivative works, such modified software should be
clearly marked, so as not to confuse it with the version available
from LANL. Full text of the BSD-3 License can be found in the
License file of the repository.
/----------------------------------------------------------------*/

// build geometric quantities
void mesh3Dv_builder :: setup_geom_factors() {
  // setup of mesh arrays
  // regions:
  int nR = mesh.n_region() ;
  mesh.R_volume.setup(nR,1) ;
  mesh.R_coords.setup(nR,3) ;
  // faces:
  int nF = mesh.n_face() ;
  mesh.F_area  .setup(nF,1) ;
  mesh.F_nor   .setup(nF,3) ;
  mesh.F_coords.setup(nF,3) ;
  // ----------------------
#if 1
  set_regn_geom_factors() ;
  set_face_geom_factors() ;
#else
    cerr << ">>>>>>>>>>>REMOVED GEOMETRIC FACTORS!!!<<<<<<<<<<<<<<" << endl << flush ;
#endif
}
void mesh3Dv_builder :: set_regn_geom_factors() {
  MSG("begin set_regn_geom_factors"<<endl<<flush) ;

  int nR   = mesh.n_region() ;
  int rule = 1 ;

#if USE_TETGEN
  QuadRegion quad(mesh,rule) ;
#endif

  double xV[4], yV[4], zV[4] ;

  for ( int iR=0; iR<nR; ++iR ) {
    // reset local vars
    double vol_R(0.), xR(0.), yR(0.), zR(0.) ;
    // set region iR
    vector<int> vlist ;
    mesh.get_regn_vrtx( iR, vlist ) ;
    int nRV = vlist.size() ;
    if ( nRV==4 ) { // R is a tetrahedron
      // get the tetrahedron vertices
      for ( int ilV=0; ilV<nRV; ++ilV ) {
	int iV = vlist[ilV] ;
	xV[ilV] = mesh.coords_V(iV,0) ;
	yV[ilV] = mesh.coords_V(iV,1) ;
	zV[ilV] = mesh.coords_V(iV,2) ;
	xR += xV[ilV] ;
	yR += yV[ilV] ;
	zR += zV[ilV] ;
      }
      // get the volume
      vol_R = tetra_volume( xV, yV, zV ) ;
      xR /= double(nRV) ;
      yR /= double(nRV) ;
      zR /= double(nRV) ;
    } else { // R is NOT a tetrahedron
#if USE_TETGEN
      quad.setup(iR) ;
      quad.region_intg_values( vol_R, xR, yR, zR ) ;
      xR /= vol_R ; 
      yR /= vol_R ;
      zR /= vol_R ;
#else // working only for convex cells with all convex faces
      // -----------------------------------------------------
      new_region_intg_values( iR, vlist, vol_R, xR, yR, zR ) ;
      // -----------------------------------------------------
#endif
    }
    // set geometric values
    mesh.R_volume(iR,0) = vol_R ;
    mesh.R_coords(iR,0) = xR ;
    mesh.R_coords(iR,1) = yR ;
    mesh.R_coords(iR,2) = zR ;
    //VAL(iR) ; PRT(vol_R) ;
  }
  MSG("end-->set_regn_geom_factors"<<endl<<flush) ;
}
void mesh3Dv_builder :: set_face_geom_factors() {
  MSG("begin set_face_geom_factors"<<endl<<flush) ;

  int nF   = mesh.n_face() ;
  int rule = 1 ;
  double face_err = 0. ;

  for ( int iF=0; iF<nF; ++iF ) {

    vector<int> vlist ;
    mesh.get_face_vrtx( iF, vlist ) ;

    int nFV = vlist.size() ;

    double xc[3] = { 0., 0., 0. } ; // barycenter of the face
    double v0(0.), v1(0.), v2(0.), ar(0.) ;
    
    if ( nFV==3 ) { // F is a triangle

      double x[2], y[2], z[2] ;
      int iV0 = vlist[0] ;
      for ( int il=0; il<2; ++il ) {
	int iV = vlist[il+1] ;
	x[il] = mesh.coords_V(iV,0)-mesh.coords_V(iV0,0) ;
	y[il] = mesh.coords_V(iV,1)-mesh.coords_V(iV0,1) ;
	z[il] = mesh.coords_V(iV,2)-mesh.coords_V(iV0,2) ;
      }
      v0 = y[0]*z[1]-y[1]*z[0] ;
      v1 = x[1]*z[0]-x[0]*z[1] ;
      v2 = x[0]*y[1]-x[1]*y[0] ;
      ar = sqrt( v0*v0+v1*v1+v2*v2 ) ;
      
      // be careful that ar is twice the triangle area!
      for ( int ilV=0; ilV<nFV; ++ilV ) {
	int iV = vlist[ilV] ;
	for ( int s=0; s<3; ++s ) {
	  xc[s] += mesh.coords_V(iV,s) * ar / 3. ;
	}
      }

    } else { // F is a planar polygon 

      // get sub-triangle vertices
      double xF = mesh.ari_coords_F( iF, 0 ) ;
      double yF = mesh.ari_coords_F( iF, 1 ) ;
      double zF = mesh.ari_coords_F( iF, 2 ) ;

      for ( int ilV=0; ilV<nFV; ++ilV ) {
	int iV  = vlist[ ilV ] ;
	int iV1 = vlist[ (ilV+1)%nFV ] ;
	
	double x[2] = { mesh.coords_V(iV,0)-xF, mesh.coords_V(iV1,0)-xF } ;
	double y[2] = { mesh.coords_V(iV,1)-yF, mesh.coords_V(iV1,1)-yF } ;
	double z[2] = { mesh.coords_V(iV,2)-zF, mesh.coords_V(iV1,2)-zF } ;

	double tmp_v0 = y[0]*z[1]-y[1]*z[0] ;
	double tmp_v1 = x[1]*z[0]-x[0]*z[1] ;
	double tmp_v2 = x[0]*y[1]-x[1]*y[0] ;
	double tmp_ar = sqrt( tmp_v0*tmp_v0+tmp_v1*tmp_v1+tmp_v2*tmp_v2 ) ;
	
	v0 += tmp_v0 ;
	v1 += tmp_v1 ;
	v2 += tmp_v2 ;
	ar += tmp_ar ;

	// be careful that: tmp_ar is twice the triangle area!
	// (the coefficient should be sub_area_triangle/6)
	for ( int s=0; s<3; ++s ) {
	  xc[s] += ( mesh.ari_coords_F(iF,s) + 
		     mesh.coords_V(iV,s) + mesh.coords_V(iV1,s) )*tmp_ar/3 ;
	}
      }
    }
    // --------------------------
    mesh.F_area  (iF,0) = ar/2. ;
    // --------------------------
    mesh.F_nor   (iF,0) = v0/ar ;
    mesh.F_nor   (iF,1) = v1/ar ;
    mesh.F_nor   (iF,2) = v2/ar ;
    // -----------------------------
    mesh.F_coords(iF,0) = xc[0]/ar ;
    mesh.F_coords(iF,1) = xc[1]/ar ;
    mesh.F_coords(iF,2) = xc[2]/ar ;
    // -----------------------------
    face_err = max( face_err, abs(sqrt( pow(v0/ar,2)+pow(v1/ar,2)+pow(v2/ar,2) )-1.) ) ;
    //assert( face_err<1.e-14 ) ;
  } // end of --> for ( int iF=0; iF<nF; ++iF ) {...
  
  PRT(face_err) ;
  MSG("end-->set_face_geom_factors"<<endl<<flush) ;
}

void mesh3Dv_builder :: change_bbox( double new_xmin, double new_ymin, double new_zmin,     // min vertex
				     double new_xmax, double new_ymax, double new_zmax ) {  // max vertex
  assert( new_xmax>new_xmin ) ;
  assert( new_ymax>new_ymin ) ;
  assert( new_zmax>new_zmin ) ;
  assert( DIM==3 ) ;

  double new_max[DIM] = { new_xmax, new_ymax, new_zmax } ;
  double new_min[DIM] = { new_xmin, new_ymin, new_zmin } ;

  double vol_scal = 1. ;
  for ( int s=0; s<DIM; ++s ) {
    vol_scal *= new_max[s]-new_min[s] ;
  }

  // set bounding box
  if ( !mesh.bb_status ) { mesh.eval_bbox() ; }

  double bb_den[DIM] ;
  for ( int s=0; s<DIM; ++s ) {
    bb_den[s] = mesh.bb_max[s]-mesh.bb_min[s] ;
  }

  // set new vertex coords
  for ( int iV=0; iV<mesh.nV; ++iV ) {
    for ( int s=0; s<DIM; ++s ) {
      mesh.V_coords(iV,s) = ( ( mesh.V_coords(iV,s)-mesh.bb_min[s] ) * new_max[s] + ( mesh.bb_max[s]-mesh.V_coords(iV,s) ) * new_min[s] ) / bb_den[s] ;
    }
  } 
  
  // set new region coords
  for ( int iR=0; iR<mesh.nR; ++iR ) {
    for ( int s=0; s<DIM; ++s ) {
      mesh.R_coords(iR,s) = ( ( mesh.R_coords(iR,s)-mesh.bb_min[s] ) * new_max[s] + ( mesh.bb_max[s]-mesh.R_coords(iR,s) ) * new_min[s] ) / bb_den[s] ;
    }
    mesh.R_volume(iR,0) *= vol_scal ;
  }
  
  // set face's geometric factors
  set_face_geom_factors() ;
}

void mesh3Dv_builder :: new_region_intg_values( int iR, vector<int> & vlist, double & vol_R, double & xR, double & yR, double & zR ) {
  assert( 0<=iR && iR<mesh.n_region() ) ;

  // calculate an internal point (xR0,yR0,zR0)
  double xR0(0.), yR0(0.), zR0(0.) ;
  int nRV = vlist.size() ;
  for ( int ilV=0; ilV<nRV; ++ilV ) {
    int iV = vlist[ilV] ;
    xR0 += mesh.coords_V(iV,0) ;
    yR0 += mesh.coords_V(iV,1) ;
    zR0 += mesh.coords_V(iV,2) ;
  }
  xR0 /= double(nRV) ;
  yR0 /= double(nRV) ;
  zR0 /= double(nRV) ;

  // set flist
  vector<int> flist ;
  mesh.get_regn_face( iR, flist ) ;
  int nRF = flist.size() ;

  // calculate an internal point of F (xF,yF,zF)
  vector<double> xF(nRF), yF(nRF), zF(nRF) ; 
  for ( int ilF=0; ilF<nRF; ++ilF ) {
    int iF = flist[ilF] ;
    vector<int> F_vlist ;
    mesh.get_face_vrtx( iF, F_vlist ) ;
    int nFV = F_vlist.size() ;
    for ( int ilV=0; ilV<nFV; ++ilV ) {
      int iV = F_vlist[ilV] ;
      xF[ilF] += mesh.coords_V(iV,0) ;
      yF[ilF] += mesh.coords_V(iV,1) ;
      zF[ilF] += mesh.coords_V(iV,2) ;
    }
    xF[ilF] /= double(nFV) ;
    yF[ilF] /= double(nFV) ;
    zF[ilF] /= double(nFV) ;
  }

  // loop on internal diamonds and accumulate
  xR = yR = zR = 0. ;
  double vR = 0. ;
  double xV[4], yV[4], zV[4] ;
  for ( int ilF=0; ilF<nRF; ++ilF ) {
    int iF = flist[ilF] ;
    vector<int> F_vlist ;
    mesh.get_face_vrtx( iF, F_vlist ) ;
    int nFV = F_vlist.size() ;
    for ( int ilV=0; ilV<nFV; ++ilV ) {
      for ( int s=0; s<2; ++s ) {
	int iV = F_vlist[ (ilV+s)%nFV ] ;
	xV[s] = mesh.coords_V(iV,0) ;
	yV[s] = mesh.coords_V(iV,1) ;
	zV[s] = mesh.coords_V(iV,2) ;
      }
      xV[2] = xF[ilF] ;
      yV[2] = yF[ilF] ;
      zV[2] = zF[ilF] ;
      // --------------
      xV[3] = xR0 ;
      yV[3] = yR0 ;
      zV[3] = zR0 ;
      // --------------
      double vol_tetra = tetra_volume( xV, yV, zV ) ;
      vR += vol_tetra ;
      xR += ( xV[0]+xV[1]+xV[2]+xV[3] )/4. * vol_tetra ;
      yR += ( yV[0]+yV[1]+yV[2]+yV[3] )/4. * vol_tetra ;
      zR += ( zV[0]+zV[1]+zV[2]+zV[3] )/4. * vol_tetra ;
    }
  }

  // final setting
  vol_R = vR ;
  xR /= vR ;
  yR /= vR ;
  zR /= vR ;
}
