/* ----------------------------------------------------------------/
This is open source software; you can redistribute it and/or modify
it under the terms of the BSD-3 License. If software is modified
to produce derivative works, such modified software should be
clearly marked, so as not to confuse it with the version available
from LANL. Full text of the BSD-3 License can be found in the
License file of the repository.
/---------------------------------------------------------------- */

#ifndef _MESH_DIRS
#define _MESH_DIRS

const int TETGEN_FMT    = 0 ;
const int REGN_FACE_FMT = 1 ;

class MeshDirsBase {
protected:
  int    fmt_flag  ;
  string mesh_dir  ;
  string mesh_name ;
  int    vec_size  ;
  vector<int>    vec_offset ;
  vector<string> vec_fname  ;
public:
  MeshDirsBase( int _vec_size ) : vec_size (_vec_size), 
				  fmt_flag (-1), 
				  mesh_dir ("<undefined>"), 
				  mesh_name("<undefined>") {
    vec_offset.resize( vec_size ) ;
    vec_fname .resize( vec_size ) ;
    for ( int i=0; i<vec_size; ++i ) { 
      vec_offset[i] = 0 ; 
      vec_fname [i] = string("<undefined>") ; 
    }
  } ;
  ~MeshDirsBase() {} ;
  virtual void setup() = 0 ;
  int get_offset( int ifile ) {
    PRT(ifile) ;
    assert( 0<=ifile && ifile<vec_size ) ;
    return vec_offset[ifile] ;
  }
  string get_fname( int ifile ) {
    assert( 0<=ifile && ifile<vec_size ) ;
    string retval = mesh_dir ;
    switch ( fmt_flag ) {
    case TETGEN_FMT:    retval += string("TG_fmt/") ; break ;
    case REGN_FACE_FMT: retval += string("RF_fmt/") ; break ;
    default: assert(false) ;
    }
    retval += vec_fname[ifile] ;
    return retval ;
  }
  int    get_fmt_flag () { return fmt_flag  ; }
  string get_mesh_dir () { return mesh_dir  ; }
  string get_mesh_name() { return mesh_name ; } 
  string get_mesh_file( int ifile ) { return vec_fname[ifile] ; }
} ;
// --------------------------------------------------------------------------------------------
// mesh index 0
// --------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------
// previous class  : Mesh_M1
// previous version: Layered/Quad-2/
// current  version: STANDARD_SET/Cubic-Cells
// mesh index: 1
class Mesh_CubicCells : public MeshDirsBase {
public:
  Mesh_CubicCells() : MeshDirsBase(7) { setup() ; }
  ~Mesh_CubicCells() {}
  virtual void setup() { // mesh M1
    fmt_flag  = REGN_FACE_FMT ;
    mesh_dir  = string("./MeshDataSets/Mesh-3D/STANDARD-SET/Cubic-Cells/") ;
    mesh_name = string("Cubic cells") ;
    vec_fname[0] = string("gcube_1x1x1") ;      // TESTING VERSION
    vec_fname[1] = string("gcube_2x2x2") ;
    vec_fname[2] = string("gcube_4x4x4") ;
    vec_fname[3] = string("gcube_8x8x8") ;
    vec_fname[4] = string("gcube_16x16x16") ;
    vec_fname[5] = string("gcube_32x32x32") ;
    vec_fname[6] = string("gcube_64x64x64") ;
  }
} ; 
// --------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------
class Mesh_DoubleCubes_RandomHex : public MeshDirsBase {
public:
  Mesh_DoubleCubes_RandomHex() : MeshDirsBase(7) { setup() ; }
  ~Mesh_DoubleCubes_RandomHex() {}
  virtual void setup() { // mesh M4
    fmt_flag  = REGN_FACE_FMT ;
    mesh_dir  = string("./MeshDataSets/Mesh-3D/STANDARD-SET/DoubleCubes/") ;
    mesh_name = string("Double Cubes - Random Hexahedra") ;
    vec_fname[0] = string("gdouble_cube.1") ;
    vec_fname[1] = string("gdouble_cube.2") ;
    vec_fname[2] = string("gdouble_cube.3") ;
    vec_fname[3] = string("gdouble_cube.4") ;
    vec_fname[4] = string("gdouble_cube.5") ;
    vec_fname[5] = string("gdouble_cube.6") ;
    vec_fname[6] = string("gdouble_cube.7") ;
  }
} ;
class Mesh_TiltedCubicCells_Y_Z : public MeshDirsBase {
public:
  Mesh_TiltedCubicCells_Y_Z() : MeshDirsBase(7) { setup() ; }
  ~Mesh_TiltedCubicCells_Y_Z() {}
  virtual void setup() { // mesh M1
    fmt_flag  = REGN_FACE_FMT ;
    mesh_dir  = string("./MeshDataSets/Mesh-3D/DOUBLE-PERM/TiltedCubicMesh-Y-Z/") ;
    mesh_name = string("Tilted Cubic cells along Y planes") ;
    vec_fname[0] = string("tgcube_1x1x1") ;      // TESTING VERSION
    vec_fname[1] = string("tgcube_2x2x2") ;
    vec_fname[2] = string("tgcube_4x4x4") ;
    vec_fname[3] = string("tgcube_8x8x8") ;
    vec_fname[4] = string("tgcube_16x16x16") ;
    vec_fname[5] = string("tgcube_32x32x32") ;
    vec_fname[6] = string("tgcube_64x64x64") ;
  }
} ; 
class Mesh_TiltedCubicCells_Y_ZX : public MeshDirsBase {
public:
  Mesh_TiltedCubicCells_Y_ZX() : MeshDirsBase(7) { setup() ; }
  ~Mesh_TiltedCubicCells_Y_ZX() {}
  virtual void setup() { // mesh M1
    fmt_flag  = REGN_FACE_FMT ;
    mesh_dir  = string("./MeshDataSets/Mesh-3D/DOUBLE-PERM/TiltedCubicMesh-Y-Z/") ;
    mesh_name = string("Tilted Cubic cells along Y planes") ;
    vec_fname[0] = string("tgcube_1x1x1") ;      // TESTING VERSION
    vec_fname[1] = string("tgcube_2x2x2") ;
    vec_fname[2] = string("tgcube_4x4x4") ;
    vec_fname[3] = string("tgcube_8x8x8") ;
    vec_fname[4] = string("tgcube_16x16x16") ;
    vec_fname[5] = string("tgcube_32x32x32") ;
    vec_fname[6] = string("tgcube_64x64x64") ;
  }
} ; 
class Mesh_TiltedCubicCells_Y_ZXr : public MeshDirsBase {
public:
  Mesh_TiltedCubicCells_Y_ZXr() : MeshDirsBase(7) { setup() ; }
  ~Mesh_TiltedCubicCells_Y_ZXr() {}
  virtual void setup() { // mesh M1
    fmt_flag  = REGN_FACE_FMT ;
    mesh_dir  = string("./MeshDataSets/Mesh-3D/DOUBLE-PERM/TiltedCubicMesh-Y-Z/") ;
    mesh_name = string("Tilted Cubic cells along Y planes") ;
    vec_fname[0] = string("tgcube_1x1x1") ;      // TESTING VERSION
    vec_fname[1] = string("tgcube_2x2x2") ;
    vec_fname[2] = string("tgcube_4x4x4") ;
    vec_fname[3] = string("tgcube_8x8x8") ;
    vec_fname[4] = string("tgcube_16x16x16") ;
    vec_fname[5] = string("tgcube_32x32x32") ;
    vec_fname[6] = string("tgcube_64x64x64") ;
  }
} ; 

class Mesh_DoubleCubeRandomHex_0 : public MeshDirsBase {
public:
  Mesh_DoubleCubeRandomHex_0() : MeshDirsBase(6) { setup() ; }
  ~Mesh_DoubleCubeRandomHex_0() {}
  virtual void setup() { // mesh M1
    fmt_flag  = REGN_FACE_FMT ;
    mesh_dir  = string("./MeshDataSets/Mesh-3D/DOUBLE-PERM/DoubleRandomHex/") ;
    mesh_name = string("Double Cubes - Random Hexahedra - 1") ;
    vec_fname[0] = string("gcube_1x1x1") ;      // TESTING VERSION
    vec_fname[1] = string("gcube_2x2x2") ;
    vec_fname[2] = string("gcube_4x4x4") ;
    vec_fname[3] = string("gcube_8x8x8") ;
    vec_fname[4] = string("gcube_16x16x16") ;
    vec_fname[5] = string("gcube_32x32x32") ;
  }
};

class Mesh_DoubleCubeRandomHex_1 : public MeshDirsBase {
public:
  Mesh_DoubleCubeRandomHex_1() : MeshDirsBase(7) { setup() ; }
  ~Mesh_DoubleCubeRandomHex_1() {}
  virtual void setup() { // mesh M1
    fmt_flag  = REGN_FACE_FMT ;
    mesh_dir  = string("./MeshDataSets/Mesh-3D/DOUBLE-PERM/DoubleRandomHex/") ;
    mesh_name = string("Double Cubes - Random Hexahedra - 1") ;
    vec_fname[0] = string("gcube_2x2x2") ;      // TESTING VERSION
    vec_fname[1] = string("gcube_4x4x4") ;
    vec_fname[2] = string("gcube_8x8x8") ;
    vec_fname[3] = string("gcube_12x12x12") ;
    vec_fname[4] = string("gcube_16x16x16") ;
    vec_fname[5] = string("gcube_20x20x20") ;
    vec_fname[6] = string("gcube_24x24x24") ;
  }
};

// --------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------
//---- set mesh dirs
// --------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------
MeshDirsBase & set_mesh_dirs( int mesh_flag ) {
  MeshDirsBase * p_mdir ;
  switch ( mesh_flag ) {
  case  1 : p_mdir = new Mesh_CubicCells             ; break ;
  case 11 : p_mdir = new Mesh_DoubleCubes_RandomHex  ; break ;
  case 12 : p_mdir = new Mesh_TiltedCubicCells_Y_Z   ; break ;
  case 13 : p_mdir = new Mesh_TiltedCubicCells_Y_ZX  ; break ;
  case 14 : p_mdir = new Mesh_TiltedCubicCells_Y_ZXr ; break ;
  case 15 : p_mdir = new Mesh_DoubleCubeRandomHex_1  ; break ;
  default: assert( bool(false) ) ;
  }
  return *p_mdir ;
}

#endif // end of  _MESH_DIRS
