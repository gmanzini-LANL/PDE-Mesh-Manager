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
    mesh_dir  = string("./Mesh-3D/Meshes/STANDARD-SET/Cubic-Cells/") ;
    mesh_name = string("Cubic cells") ;
    vec_fname[0] = string("gcube_1x1x1") ;
    vec_fname[1] = string("gcube_2x2x2") ;
    vec_fname[2] = string("gcube_4x4x4") ;
    vec_fname[3] = string("gcube_8x8x8") ;
    vec_fname[4] = string("gcube_16x16x16") ;
    vec_fname[5] = string("gcube_32x32x32") ;
    vec_fname[6] = string("gcube_64x64x64") ;
  }
} ;

// --------------------------------------------------------------------------------------------
// previous class   : Mesh_M3
// previous version : Layered/Dual-1b
// current  version : STANDARD-SET/Prysmatic-Cells
// mesh index: 2
class Mesh_PrysmaticCells : public MeshDirsBase {
public:
  Mesh_PrysmaticCells() : MeshDirsBase(8) { setup() ; }
  ~Mesh_PrysmaticCells() {}
  virtual void setup() { // mesh M3
    fmt_flag  = REGN_FACE_FMT ;
    mesh_dir  = string("./Mesh-3D/Meshes/STANDARD-SET/Prysmatic-Cells/") ;
    mesh_name = string("Prysmatic cells with polygonal basis") ;
    vec_fname[0] = string("gdual_5x5x5") ;
    vec_fname[1] = string("gdual_10x10x10") ;
    vec_fname[2] = string("gdual_15x15x15") ;
    vec_fname[3] = string("gdual_20x20x20") ;
    vec_fname[4] = string("gdual_25x25x25") ;
    vec_fname[5] = string("gdual_30x30x30") ;
    vec_fname[6] = string("gdual_35x35x35") ;
    vec_fname[7] = string("gdual_40x40x40") ;
  }
} ;

// --------------------------------------------------------------------------------------------
// previous class   : Mesh_M4
// previous version : Hexa-test-2
// current  version : STANDARD_SET/Random-Hexahedra
// mesh index: 3
class Mesh_RandomHexahedra : public MeshDirsBase {
public:
  Mesh_RandomHexahedra() : MeshDirsBase(5) { setup() ; }
  ~Mesh_RandomHexahedra() {}
  virtual void setup() { // mesh M4
    fmt_flag  = REGN_FACE_FMT ;
    mesh_dir  = string("./Mesh-3D/Meshes/STANDARD-SET/Random-Hexahedra/") ;
    mesh_name = string("Random Hexahedra") ;
    vec_fname[0] = string("gcube.1") ;
    vec_fname[1] = string("gcube.2") ;
    vec_fname[2] = string("gcube.3") ;
    vec_fname[3] = string("gcube.4") ;
    vec_fname[4] = string("gcube.5") ;
  }
} ;

class Mesh_Cube_Tetrahedra_0 : public MeshDirsBase {
public:
  Mesh_Cube_Tetrahedra_0() : MeshDirsBase(6) { setup() ; }
  ~Mesh_Cube_Tetrahedra_0() {}
  virtual void setup() { // mesh
    fmt_flag = TETGEN_FMT ;
    mesh_dir = string("./Mesh-3D/Meshes/STANDARD-SET/Tetgen-Cube-0/") ;
    mesh_name = string("Cube - Tetrahedra") ;
    vec_fname[0] = string("cube.1") ;
    vec_fname[1] = string("cube.2") ;
    vec_fname[2] = string("cube.3") ;
    vec_fname[3] = string("cube.4") ;
    vec_fname[4] = string("cube.5") ;
    vec_fname[5] = string("cube.6") ;
    for ( int i=0; i<vec_offset.size(); ++i ) { vec_offset[i] = 1 ; }
  }
} ;

class Mesh_Cube_Tetrahedra_1 : public MeshDirsBase {
public:
  Mesh_Cube_Tetrahedra_1() : MeshDirsBase(5) { setup() ; }
  ~Mesh_Cube_Tetrahedra_1() {}
  virtual void setup() { // mesh 
    fmt_flag = TETGEN_FMT ;
    mesh_dir = string("./Mesh-3D/Meshes/STANDARD-SET/Tetgen-Cube-1/") ;
    mesh_name = string("Cube - Tetrahedra") ;
    vec_fname[0] = string("cube.1") ;
    vec_fname[1] = string("cube.2") ;
    vec_fname[2] = string("cube.3") ;
    vec_fname[3] = string("cube.4") ;
    vec_fname[4] = string("cube.5") ;
    for ( int i=0; i<vec_offset.size(); ++i ) { vec_offset[i] = 1 ; }
  }
} ;

// --------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------
//---- set mesh dirs
// --------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------
MeshDirsBase & set_mesh_dirs( int mesh_flag ) {
  MeshDirsBase * p_mdir ;
  switch ( mesh_flag ) {
  case  0 : p_mdir = new Mesh_CubicCells        ; break ;
  case  1 : p_mdir = new Mesh_PrysmaticCells    ; break ;
  case  2 : p_mdir = new Mesh_RandomHexahedra   ; break ;
  case  3 : p_mdir = new Mesh_Cube_Tetrahedra_0 ; break ;
  case  4 : p_mdir = new Mesh_Cube_Tetrahedra_1 ; break ;
  default: assert( bool(false) ) ;
  }
  return *p_mdir ;
}

#endif // end of  _MESH_DIRS
