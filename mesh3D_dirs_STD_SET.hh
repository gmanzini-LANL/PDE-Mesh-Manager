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
class Mesh_Cube_Tetrahedra_2 : public MeshDirsBase {
public:
  Mesh_Cube_Tetrahedra_2() : MeshDirsBase(6) { setup() ; }
  ~Mesh_Cube_Tetrahedra_2() {}
  virtual void setup() { // mesh
    fmt_flag = TETGEN_FMT ;
    mesh_dir = string("./MeshDataSets/Mesh-3D/STANDARD-SET/CUBE-test-0/") ;
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
    mesh_dir = string("./MeshDataSets/Mesh-3D/STANDARD-SET/Tetgen-Cube-1/") ;
    mesh_name = string("Cube - Tetrahedra") ;
    vec_fname[0] = string("cube.1") ;
    vec_fname[1] = string("cube.2") ;
    vec_fname[2] = string("cube.3") ;
    vec_fname[3] = string("cube.4") ;
    vec_fname[4] = string("cube.5") ;
    for ( int i=0; i<vec_offset.size(); ++i ) { vec_offset[i] = 1 ; }
  }
} ;
class Mesh_Cube_Tetrahedra_0 : public MeshDirsBase {
public:
  Mesh_Cube_Tetrahedra_0() : MeshDirsBase(6) { setup() ; }
  ~Mesh_Cube_Tetrahedra_0() {}
  virtual void setup() { // mesh
    fmt_flag = TETGEN_FMT ;
    mesh_dir = string("./MeshDataSets/Mesh-3D/STANDARD-SET/Tetgen-Cube-0/") ;
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
// --------------------------------------------------------------------------------------------
// previous class  : Mesh_M1
// previous version: Layered/Quad-2/
// current  version: STANDARD_SET/Cubic-Cells
// mesh index: 1
#if 0
class Mesh_CubicCells : public MeshDirsBase {
public:
  Mesh_CubicCells() : MeshDirsBase(6) { setup() ; }
  ~Mesh_CubicCells() {}
  virtual void setup() { // mesh M1
    fmt_flag  = REGN_FACE_FMT ;
    mesh_dir  = string("./MeshDataSets/Mesh-3D/STANDARD-SET/Cubic-Cells/") ;
    mesh_name = string("Cubic cells") ;
    vec_fname[0] = string("gcube_2x2x2") ;
    vec_fname[1] = string("gcube_4x4x4") ;
    vec_fname[2] = string("gcube_8x8x8") ;
    vec_fname[3] = string("gcube_16x16x16") ;
    vec_fname[4] = string("gcube_32x32x32") ;
    vec_fname[5] = string("gcube_64x64x64") ;
  }
} ; 
#endif
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
// previous class   : Mesh_M3
// previous version : Layered/Dual-1b
// current  version : STANDARD-SET/Prysmatic-Cells
// mesh index: 2
// no hanging nodes (they have been removed from the 2D base mesh)
class Mesh_PrysmaticCells : public MeshDirsBase {
public:
  Mesh_PrysmaticCells() : MeshDirsBase(9) { setup() ; }
  ~Mesh_PrysmaticCells() {}
  virtual void setup() { // mesh M3
    fmt_flag  = REGN_FACE_FMT ;
    mesh_dir  = string("./MeshDataSets/Mesh-3D/STANDARD-SET/Prysmatic-Cells/") ;
    mesh_name = string("Prysmatic cells with polygonal basis - no hanging nodes") ;
    vec_fname[0] = string("vprism_2x2x2") ;
    vec_fname[1] = string("vprism_4x4x4") ;
    vec_fname[2] = string("vprism_8x8x8") ;
    vec_fname[3] = string("vprism_16x16x16") ;
    vec_fname[4] = string("vprism_32x32x32") ;
    vec_fname[5] = string("vprism_64x64x64") ;
    vec_fname[6] = string("vprism_128x128x128") ;
  }
} ;
#if 0
// these cells have hanging nodes
class Mesh_PrysmaticCells : public MeshDirsBase {
public:
  Mesh_PrysmaticCells() : MeshDirsBase(9) { setup() ; }
  ~Mesh_PrysmaticCells() {}
  virtual void setup() { // mesh M3
    fmt_flag  = REGN_FACE_FMT ;
    mesh_dir  = string("./MeshDataSets/Mesh-3D/STANDARD-SET/Prysmatic-Cells-1/") ;
    mesh_name = string("Prysmatic cells with polygonal basis") ;
    vec_fname[0] = string("gdual_1x1x1") ;       // TESTING VERSION
    vec_fname[1] = string("gdual_5x5x5") ;
    vec_fname[2] = string("gdual_10x10x10") ;
    vec_fname[3] = string("gdual_15x15x15") ;
    vec_fname[4] = string("gdual_20x20x20") ;
    vec_fname[5] = string("gdual_25x25x25") ;
    vec_fname[6] = string("gdual_30x30x30") ;
    vec_fname[7] = string("gdual_35x35x35") ;
    vec_fname[8] = string("gdual_40x40x40") ;
  }
} ;
#endif
// --------------------------------------------------------------------------------------------
// previous class   : Mesh_M4
// previous version : Hexa-test-2
// current  version : STANDARD_SET/Random-Hexahedra
// mesh index: 3
class Mesh_RandomHexahedra : public MeshDirsBase {
public:
  Mesh_RandomHexahedra() : MeshDirsBase(6) { setup() ; }
  ~Mesh_RandomHexahedra() {}
  virtual void setup() { // mesh M4
    fmt_flag  = REGN_FACE_FMT ;
    mesh_dir  = string("./MeshDataSets/Mesh-3D/STANDARD-SET/Random-Hexahedra/") ;
    mesh_name = string("Random Hexahedra") ;
    vec_fname[0] = string("gcube.1") ;
    vec_fname[1] = string("gcube.2") ;
    vec_fname[2] = string("gcube.3") ;
    vec_fname[3] = string("gcube.4") ;
    vec_fname[4] = string("gcube.5") ;
    vec_fname[5] = string("gcube.6") ;
  }
} ;
// --------------------------------------------------------------------------------------------
// previous class   : Mesh_M7
// previous version : Voronoi-SuperMesh/Voro-BIG
// current  version : STANDARD_SET/Voronoi-SuperMesh/Voro-BIG
// mesh index: 4
class Mesh_VoronoiBIG : public MeshDirsBase {
public:
  Mesh_VoronoiBIG() : MeshDirsBase(6) { setup() ; }
  ~Mesh_VoronoiBIG() {}
  virtual void setup() { // mesh M7
    fmt_flag  = REGN_FACE_FMT ;
    mesh_dir  = string("./MeshDataSets/Mesh-3D/STANDARD-SET/Voro-BIG/") ;
    mesh_name = string("Voronoi (BIG)") ;
    vec_fname[0] = string("voro-05") ;
    vec_fname[1] = string("voro-10") ;
    vec_fname[2] = string("voro-15") ;
    vec_fname[3] = string("voro-20") ;
    vec_fname[4] = string("voro-25") ;
    vec_fname[5] = string("voro-30") ;
  }
} ;
// --------------------------------------------------------------------------------------------
// previous class   : Mesh_M8
// previous version : Voronoi-SuperMesh/Voro-small-0
// current  version : STANDARD_SET/Voronoi-SuperMesh/Voro-small-0
// mesh index: 5
class Mesh_VoronoiSmall : public MeshDirsBase {
public:
  Mesh_VoronoiSmall() : MeshDirsBase(8) { setup() ; }
  ~Mesh_VoronoiSmall() {}
  virtual void setup() { // mesh M8
    fmt_flag = REGN_FACE_FMT ;
    mesh_dir = string("./MeshDataSets/Mesh-3D/STANDARD-SET/Voro-small-0/") ;
    mesh_name = string("Voronoi (small-0)") ;
    vec_fname[0] = string("voro-2") ;
    vec_fname[1] = string("voro-4") ;
    vec_fname[2] = string("voro-6") ;
    vec_fname[3] = string("voro-8") ;
    vec_fname[4] = string("voro-10") ;
    vec_fname[5] = string("voro-12") ;
    vec_fname[6] = string("voro-14") ;
    vec_fname[7] = string("voro-16") ;
  }
} ;
// this is the one of FVCA6 benchmark
class Mesh_VoronoiSmall1 : public MeshDirsBase {
public:
  Mesh_VoronoiSmall1() : MeshDirsBase(5) { setup() ; }
  ~Mesh_VoronoiSmall1() {}
  virtual void setup() { // mesh 
    fmt_flag = REGN_FACE_FMT ;
    mesh_dir = string("./MeshDataSets/Mesh-3D/STANDARD-SET/Voro-small-1/") ;
    mesh_name = string("Voronoi (small-1)") ;
    vec_fname[0] = string("voro.2") ;
    vec_fname[1] = string("voro.3") ;
    vec_fname[2] = string("voro.4") ;
    vec_fname[3] = string("voro.5") ;
    vec_fname[4] = string("voro.6") ;
  }
} ;
class Mesh_VoronoiSmall2 : public MeshDirsBase {
public:
  Mesh_VoronoiSmall2() : MeshDirsBase(7) { setup() ; }
  ~Mesh_VoronoiSmall2() {}
  virtual void setup() { // mesh 
    fmt_flag = REGN_FACE_FMT ;
    mesh_dir = string("./MeshDataSets/Mesh-3D/STANDARD-SET/Voro-small-2/") ;
    mesh_name = string("Voronoi (small-2)") ;
    vec_fname[0] = string("voro.2") ;
    vec_fname[1] = string("voro.4") ;
    vec_fname[2] = string("voro.6") ;
    vec_fname[3] = string("voro.8") ;
    vec_fname[4] = string("voro.3") ;
    vec_fname[5] = string("voro.5") ;
    vec_fname[6] = string("voro.7") ;
  }
} ;
// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
class Mesh_Voronoi_new_0 : public MeshDirsBase {
public:
  Mesh_Voronoi_new_0() : MeshDirsBase(4) { setup() ; }
  ~Mesh_Voronoi_new_0() {}
  virtual void setup() { // mesh 
    fmt_flag = REGN_FACE_FMT ;
    mesh_dir = string("./MeshDataSets/Mesh-3D/STANDARD-SET/Voro-Cube-0/") ;
    mesh_name = string("Voronoi (new-cube-0)") ;
    vec_fname[0] = string("voro.5") ;
    vec_fname[1] = string("voro.10") ;
    vec_fname[2] = string("voro.20") ;
    vec_fname[3] = string("voro.30") ;
  }
} ;
class Mesh_Voronoi_new_1 : public MeshDirsBase {
public:
  Mesh_Voronoi_new_1() : MeshDirsBase(6) { setup() ; }
  ~Mesh_Voronoi_new_1() {}
  virtual void setup() { // mesh 
    fmt_flag = REGN_FACE_FMT ;
    mesh_dir = string("./MeshDataSets/Mesh-3D/STANDARD-SET/Voro-Tets-1/") ;
    mesh_name = string("Voronoi (new-tets-1)") ;
    vec_fname[0] = string("voro.1") ;
    vec_fname[1] = string("voro.2") ;
    vec_fname[2] = string("voro.3") ;
    vec_fname[3] = string("voro.4") ;
    vec_fname[4] = string("voro.5") ;
    vec_fname[5] = string("voro.6") ;
  }
} ;
class Mesh_Voronoi_new_2 : public MeshDirsBase {
public:
  Mesh_Voronoi_new_2() : MeshDirsBase(6) { setup() ; }
  ~Mesh_Voronoi_new_2() {}
  virtual void setup() { // mesh 
    fmt_flag = REGN_FACE_FMT ;
    mesh_dir = string("./MeshDataSets/Mesh-3D/STANDARD-SET/Voro-Tets-2/") ;
    mesh_name = string("Voronoi (new-tets-2)") ;
    vec_fname[0] = string("voro.1") ;
    vec_fname[1] = string("voro.2") ;
    vec_fname[2] = string("voro.3") ;
    vec_fname[3] = string("voro.4") ;
    vec_fname[4] = string("voro.5") ;
    vec_fname[5] = string("voro.6") ;
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
class Mesh_TiltedCubicCells : public MeshDirsBase {
public:
  Mesh_TiltedCubicCells() : MeshDirsBase(7) { setup() ; }
  ~Mesh_TiltedCubicCells() {}
  virtual void setup() { // mesh M1
    fmt_flag  = REGN_FACE_FMT ;
    mesh_dir  = string("./MeshDataSets/Mesh-3D/TiltedCubicMesh/") ;
    mesh_name = string("Tilted Cubic cells") ;
    vec_fname[0] = string("tgcube_1x1x1") ;      // TESTING VERSION
    vec_fname[1] = string("tgcube_2x2x2") ;
    vec_fname[2] = string("tgcube_4x4x4") ;
    vec_fname[3] = string("tgcube_8x8x8") ;
    vec_fname[4] = string("tgcube_16x16x16") ;
    vec_fname[5] = string("tgcube_32x32x32") ;
    vec_fname[6] = string("tgcube_64x64x64") ;
  }
} ; 


// 
class Mesh_G : public MeshDirsBase {
public:
  Mesh_G() : MeshDirsBase(5) { setup() ; }
  ~Mesh_G() {}
  virtual void setup() { // mesh family G
    fmt_flag  = REGN_FACE_FMT ;
    mesh_dir  = string("./MeshDataSets/Mesh-3D/BENCHMARK-FVCA6/Mesh-G/") ;
    mesh_name = string("Cubic cells with corner refinement") ;
    vec_fname[0] = string("hexa_2x2x2_locraf_new") ;
    vec_fname[1] = string("hexa_4x4x4_locraf_new") ;
    vec_fname[2] = string("hexa_8x8x8_locraf_new") ;
    vec_fname[3] = string("hexa_16x16x16_locraf_new") ;
    vec_fname[4] = string("hexa_32x32x32_locraf_new") ;
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
  case  0 : p_mdir = new Mesh_Cube_Tetrahedra_0     ; break ;
  case  1 : p_mdir = new Mesh_CubicCells            ; break ;
  case  2 : p_mdir = new Mesh_PrysmaticCells        ; break ;
  case  3 : p_mdir = new Mesh_RandomHexahedra       ; break ;
  case  4 : p_mdir = new Mesh_VoronoiBIG            ; break ;
  case  5 : p_mdir = new Mesh_VoronoiSmall          ; break ;
  case  6 : p_mdir = new Mesh_VoronoiSmall1         ; break ;
  //case  7 : p_mdir = new Mesh_VoronoiSmall2         ; break ;
  case  8 : p_mdir = new Mesh_Voronoi_new_0         ; break ;
  case  9 : p_mdir = new Mesh_Voronoi_new_1         ; break ;
  case 10 : p_mdir = new Mesh_Voronoi_new_2         ; break ;
  case 11 : p_mdir = new Mesh_DoubleCubes_RandomHex ; break ;
  case 12 : p_mdir = new Mesh_TiltedCubicCells      ; break ;
  case  7 : p_mdir = new Mesh_G                     ; break ;
  default: assert( bool(false) ) ;
  }
  return *p_mdir ;
}

#endif // end of  _MESH_DIRS
