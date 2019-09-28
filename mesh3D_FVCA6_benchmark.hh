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
// --------------------------------------------------------------------------------------------
// this mesh, whichs correspond to flag 1, is not in the FVCA6 benchmark
// mesh index: 0
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
// --------------------------------------------------------------------------------------------
// mesh index: 1
#if 1 // regular hexahedra (cubes)
class Mesh_A : public MeshDirsBase {
public:
  Mesh_A() : MeshDirsBase(5) { setup() ; }
  ~Mesh_A() {}
  virtual void setup() { // mesh family A
    fmt_flag  = REGN_FACE_FMT ;
    mesh_dir  = string("./MeshDataSets/Mesh-3D/BENCHMARK-FVCA6/Mesh-A/") ;
    mesh_name = string("Regular hexahedra (2x2x2)") ;
    vec_fname[0] = string("hexa_2x2x2") ;
    vec_fname[1] = string("hexa_4x4x4") ;
    vec_fname[2] = string("hexa_8x8x8") ;
    vec_fname[3] = string("hexa_16x16x16") ;
    vec_fname[4] = string("hexa_32x32x32") ;
  }
} ;
#endif
class Mesh_RandomHexahedra : public MeshDirsBase {
public:
  Mesh_RandomHexahedra() : MeshDirsBase(5) { setup() ; }
  ~Mesh_RandomHexahedra() {}
  virtual void setup() {
    fmt_flag  = REGN_FACE_FMT ;
    mesh_dir  = string("./MeshDataSets/Mesh-3D/STANDARD-SET/Random-Hexahedra/") ;
    mesh_name = string("Random Hexahedra") ;
    vec_fname[0] = string("gcube.1") ;
    vec_fname[1] = string("gcube.2") ;
    vec_fname[2] = string("gcube.3") ;
    vec_fname[3] = string("gcube.4") ;
    //    vec_fname[4] = string("gcube.5") ; does not exist!!!
  }
} ;
// --------------------------------------------------------------------------------------------
// mesh index: 2
class Mesh_B : public MeshDirsBase {
public:
  Mesh_B() : MeshDirsBase(8) { setup() ; }
  ~Mesh_B() {}
  virtual void setup() { // mesh family B
    fmt_flag  = REGN_FACE_FMT ;
    mesh_dir  = string("./MeshDataSets/Mesh-3D/BENCHMARK-FVCA6/Mesh-B/") ;
    mesh_name = string("tetrahedra") ;
    vec_fname[0] = string("tet.0") ;
    vec_fname[1] = string("tet.1") ;
    vec_fname[2] = string("tet.2") ;
    vec_fname[3] = string("tet.3") ;
    vec_fname[4] = string("tet.4") ;
    vec_fname[5] = string("tet.5") ;
    vec_fname[6] = string("tet.6") ;
    vec_fname[7] = string("tet.00") ;
  }
} ; 
// --------------------------------------------------------------------------------------------
// mesh index: 3
class Mesh_C : public MeshDirsBase {
public:
  Mesh_C() : MeshDirsBase(5) { setup() ; }
  ~Mesh_C() {}
  virtual void setup() { // mesh family C
    fmt_flag  = REGN_FACE_FMT ;
    mesh_dir  = string("./MeshDataSets/Mesh-3D/BENCHMARK-FVCA6/Mesh-C/") ;
    mesh_name = string("Voronoi cells") ;
    vec_fname[0] = string("vmesh_1") ;
    vec_fname[1] = string("vmesh_2") ;
    vec_fname[2] = string("vmesh_3") ;
    vec_fname[3] = string("vmesh_4") ;
    vec_fname[4] = string("vmesh_5") ;
  }
} ;
// --------------------------------------------------------------------------------------------
// mesh index: 4
class Mesh_D : public MeshDirsBase {
public:
  Mesh_D() : MeshDirsBase(4) { setup() ; }
  ~Mesh_D() {}
  virtual void setup() { // mesh family D
    fmt_flag  = REGN_FACE_FMT ;
    mesh_dir  = string("./MeshDataSets/Mesh-3D/BENCHMARK-FVCA6/Mesh-D/") ;
    mesh_name = string("Double Kershaw") ;
    vec_fname[0] = string("dkershaw08") ;
    vec_fname[1] = string("dkershaw16") ;
    vec_fname[2] = string("dkershaw32") ;
    vec_fname[3] = string("dkershaw64") ;
  }
} ;
// --------------------------------------------------------------------------------------------
// mesh index: 5
class Mesh_E : public MeshDirsBase {
public:
  Mesh_E() : MeshDirsBase(4) { setup() ; }
  ~Mesh_E() {}
  virtual void setup() { // mesh family E
    fmt_flag  = REGN_FACE_FMT ;
    mesh_dir  = string("./MeshDataSets/Mesh-3D/BENCHMARK-FVCA6/Mesh-E/") ;
    mesh_name = string("Prysm with triangular base") ;
    vec_fname[0] = string("bls_10") ;
    vec_fname[1] = string("bls_20") ;
    vec_fname[2] = string("bls_30") ;
    vec_fname[3] = string("bls_40") ;
  }
} ;
// --------------------------------------------------------------------------------------------
//this is the BENCHMARK version of Layered/Quad-1, it gives the same results
// mesh index: 6
class Mesh_F : public MeshDirsBase {
public:
  Mesh_F() : MeshDirsBase(4) { setup() ; }
  ~Mesh_F() {}
  virtual void setup() { // mesh family F
    fmt_flag  = REGN_FACE_FMT ;
    mesh_dir  = string("./MeshDataSets/Mesh-3D/BENCHMARK-FVCA6/Mesh-F/") ;
    mesh_name = string("Prysm with polygonal base") ;
    vec_fname[0] = string("dbls_10") ;
    vec_fname[1] = string("dbls_20") ;
    vec_fname[2] = string("dbls_30") ;
    vec_fname[3] = string("dbls_40") ; 
  }
} ;
// --------------------------------------------------------------------------------------------
// mesh index: 7
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
// mesh index: 8
class Mesh_H : public MeshDirsBase {
public:
  Mesh_H() : MeshDirsBase(5) { setup() ; }
  ~Mesh_H() {}
  virtual void setup() { // mesh family H
    fmt_flag  = REGN_FACE_FMT ;
    mesh_dir  = string("./MeshDataSets/Mesh-3D/BENCHMARK-FVCA6/Mesh-H/") ;
    mesh_name = string("Cells with local refinement") ;
    // ----------------------------------------
    vec_fname[0]  = string("locrafgrid_1") ;
    vec_fname[1]  = string("locrafgrid_2") ;
    vec_fname[2]  = string("locrafgrid_3") ;
    vec_fname[3]  = string("locrafgrid_4") ;
    vec_fname[4]  = string("locrafgrid_5") ;
    // ----------------------------------------
    //     vec_fname[5]  = string("locrafgrid_1_new") ;
    //     vec_fname[6]  = string("locrafgrid_2_new") ;
    //     vec_fname[7]  = string("locrafgrid_3_new") ;
    //     vec_fname[8]  = string("locrafgrid_4_new") ;
    //     vec_fname[9]  = string("locrafgrid_5_new") ;
    //     vec_fname[10] = string("locrafgrid_6_new") ;
    //     // ----------------------------------------
  }
} ;
// mesh index: 81
class Mesh_H1 : public MeshDirsBase {
public:
  Mesh_H1() : MeshDirsBase(6) { setup() ; }
  ~Mesh_H1() {}
  virtual void setup() { // mesh family H
    fmt_flag  = REGN_FACE_FMT ;
    mesh_dir  = string("./MeshDataSets/Mesh-3D/BENCHMARK-FVCA6/Mesh-H/") ;
    mesh_name = string("Cells with local refinement (new)") ;
      // ----------------------------------------
      //     vec_fname[0]  = string("locrafgrid_1") ;
      //     vec_fname[1]  = string("locrafgrid_2") ;
      //     vec_fname[2]  = string("locrafgrid_3") ;
      //     vec_fname[3]  = string("locrafgrid_4") ;
      //     vec_fname[4]  = string("locrafgrid_5") ;
      //     // ----------------------------------------
    int n=5 ;
    vec_fname[5-n]  = string("locrafgrid_1") ;     // (the 1st mesh is the same)
    vec_fname[6-n]  = string("locrafgrid_2_new") ;
    vec_fname[7-n]  = string("locrafgrid_3_new") ;
    vec_fname[8-n]  = string("locrafgrid_4_new") ;
    vec_fname[9-n]  = string("locrafgrid_5_new") ;
    vec_fname[10-n] = string("locrafgrid_6_new") ;
    // ----------------------------------------
  }
} ;
// --------------------------------------------------------------------------------------------
// mesh index: 9
class Mesh_I : public MeshDirsBase {
public:
  Mesh_I() : MeshDirsBase(4) { setup() ; }
  ~Mesh_I() {}
  virtual void setup() { // mesh family I
    fmt_flag  = REGN_FACE_FMT ;
    mesh_dir  = string("./MeshDataSets/Mesh-3D/BENCHMARK-FVCA6/Mesh-I/") ;
    mesh_name = string("Checkerboard mesh") ;
    vec_fname[0] = string("damier_1") ;
    vec_fname[1] = string("damier_2") ;
    vec_fname[2] = string("damier_3") ;
    vec_fname[3] = string("damier_4") ;
  }
} ;
// --------------------------------------------------------------------------------------------
// mesh index: 10
class Mesh_AA : public MeshDirsBase {
public:
  Mesh_AA() : MeshDirsBase(4) { setup() ; }
  ~Mesh_AA() {}
  virtual void setup() { // mesh family AA
    fmt_flag  = REGN_FACE_FMT ;
    mesh_dir  = string("./MeshDataSets/Mesh-3D/BENCHMARK-FVCA6/Mesh-AA/") ;
    mesh_name = string("Random mesh") ;
    vec_fname[0] = string("RandMesh4") ;
    vec_fname[1] = string("RandMesh8") ;
    vec_fname[2] = string("RandMesh16") ;
    vec_fname[3] = string("RandMesh32") ;
  }
} ;
// --------------------------------------------------------------------------------------------
// mesh index: 11
class Mesh_BB : public MeshDirsBase {
public:
  Mesh_BB() : MeshDirsBase(7) { setup() ; }
  ~Mesh_BB() {}
  virtual void setup() { // mesh family BB
    fmt_flag  = REGN_FACE_FMT ;
    mesh_dir  = string("./MeshDataSets/Mesh-3D/BENCHMARK-FVCA6/Mesh-BB/") ;
    mesh_name = string("Well cells") ;
    vec_fname[0] = string("WellMesh05") ;
    vec_fname[1] = string("WellMesh08") ;
    vec_fname[2] = string("WellMesh12") ;
    vec_fname[3] = string("WellMesh16") ;
    vec_fname[4] = string("WellMesh21") ;
    vec_fname[5] = string("WellMesh26") ;
    vec_fname[6] = string("WellMesh32") ;
  }
} ;
// --------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------
//---- set mesh dirs
// --------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------
//MeshDirsBase & set_mesh_dirs( RPars & rpar ) {
//int mesh_flag = rpar.get_mesh_flag() ;
// --------------------------------------------
MeshDirsBase & set_mesh_dirs( int mesh_flag ) {
  MeshDirsBase * p_mdir ;
  switch ( mesh_flag ) {
  case  0 : p_mdir = new Mesh_CubicCells ; break ; //(not in FVCA6 benchmark set)
    //case  1 : p_mdir = new Mesh_A  ; break ; // regular hexa
  case  1 : p_mdir = new Mesh_RandomHexahedra ; break ; 
  case  2 : p_mdir = new Mesh_B  ; break ; // tetrahedra
  case  3 : p_mdir = new Mesh_C  ; break ; // voronoi
  case  4 : p_mdir = new Mesh_D  ; break ; // kershaw
  case  5 : p_mdir = new Mesh_E  ; break ; // prysm, triangular base
  case  6 : p_mdir = new Mesh_F  ; break ; // prysm, polygonal base
  case  7 : p_mdir = new Mesh_G  ; break ; // hexa + local refinement
  case  8 : p_mdir = new Mesh_H  ; break ; // hexa + strip refinement, 1st seq
  case 81 : p_mdir = new Mesh_H1 ; break ; // hexa + strip refinement, 2nd seq
  case  9 : p_mdir = new Mesh_I  ; break ; // checkerboard (damier)
  case 10 : p_mdir = new Mesh_AA ; break ; // random mesh
  case 11 : p_mdir = new Mesh_BB ; break ; // well mesh
  default: assert( bool(false) ) ;
  }
  string str_dir = p_mdir->get_mesh_dir() ;
  PRT(str_dir) ;
  return *p_mdir ;
}
#endif // end of  _MESH_DIRS
