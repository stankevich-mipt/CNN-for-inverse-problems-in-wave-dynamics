cmake_minimum_required(VERSION 2.6)
project(Volume)
set(cmake_verbose_makefile off)
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpic++)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mtune=native -funroll-loops -ftree-vectorize -mfpmath=sse -ffast-math -ftree-vectorizer-verbose=0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -mtune=native -funroll-loops -mfpmath=sse -ffast-math -ftree-vectorize -ftree-vectorizer-verbose=0")

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
ENDIF(NOT CMAKE_BUILD_TYPE)

string(TOUPPER ${CMAKE_SYSTEM} BuildPlatform)
message("Generating makefiles for: " ${BuildPlatform} "\n")
message("Build type: " ${CMAKE_BUILD_TYPE} "\n")

SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE ON)

FIND_PACKAGE(MPI)
if (MPI_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_FLAGS}")
  include_directories(${MPI_INCLUDE_PATH})
  message("MPI is found\n")
endif(MPI_FOUND)

#FIND_PACKAGE(OpenMP)
#if (OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
  message("OpenMP is found\n")
#endif (OPENMP_FOUND)

#ZLIB
FIND_PACKAGE(ZLIB REQUIRED)
IF (ZLIB_FOUND)
    MESSAGE(STATUS "Setting up ZLIB")    
    include_directories( ${ZLIB_INCLUDE_DIRS} )    
ENDIF (ZLIB_FOUND)

set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/.bin)

MAKE_DIRECTORY(${CMAKE_SOURCE_DIR}/.bin)

add_definitions(-DTIXML_USE_STL)

set(VtkSrc
  ${CMAKE_SOURCE_DIR}/src/IO/Vtk/BasicVtkWriter.h
  ${CMAKE_SOURCE_DIR}/src/IO/Vtk/ContactVtkWriter.h
  ${CMAKE_SOURCE_DIR}/src/IO/Vtk/ContactVtkWriter.inl
  ${CMAKE_SOURCE_DIR}/src/IO/Vtk/MeshVtkWriter.h
  ${CMAKE_SOURCE_DIR}/src/IO/Vtk/SnapshotVtkWriter.h
  ${CMAKE_SOURCE_DIR}/src/IO/Vtk/VtkReader.h
  ${CMAKE_SOURCE_DIR}/src/IO/Vtk/CellInfoVtkWriter.h
  ${CMAKE_SOURCE_DIR}/src/IO/Vtk/CellInfoVtkWriter.inl
  ${CMAKE_SOURCE_DIR}/src/IO/Segy/SegySeismo.h
  ${CMAKE_SOURCE_DIR}/src/IO/Segy/SegySeismo.inl
)

set(UtilsSrc
  ${CMAKE_SOURCE_DIR}/src/Utils/Base64.h
  ${CMAKE_SOURCE_DIR}/src/Utils/Utils.h
)

set(MathsSrc
  ${CMAKE_SOURCE_DIR}/src/Maths/AABB2.h
  ${CMAKE_SOURCE_DIR}/src/Maths/AABB3.h
  ${CMAKE_SOURCE_DIR}/src/Maths/Coordinates2.h
  ${CMAKE_SOURCE_DIR}/src/Maths/Coordinates3.h
  ${CMAKE_SOURCE_DIR}/src/Maths/Matrix2x2.h
  ${CMAKE_SOURCE_DIR}/src/Maths/Matrix3x3.h
  ${CMAKE_SOURCE_DIR}/src/Maths/MatrixInv.h
  ${CMAKE_SOURCE_DIR}/src/Maths/MatrixMaths.h
  ${CMAKE_SOURCE_DIR}/src/Maths/Random.h
  ${CMAKE_SOURCE_DIR}/src/Maths/Spaces.h
  ${CMAKE_SOURCE_DIR}/src/Maths/Util.h
  ${CMAKE_SOURCE_DIR}/src/Maths/Vector2.h
  ${CMAKE_SOURCE_DIR}/src/Maths/Vector3.h
  ${CMAKE_SOURCE_DIR}/src/Maths/Tensor2.h
  ${CMAKE_SOURCE_DIR}/src/Maths/Tensor3.h
  ${CMAKE_SOURCE_DIR}/src/Maths/VectorMaths.h
  ${CMAKE_SOURCE_DIR}/src/Maths/Collisions.h
)

set(TinyXmlSrc
  ${CMAKE_SOURCE_DIR}/3rdparty/tinyxml/tinystr.cpp
  ${CMAKE_SOURCE_DIR}/3rdparty/tinyxml/tinyxml.cpp
  ${CMAKE_SOURCE_DIR}/3rdparty/tinyxml/tinyxmlerror.cpp
  ${CMAKE_SOURCE_DIR}/3rdparty/tinyxml/tinyxmlparser.cpp
  ${CMAKE_SOURCE_DIR}/3rdparty/tinyxml/tinystr.h
  ${CMAKE_SOURCE_DIR}/3rdparty/tinyxml/tinyxml.h
)

set(DifferentialSolversSrc
  ${CMAKE_SOURCE_DIR}/src/DifferentialSolvers/DifferentialSolver.h
  ${CMAKE_SOURCE_DIR}/src/DifferentialSolvers/DifferentialSystem.h
  ${CMAKE_SOURCE_DIR}/src/DifferentialSolvers/Solvers/AdaptiveRKSolver.h
  ${CMAKE_SOURCE_DIR}/src/DifferentialSolvers/Solvers/DormandPrinceSolver.h
  ${CMAKE_SOURCE_DIR}/src/DifferentialSolvers/Solvers/EulerSolver.h
  ${CMAKE_SOURCE_DIR}/src/DifferentialSolvers/Solvers/RK2Solver.h
  ${CMAKE_SOURCE_DIR}/src/DifferentialSolvers/Solvers/RungeSolver4.h
  ${CMAKE_SOURCE_DIR}/src/DifferentialSolvers/Solvers/TrapezoidSolver.h
  ${CMAKE_SOURCE_DIR}/src/DifferentialSolvers/Solvers/RungeSolver8_10.h
  ${CMAKE_SOURCE_DIR}/src/DifferentialSolvers/Solvers/RungeSolver12_14.h
)

set(NonlinearSolversSrc
  ${CMAKE_SOURCE_DIR}/src/NonlinearSolvers/NonlinearSystemSolver.h
  ${CMAKE_SOURCE_DIR}/src/NonlinearSolvers/NonlinearFunctionSystem.h
  ${CMAKE_SOURCE_DIR}/src/NonlinearSolvers/Solvers/SimpleIterations.h
  ${CMAKE_SOURCE_DIR}/src/NonlinearSolvers/Solvers/KrylovProjector.h
)

set(NetworkSrc
  ${CMAKE_SOURCE_DIR}/src/Network/NetworkInterface.h
  ${CMAKE_SOURCE_DIR}/src/Network/NetworkInterface.cpp
)

set(MetisGKlibSrc
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/b64.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/blas.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/csr.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/error.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/evaluate.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/fkvkselect.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/fs.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/getopt.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/gk_arch.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/gk_defs.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/gk_externs.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/gk_getopt.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/GKlib.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/gk_macros.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/gk_mkblas.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/gk_mkmemory.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/gk_mkpqueue2.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/gk_mkpqueue.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/gk_mkrandom.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/gk_mksort.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/gk_mkutils.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/gk_proto.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/gkregex.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/gkregex.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/gk_struct.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/gk_types.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/gk_graph.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/htable.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/io.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/itemsets.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/mcore.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/memory.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/ms_inttypes.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/ms_stat.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/ms_stdint.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/omp.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/pdb.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/pqueue.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/random.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/rw.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/seq.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/sort.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/string.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/timers.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/tokenizer.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib/gk_util.c
)

set(MetisIncludeSrc
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/include/metis.h
)

set(QuadratureIntegrationSrc
  ${CMAKE_SOURCE_DIR}/3rdparty/quadrature_integration/tetrahedron_arbq_rule.cpp
  ${CMAKE_SOURCE_DIR}/3rdparty/quadrature_integration/triangle_fekete_rule.cpp
)

set(MetisLibSrc
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/auxapi.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/balance.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/bucketsort.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/checkgraph.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/coarsen.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/compress.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/contig.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/debug.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/defs.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/fm.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/fortran.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/frename.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/gklib.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/gklib_defs.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/gklib_rename.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/graph.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/initpart.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/kmetis.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/kwayfm.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/kwayrefine.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/macros.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/mcutil.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/mesh.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/meshpart.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/metislib.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/minconn.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/mincover.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/mmd.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/ometis.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/options.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/parmetis.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/pmetis.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/proto.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/refine.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/rename.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/separator.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/sfm.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/srefine.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/stat.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/stdheaders.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/struct.h
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/timing.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/util.c
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis/wspace.c
)

set(MetisProgramSrc
  ${CMAKE_SOURCE_DIR}/3rdparty/metis/programs/metisbin.h
)

include_directories(${CMAKE_SOURCE_DIR}/3rdparty/metis/GKlib)
include_directories(${CMAKE_SOURCE_DIR}/3rdparty/metis/libmetis)
include_directories(${CMAKE_SOURCE_DIR}/3rdparty/metis/include)
include_directories(${CMAKE_SOURCE_DIR}/3rdparty/metis/programs)
include_directories(${CMAKE_SOURCE_DIR}/3rdparty/eigen)

add_definitions(-DUSE_GKREGEX)
add_definitions(-DNDEBUG)
add_definitions(-DNDEBUG2)
add_definitions(-DEIGEN_DEFAULT_TO_ROW_MAJOR)
add_definitions(-DEIGEN_DONT_PARALLELIZE)

add_definitions(-DEIGEN_NO_DEBUG)
add_definitions(-DEIGEN_NO_AUTOMATIC_RESIZING)
add_definitions(-DEIGEN_FAST_MATH)
add_definitions(-DEIGEN_NO_MALLOC)

add_subdirectory(${CMAKE_SOURCE_DIR}/src/Task/Task)
add_subdirectory(${CMAKE_SOURCE_DIR}/src/Task/MeshBuilder)
add_subdirectory(${CMAKE_SOURCE_DIR}/src/Task/ResultCombiner)