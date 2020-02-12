set(MeshBuilderSrc
  ${CMAKE_SOURCE_DIR}/src/Task/MeshBuilder/MeshBuilder.cpp
  ${CMAKE_SOURCE_DIR}/src/Task/MeshBuilder/MeshDistributor.h
  ${CMAKE_SOURCE_DIR}/src/Task/MeshBuilder/MeshSplitter/MeshSplitter.h
  ${CMAKE_SOURCE_DIR}/src/Task/MeshBuilder/MeshSplitter/MeshSplitterCommon.h
  ${CMAKE_SOURCE_DIR}/src/Task/MeshBuilder/MeshSplitter/MeshSplitter2.inl
  ${CMAKE_SOURCE_DIR}/src/Task/MeshBuilder/MeshSplitter/MeshSplitter3.inl
  ${CMAKE_SOURCE_DIR}/src/Task/MeshBuilder/MeshSplitter/MeshSplitterCommon.inl
)

set(SpecificMeshBuildersSrc
  ${CMAKE_SOURCE_DIR}/src/Task/MeshBuilder/SpecificMeshBuilders/BasicMeshBuilder.h
  ${CMAKE_SOURCE_DIR}/src/Task/MeshBuilder/SpecificMeshBuilders/SalomeMeshBuilder.h
  ${CMAKE_SOURCE_DIR}/src/Task/MeshBuilder/SpecificMeshBuilders/SpecificMeshBuilders.h
  ${CMAKE_SOURCE_DIR}/src/Task/MeshBuilder/SpecificMeshBuilders/UnstructuredMeshBuilder.h
)

set(ElasticSystemSrc 
  ${CMAKE_SOURCE_DIR}/src/Task/ElasticSystem/ElasticSystem.h
  ${CMAKE_SOURCE_DIR}/src/Task/ElasticSystem/IniStates.h
  ${CMAKE_SOURCE_DIR}/src/Task/ElasticSystem/VectorFunctors.h
)

set(TriGeomSrc
  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/MeshIO/Distributed/DistributedMeshIO.h
  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/MeshIO/Local/MeshIO.h
  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/GeomMesh/GeomMesh.h
)

set(BINDIR ${CMAKE_SOURCE_DIR}/.bin)

MAKE_DIRECTORY(${BINDIR})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BINDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${BINDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${BINDIR})

IF (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    SET(BinaryFileName meshbuilderd)
ENDIF()

IF (${CMAKE_BUILD_TYPE} STREQUAL "Release")
    SET(BinaryFileName meshbuilder)
ENDIF()

add_executable(${BinaryFileName}
  ${MeshBuilderSrc}
  ${SpecificMeshBuildersSrc}
  ${VtkSrc}
  ${MathsSrc}
  ${TinyXmlSrc}
  ${ElasticSystemSrc}
  ${TriGeomSrc}
  ${TriangleSrc}
  ${MetisGKlibSrc}
  ${MetisIncludeSrc}
  ${MetisLibSrc}
  ${MetisProgramSrc}
)

target_link_libraries(${BinaryFileName} ${ZLIB_LIBRARIES}) 
			
