set(TaskSrc
  ${CMAKE_SOURCE_DIR}/src/Task/Task/Task.cpp
  ${CMAKE_SOURCE_DIR}/src/Task/Task/Task.h
  ${CMAKE_SOURCE_DIR}/src/Task/Task/SettingsParser/MeshSettingsParser.h
  ${CMAKE_SOURCE_DIR}/src/Task/Task/SettingsParser/ParserUtil.h
  ${CMAKE_SOURCE_DIR}/src/Task/Task/SettingsParser/ScheduleSettingsParser.h
  ${CMAKE_SOURCE_DIR}/src/Task/Task/SettingsParser/SettingsParser.h
  ${CMAKE_SOURCE_DIR}/src/Task/Task/SettingsParser/SnapshotSettingsParser.h
  ${CMAKE_SOURCE_DIR}/src/Task/Task/SettingsParser/SolverSettingsParser.h
  ${CMAKE_SOURCE_DIR}/src/Task/Task/SettingsParser/TaskSettingsParser.h
)

set(ElasticSystemSrc 
  ${CMAKE_SOURCE_DIR}/src/Task/ElasticSystem/ElasticSystem.h
  ${CMAKE_SOURCE_DIR}/src/Task/ElasticSystem/IniStates.h
  ${CMAKE_SOURCE_DIR}/src/Task/ElasticSystem/VectorFunctors.h

  ${CMAKE_SOURCE_DIR}/src/Task/ElasticSystem/ElasticSystem.inl
  ${CMAKE_SOURCE_DIR}/src/Task/ElasticSystem/ElasticSystemCommon.inl
  ${CMAKE_SOURCE_DIR}/src/Task/ElasticSystem/SourceFunctors.h

  ${CMAKE_SOURCE_DIR}/src/Task/ElasticSystem/BoundaryInfos.h
  ${CMAKE_SOURCE_DIR}/src/Task/ElasticSystem/ContactInfos.h
  ${CMAKE_SOURCE_DIR}/src/Task/ElasticSystem/ElasticSystemCommon.h
)

set(ElasticVolumeMeshSrc
  ${CMAKE_SOURCE_DIR}/src/Task/ElasticVolumeMesh/Local/ElasticVolumeMesh.h
  ${CMAKE_SOURCE_DIR}/src/Task/ElasticVolumeMesh/Local/ElasticVolumeMesh2.inl
  ${CMAKE_SOURCE_DIR}/src/Task/ElasticVolumeMesh/Local/ElasticVolumeMesh3.inl
  ${CMAKE_SOURCE_DIR}/src/Task/ElasticVolumeMesh/Local/ElasticVolumeMeshCommon.inl

  ${CMAKE_SOURCE_DIR}/src/Task/ElasticVolumeMesh/Distributed/DistributedElasticVolumeMesh.h

  ${CMAKE_SOURCE_DIR}/src/Task/ElasticVolumeMesh/Distributed/DistributedElasticVolumeMesh2.inl
  ${CMAKE_SOURCE_DIR}/src/Task/ElasticVolumeMesh/Distributed/DistributedElasticVolumeMesh3.inl
  ${CMAKE_SOURCE_DIR}/src/Task/ElasticVolumeMesh/Distributed/DistributedElasticVolumeMeshCommon.h
  ${CMAKE_SOURCE_DIR}/src/Task/ElasticVolumeMesh/Distributed/DistributedElasticVolumeMeshCommon.inl
)

set(GeomMeshSrc
  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/AABBTree.h
  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/BucketStorage.h
  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/DisjointSetBuilder.h
  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/TimeHierarchyLevelsManager.h
  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/TopologyReconstructor.h

  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/MeshIO/Distributed/DistributedMeshIO.h
  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/MeshIO/Distributed/DistributedMeshIO.inl
  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/MeshIO/Distributed/DistributedMeshIO2.inl
  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/MeshIO/Distributed/DistributedMeshIO3.inl
  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/MeshIO/Distributed/TransitionInfo.h

  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/MeshIO/Local/MeshIO.h
  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/MeshIO/Local/MeshIO.inl
  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/MeshIO/Local/MeshIO2.inl
  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/MeshIO/Local/MeshIO3.inl
  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/MeshIO/Local/MeshIOBase.h

  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/GeomMesh/GeomMeshCommon.inl
  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/GeomMesh/GeomMeshCommon.h
  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/GeomMesh/GeomMeshBase.inl
  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/GeomMesh/GeomMeshBase.h
  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/GeomMesh/GeomMesh3.inl
  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/GeomMesh/GeomMesh2.inl
  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/GeomMesh/GeomMesh.h
  ${CMAKE_SOURCE_DIR}/src/Task/GeomMesh/GeomMesh/AdditionalCellInfo.h
)

set(VolumeMethodSrc 
  ${CMAKE_SOURCE_DIR}/src/Task/VolumeMethod/Spaces/LagrangeSpace.h
  ${CMAKE_SOURCE_DIR}/src/Task/VolumeMethod/Spaces/PolynomialPrecomputer.h
  ${CMAKE_SOURCE_DIR}/src/Task/VolumeMethod/Spaces/PolynomialSpace.h

  ${CMAKE_SOURCE_DIR}/src/Task/VolumeMethod/FunctionGetters/SourceFunctionGetter.h
  ${CMAKE_SOURCE_DIR}/src/Task/VolumeMethod/FunctionGetters/GhostFunctionGetter.h
  ${CMAKE_SOURCE_DIR}/src/Task/VolumeMethod/FunctionGetters/FunctionGetters.h
  ${CMAKE_SOURCE_DIR}/src/Task/VolumeMethod/FunctionGetters/BoundaryFunctionGetter.h

  ${CMAKE_SOURCE_DIR}/src/Task/VolumeMethod/VolumeMesh/VolumeMesh.h
  ${CMAKE_SOURCE_DIR}/src/Task/VolumeMethod/VolumeMesh/VolumeMesh2.inl
  ${CMAKE_SOURCE_DIR}/src/Task/VolumeMethod/VolumeMesh/VolumeMesh3.inl
  ${CMAKE_SOURCE_DIR}/src/Task/VolumeMethod/VolumeMesh/VolumeMeshCommon.inl
  ${CMAKE_SOURCE_DIR}/src/Task/VolumeMethod/Cell.h
)

set(BINDIR ${CMAKE_SOURCE_DIR}/.bin)

MAKE_DIRECTORY(${BINDIR})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BINDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${BINDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${BINDIR})

IF (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    SET(BinaryFileName taskd)
ENDIF()

IF (${CMAKE_BUILD_TYPE} STREQUAL "Release")
    SET(BinaryFileName task)
ENDIF()

add_executable(${BinaryFileName}
  ${TaskSrc}
  ${NetworkSrc}
  ${MathsSrc}
  ${TinyXmlSrc}
  ${QuadratureIntegrationSrc}
  ${DifferentialSolversSrc}
  ${VtkSrc}
  ${ElasticSystemSrc}
  ${ElasticVolumeMeshSrc}
  ${GeomMeshSrc}
  ${VolumeMethodSrc}
)

target_link_libraries(${BinaryFileName} ${ZLIB_LIBRARIES} ${MPI_LIBRARIES}) 
