import math as m
import salome
import SMESH
import GEOM
import numpy as np
from salome.geom import geomBuilder
from salome.smesh import smeshBuilder
from numpy import genfromtxt


def resize_list(l, size):
    l += [0] * (max(0, size - len(l)))

#preset for size of a mesh
SIZE_X = 2000.0
SIZE_Y = 2000.0

for j in range(8577, 8580):

    salome.salome_init()
    geompy = geomBuilder.New()
    smesh = smeshBuilder.New()
    #gg = salome.ImportComponentGUI("GEOM")

    # ===============================================================
    # Create a rectangle face with convex heterogenity at it's random point

    pnt1 = geompy.MakeVertex (0.0, 0.0, 0.0)
    pnt2 = geompy.MakeVertex (SIZE_X, 0.0, 0.0)
    pnt3 = geompy.MakeVertex (SIZE_X, SIZE_Y, 0.0)
    pnt4 = geompy.MakeVertex (0.0, SIZE_Y,0)

    geompy.addToStudy(pnt1, "pnt1")
    geompy.addToStudy(pnt2, "pnt2")
    geompy.addToStudy(pnt3, "pnt3")
    geompy.addToStudy(pnt4, "pnt4")

    line1 = geompy.MakeEdge(pnt1, pnt2)
    line2 = geompy.MakeEdge(pnt2, pnt3)
    line3 = geompy.MakeEdge(pnt3, pnt4)
    line4 = geompy.MakeEdge(pnt4, pnt1)
    id_line1 = geompy.addToStudy(line1, "line1")
    id_line2 = geompy.addToStudy(line2, "line2")
    id_line3 = geompy.addToStudy(line3, "line3")
    id_line4 = geompy.addToStudy(line4, "line4")

    base = geompy.MakeFaceWires([line1,line2,line3,line4], 1)
    id_base = geompy.addToStudy(base, "base")
    
  
    #create convex heterogenity 

    path_to_shape = 'shapes/'
    heterogenity = genfromtxt(path_to_shape + f'shape_batch2_{j}.csv', delimiter=',')
    
    pts = list()
    for k in range(heterogenity.shape[0]):
        tmp = geompy.MakeVertex(heterogenity[k][0], heterogenity[k][1], 0.0)
        id_tmp = geompy.addToStudy(tmp, f"sp{k}") # shape point k
        pts.append(tmp)

    edges = list()
    for k in range(len(pts) - 1):
        tmp = geompy.MakeEdge(pts[k], pts[k+1])
        id_tmp = geompy.addToStudy(tmp, f"se{k}")
        edges.append(tmp)

    closure = geompy.MakeEdge(pts[len(pts) - 1], pts[0])
    id_closure = geompy.addToStudy(tmp, f"se{len(pts) - 1}")
    edges.append(closure)

    heterogenity = geompy.MakeFaceWires(edges, 1)
    id_heterogenity = geompy.addToStudy(heterogenity, "heterogeneity")

                  
    partition = geompy.MakePartition([base, heterogenity])
    id_partition = geompy.addToStudy(partition,"Partition")

    # ===============================================================
    # Create a separate geometry for detectors

    pnt5 = geompy.MakeVertex (0.0, 0.0, 0.0)
    pnt6 = geompy.MakeVertex (SIZE_X, 0.0, 0.0)
    detectors = geompy.MakeEdge(pnt5, pnt6)
    geompy.addToStudy(detectors, f"detectors")
    

    # ===============================================================
    # Create a separate geometry group for each of faces and edges

    faceList = geompy.SubShapeAll(partition,geompy.ShapeType["FACE"])

    for face in faceList:
        name = geompy.SubShapeName(face, partition)
        geompy.addToStudyInFather(partition, face, name)

    edgeList = geompy.SubShapeAll(partition ,geompy.ShapeType["EDGE"])

    for edge in edgeList:
        name = geompy.SubShapeName(edge, partition)
        geompy.addToStudyInFather(partition, edge, name) 

    groups_faces = {}
    for i, face in enumerate(faceList):
        groups_faces[i] = geompy.CreateGroup(partition, geompy.ShapeType["FACE"])
        geompy.UnionList(groups_faces[i], [face])

    groups_edges = {}
    contacts = []

    for i, edge in enumerate(edgeList):
        if i < 4:
            groups_edges[i] = geompy.CreateGroup(partition, geompy.ShapeType["EDGE"])
            geompy.UnionList(groups_edges[i], [edge])
        else:
            contacts.append(edge) 
      
    group_contacts = geompy.CreateGroup(partition, geompy.ShapeType["EDGE"])
    geompy.UnionList(group_contacts, contacts)
    # ===============================================================
    # Create a main mesh over the partition, submeshes for each
    # edge, and a submesh for heterogenity with it's wire

         
    mesh = smesh.Mesh(partition)

    for i, group in enumerate(groups_faces.values()):
        SmeshGroup = mesh.GroupOnGeom(group, "faces<{}>".format(i))
    for i, group in enumerate(groups_edges.values()):
        SmeshGroup = mesh.GroupOnGeom(group, "edges{{{}}}".format(i))
    SmeshGroup = mesh.GroupOnGeom(group_contacts, "contacts[0]")

    Regular_1D = mesh.Segment()
    Nb_Segments = Regular_1D.NumberOfSegments(8)
    Nb_Segments.SetDistrType(0)
    mefisto = mesh.Triangle(smeshBuilder.MEFISTO)
    mefisto.LengthFromEdges()

    Regular_1D_1 = mesh.Segment(geom=groups_faces[1])
    Nb_Segments_2 = Regular_1D_1.Adaptive(50, 200, 0.1)
    isDone = mesh.Compute()
    if not isDone : print("Mesh is not computed")

    # ===============================================================
    # Create a detector mesh over the detectors geometry

    detectorMesh = smesh.Mesh(detectors)   
    Regular_1D = detectorMesh.Segment()
    Nb_Segments = Regular_1D.NumberOfSegments(49)
    Nb_Segments.SetDistrType(0)
    isDone = detectorMesh.Compute()
    if not isDone : print("Mesh for detectors is not computed")

    # ===============================================================
    # Export main mesh and detector mesh to a solver's format

    detectorMeshes = [detectorMesh]
    outFolder = 'OutMeshes'

    salomeCellType = None
    salomeBoundaryType = None

    outFileName = outFolder + '/%(meshName)s' % {'meshName': f'heterogeneity_{j}'}
    outFile = open(outFileName + '.mesh', 'w')
    paramsFile = open(outFileName + '.params', 'wb')

    dimension = mesh.MeshDimension()

    if dimension == 2:
        salomeCellType = SMESH.FACE
        salomeBoundaryType = SMESH.EDGE
    elif dimension == 3:
        salomeCellType = SMESH.VOLUME
        salomeBoundaryType = SMESH.FACE

    print("Mesh: ", mesh.GetName())
    print("Dimension: ", dimension)
    nodeIds = mesh.GetNodesId()
    outFile.write('%(nodesCount)d\n' % {'nodesCount': len(nodeIds)})

    for nodeId in nodeIds:
        xyz = mesh.GetNodeXYZ(nodeId)
        for coordIndex in range(dimension):
            outFile.write('%(coord)f ' % {'coord': xyz[coordIndex]})
        outFile.write('\n')

    cells = mesh.GetElementsByType(salomeCellType)
    outFile.write('\n')
    outFile.write('%(indicesCount)d\n' % {'indicesCount': len(cells) * (dimension + 1)})

    for cell in cells:
        cellNodes = mesh.GetElemNodes(cell)
        for nodeIndex in cellNodes:
            outFile.write('%(node)d ' % {'node': nodeIndex - 1})
        outFile.write('\n')

    maxContactTypeNumber = 0
    edgesOfContactTypeCount = []

    maxBoundaryTypeNumber = 0
    edgesOfBoundaryTypeCount = []

    groupTypeNumber = []
    groupTypes = []  # 0 is contact 1 is boundary

    totalContactEdges = 0
    totalBoundaryEdges = 0

    for groupIndex in range(len(mesh.GetGroups())):
        
        group = mesh.GetGroups()[groupIndex]

        groupTypeNumber += [-1]
        groupTypes += [-1]

        if group.GetType() == salomeBoundaryType:
            name = group.GetName()
            print(name)
            number = int(-1)
            try:
                number = int(name.split('[')[1].split(']')[0])
            except:
                pass

            if number >= 0:
                # print "contact %d" % groupIndex
                if number > maxContactTypeNumber:
                    maxContactTypeNumber = number
                groupTypeNumber[groupIndex] = number
                groupTypes[groupIndex] = 0  # contact
                resize_list(edgesOfContactTypeCount, number + 1)
                edgesOfContactTypeCount[number] += len(group.GetIDs())
                totalContactEdges += len(group.GetIDs())

            number = int(-1)
            try:
                number = int(name.split('{')[1].split('}')[0])
            except:
                pass

            if number >= 0:
                # print "boundary %d" % groupIndex
                if number > maxBoundaryTypeNumber:
                    maxBoundaryTypeNumber = number
                groupTypeNumber[groupIndex] = number
                groupTypes[groupIndex] = 1  # boundary
                resize_list(edgesOfBoundaryTypeCount, number + 1)
                edgesOfBoundaryTypeCount[number] += len(group.GetIDs())
                totalBoundaryEdges += len(group.GetIDs())

    outFile.write('\n%(totalContactEdges)d\n' % {'totalContactEdges': totalContactEdges})
    for typeIndex in range(len(edgesOfContactTypeCount)):
        for groupIndex in range(len(mesh.GetGroups())):
            group = mesh.GetGroups()[groupIndex]
            if groupTypeNumber[groupIndex] == typeIndex and groupTypes[groupIndex] == 0:
                for edgeId in group.GetIDs():
                    edgeNodes = mesh.GetElemNodes(edgeId)
                    # double because left side == right side and will be duplicated in mesh builder
                    for num in range(2):
                        for nodeIndex in range(dimension):
                            outFile.write('%(node)d ' % {'node': edgeNodes[nodeIndex] - 1})

    print("Contact types: %d " % len(edgesOfContactTypeCount))
    outFile.write('\n')
    outFile.write('%d ' % len(edgesOfContactTypeCount))
    outFile.write('\n')
    for typeIndex in range(len(edgesOfContactTypeCount)):
        outFile.write('%d ' % edgesOfContactTypeCount[typeIndex])

    outFile.write('\n%(totalBoundaryEdges)d\n' % {'totalBoundaryEdges': totalBoundaryEdges})

    for typeIndex in range(len(edgesOfBoundaryTypeCount)):
        for groupIndex in range(len(mesh.GetGroups())):
            group = mesh.GetGroups()[groupIndex]
            if groupTypeNumber[groupIndex] == typeIndex and groupTypes[groupIndex] == 1:
                for edgeId in group.GetIDs():
                    edgeNodes = mesh.GetElemNodes(edgeId)
                    for nodeIndex in range(dimension):
                        outFile.write('%(node)d ' % {'node': edgeNodes[nodeIndex] - 1})

    print("Boundary types: %d " % len(edgesOfBoundaryTypeCount))
    outFile.write('\n')
    outFile.write('%d ' % len(edgesOfBoundaryTypeCount))
    outFile.write('\n')
    for typeIndex in range(len(edgesOfBoundaryTypeCount)):
        outFile.write('%d ' % edgesOfBoundaryTypeCount[typeIndex])
        print("Boundary edges of type %(type)d: %(count)d " % {'type' : typeIndex, 'count' : edgesOfBoundaryTypeCount[typeIndex]})


    totalDetectorsCount = 0
    for detectorMesh in detectorMeshes:
        totalDetectorsCount += len(detectorMesh.GetNodesId())

    print("Detectors count: %d\n" % totalDetectorsCount)
    outFile.write('%d ' % totalDetectorsCount)

    for detectorMesh in detectorMeshes:
        for nodeId in detectorMesh.GetNodesId():
            xyz = detectorMesh.GetNodeXYZ(nodeId)
            for coordIndex in range(dimension):
                outFile.write('%(coord)f ' % {'coord': xyz[coordIndex]})
            outFile.write('\n')

    elements = mesh.GetElementsId()
    elementSubmeshIndices = len(elements) * [0]  # we'll have edge elements as dummies in this array as well but they won't be used

    print("Total cells count: %d" % len(cells))

    for groupIndex in range(len(mesh.GetGroups())):
        group = mesh.GetGroups()[groupIndex]
        groupType = -1

        if group.GetType() == salomeCellType:
            name = group.GetName()
            try:
                groupType = int(name.split('<')[1].split('>')[0])
            except:
                pass
            print("Cell group" + name + "is index %d" % groupType)

        if groupType >= 0:
            for elementIndex in group.GetIDs():
                if mesh.GetElementType(elementIndex, True) == salomeCellType:
                    elementSubmeshIndices[elementIndex - 1] = groupType

    cellSubmeshIndices = []
    for elementIndex in range(len(elementSubmeshIndices)):
        if mesh.GetElementType(elementIndex + 1, True) == salomeCellType:
            cellSubmeshIndices.append(elementSubmeshIndices[elementIndex])

    print("Cell params count %d" % len(cellSubmeshIndices))

    paramsFileByteArray = bytearray(cellSubmeshIndices)

    print("Cell params array size %d" % len(paramsFileByteArray))
    paramsFile.write(paramsFileByteArray)

    outFile.close()
    paramsFile.close()

    print("Exporting done")

    #clean up everything
    salome.salome_close()

