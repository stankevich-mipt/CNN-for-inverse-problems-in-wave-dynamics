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
SIZE_X = 100.
SIZE_Y = 50.

for j in range(10000):

    salome.salome_init()
    geompy = geomBuilder.New()
    smesh = smeshBuilder.New()
    
    # base box

    pnt1 = geompy.MakeVertex(0.0, 0.0,  0.0)
    pnt2 = geompy.MakeVertex(SIZE_X, 0.0, 0.0)
    pnt3 = geompy.MakeVertex(SIZE_X, SIZE_Y, 0.0)
    pnt4 = geompy.MakeVertex(0.0, SIZE_Y, 0)
    pnt5 = geompy.MakeVertex(0.0, 0.0,  50.0)
    pnt6 = geompy.MakeVertex(0.0, 0.0, -50.0)

    v1   = geompy.MakeVector(pnt5, pnt6)

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
    
    #add boundary between two surfaces 

    path_to_shape = 'shapes/'
    heterogenity = genfromtxt(path_to_shape + f'bone_{j}.csv', delimiter=',')
    
    pts = list()
    for k in range(heterogenity.shape[0]):
        tmp = geompy.MakeVertex(heterogenity[k][0], heterogenity[k][1], 0.0)
        id_tmp = geompy.addToStudy(tmp, f"cp{k}") # shape point k
        pts.append(tmp)

    edges = list()
    for k in range(len(pts) - 1):
        tmp = geompy.MakeEdge(pts[k], pts[k+1])
        id_tmp = geompy.addToStudy(tmp, f"ce{k}")
        edges.append(tmp)

    curve    = geompy.MakeWire(edges)
    id_curve = geompy.addToStudy(curve, "curve") 

    extr    = geompy.MakePrismVecH2Ways(curve, v1, 10)
    id_extr = geompy.addToStudy(extr, "extr")   

    # base surface 

    partition    = geompy.MakePartition([base], [extr])
    id_partition = geompy.addToStudy(partition, "partition")

    # Create a separate geometry group for each of faces and edges

    faceList = geompy.SubShapeAll(partition,geompy.ShapeType["FACE"])

    for face in faceList:
        name = geompy.SubShapeName(face, partition)
        geompy.addToStudyInFather(partition, face, name)


    wires = geompy.SubShapeAll(partition, geompy.ShapeType["WIRE"])
    edges = geompy.SubShapeAll(partition, geompy.ShapeType["EDGE"])

    for edge in edges:
        name = geompy.SubShapeName(edge, partition)
        geompy.addToStudyInFather(partition, edge, name) 

    groups_wires = [
        geompy.SubShapeAll(wire, geompy.ShapeType["EDGE"]) for wire in wires]
   
    upperWireEdges = geompy.CreateGroup(partition, geompy.ShapeType["EDGE"])
    geompy.UnionList(upperWireEdges, edges)
    geompy.DifferenceList(upperWireEdges, groups_wires[0])

    lowerWireEdges = geompy.CreateGroup(partition, geompy.ShapeType["EDGE"])
    geompy.UnionList(lowerWireEdges, edges)
    geompy.DifferenceList(lowerWireEdges, groups_wires[1])

    foo = geompy.CreateGroup(partition, geompy.ShapeType["EDGE"])
    geompy.UnionList(foo, groups_wires[0])
    bar = geompy.CreateGroup(partition, geompy.ShapeType["EDGE"])
    geompy.UnionList(bar, groups_wires[1])
    
    contactEdges = geompy.IntersectGroups(foo, bar)
    geompy.addToStudyInFather(partition, contactEdges, 'contactEdges')

    boundaryEdges = geompy.UnionGroups(upperWireEdges, lowerWireEdges)
    geompy.addToStudyInFather(partition, boundaryEdges, 'boundaryEdges')

    lb, rb, tb, bb = [], [], [], []

    for edge in edges:
        pts = geompy.SubShapeAll(edge, geompy.ShapeType["VERTEX"])
        coord0 = geompy.PointCoordinates(pts[0])
        coord1 = geompy.PointCoordinates(pts[1])
        cmx = (coord0[0] + coord1[0]) / 2
        cmy = (coord0[1] + coord1[1]) / 2
        
        if cmx == SIZE_X : rb.append(edge)
        if cmx == 0.0    : lb.append(edge)
        if cmy == SIZE_Y : tb.append(edge)
        if cmy == 0.0    : bb.append(edge)


    leftBoundaryEdges = geompy.CreateGroup(partition, geompy.ShapeType["EDGE"])
    geompy.UnionList(leftBoundaryEdges, lb)
    geompy.addToStudyInFather(partition, leftBoundaryEdges, 'leftBoundaryEdges')

    rightBoundaryEdges = geompy.CreateGroup(partition, geompy.ShapeType["EDGE"])
    geompy.UnionList(rightBoundaryEdges, rb)
    geompy.addToStudyInFather(partition, rightBoundaryEdges, 'rightBoundaryEdges')

    topBoundaryEdges = geompy.CreateGroup(partition, geompy.ShapeType["EDGE"])
    geompy.UnionList(topBoundaryEdges, tb)
    geompy.addToStudyInFather(partition, topBoundaryEdges, 'topBoundaryEdges')

    botBoundaryEdges = geompy.CreateGroup(partition, geompy.ShapeType["EDGE"])
    geompy.UnionList(botBoundaryEdges, bb)
    geompy.addToStudyInFather(partition, botBoundaryEdges, 'botBoundaryEdges')

    groups_faces = {}
    for i, face in enumerate(faceList):
        groups_faces[i] = geompy.CreateGroup(partition, geompy.ShapeType["FACE"])
        geompy.UnionList(groups_faces[i], [face])

    groups_edges = {
        0: topBoundaryEdges,
        1: botBoundaryEdges, 
        2: leftBoundaryEdges, 
        3: rightBoundaryEdges
    }

    # generate a mesh

    mesh = smesh.Mesh(partition)

    print(smeshBuilder.__dict__)

    for i, group in enumerate(groups_faces.values()):
        SmeshGroup = mesh.GroupOnGeom(group, "faces<{}>".format(i))
    for i, group in enumerate(groups_edges.values()):
        SmeshGroup = mesh.GroupOnGeom(group, "edges{{{}}}".format(i))

    SmeshGroup = mesh.GroupOnGeom(contactEdges, "contacts[0]")
    
    Regular_1D = mesh.Segment()
    Nb_Segments = Regular_1D.Adaptive(5.0, 10.0, 0.1)
    algo = mesh.Triangle(smeshBuilder.NETGEN_2D)
    algo.LengthFromEdges()

    Regular_1D_bb = mesh.Segment(geom=leftBoundaryEdges)
    Nb_Segments = Regular_1D_bb.Adaptive(5.0, 10.0, 0.1)

    Regular_1D_tb = mesh.Segment(geom=rightBoundaryEdges)
    Nb_Segments = Regular_1D_tb.Adaptive(5.0, 10.0, 0.1)

    Regular_1D_contacts = mesh.Segment(geom=contactEdges)
    Nb_Segments = Regular_1D_contacts.NumberOfSegments(1)
    Nb_Segments.SetDistrType(0)


    isDone = mesh.Compute()
    if not isDone : print("Mesh is not computed")
    
    # Export mesh to solver's format 

    detectorMeshes = []
    outFolder = 'OutMeshes'

    salomeCellType = None
    salomeBoundaryType = None

    outFileName = outFolder + '/%(meshName)s' % {'meshName': f'bone_{j}'}
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
    

    


