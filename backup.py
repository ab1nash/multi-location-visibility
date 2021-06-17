def dot(vA, vB):
    return vA[0]*vB[0]+vA[1]*vB[1]
def ang(lineA, lineB):
    # Get nicer vector form
    vA = [(lineA[0][0]-lineA[1][0]), (lineA[0][1]-lineA[1][1])]
    vB = [(lineB[0][0]-lineB[1][0]), (lineB[0][1]-lineB[1][1])]
    # Get dot prod
    dot_prod = dot(vA, vB)
    # Get magnitudes
    magA = dot(vA, vA)**0.5
    magB = dot(vB, vB)**0.5
    # Get cosine value
    cos_ = dot_prod/magA/magB
    # Get angle in radians and then convert to degrees
    angle = math.acos(dot_prod/magB/magA)
    # Basically doing angle <- angle mod 360
    ang_deg = math.degrees(angle)%360

    if ang_deg-180>=0:
        # As in if statement
        return 360 - ang_deg
    else: 
        return ang_deg

# ------------------------------------------------------------

def getLength(vrset:list, queryPoint:Point):
    '''
        Plain visible region calculation
        @args:
            visible region:[segments]
            query point:Point
        @returns:
            length:int
    '''
    lineA = [[queryPoint.x,queryPoint.y],[]]
    lineB = [[],[]]
    dist=0
    for vr in vrset:
        if(vr):
            # print(vr)
            src = asPoint(vr.source())
            tgt = asPoint(vr.target())
            start = [src.x,src.y]
            end = [tgt.x,tgt.y]
            lineA[1] = [(start[0]+end[0])/2,(start[1]+end[1])/2]
            lineB[0] = start
            lineB[1] = end
            angle = ang(lineA,lineB)
            dist+=(angle/90)*src.distance(tgt)
    
    # print(dist)
    return dist


# ------------------------------------------------------------

def buildVRPolygons(vis:list, queryPoint:Point):
    '''
        Return set of VR Polygons
        @args:
            visible region:[segments]
            query point:Point
        @returns:
            vrSet:int
    '''
    vrSet = []
    for seg in vis:
        VRvertices = [queryPoint]
        if(seg.source() not in VRvertices):
            VRvertices.append((seg.source().x(),seg.source().y()))
        if(seg.target() not in VRvertices):
            VRvertices.append((seg.target().x(),seg.target().y()))
        vrSet.append(Polygon(VRvertices))
    
    return vrSet

# ------------------------------------------------------------

def insideVisibleRegion(visibleRegionSet:list, node:Polygon):
    '''
        Test if VR and (skgeom)polygon intersect
            @args:
                visibilityRegionSet:[Polygons] --> VR from query point
                node:[segments] --> obstacle
            @returns:
                bool --> True when polygon and VR intersect, False o.w.
    '''
    # node scaling
    
    for vrPolygon in visibleRegionSet:
        intersection = node.boundary.intersection(vrPolygon.boundary)
        print(intersection)
        if isinstance(intersection,MultiPoint):
            return True
    
    return False

# ------------------------------------------------------------

def kMaximumVisibility(T,Q,k, gdf):
    '''
        Main algorithm from MV paper
    '''
    L = np.empty(len(Q.index))
    VQ = PriorityQueue()
    LO = np.full(len(Q.index),PriorityQueue())
    CO = set()
    vrset = np.empty(len(Q.index),dtype=object)
    vrPoly = np.empty(len(Q.index),dtype=object)
    end = False
    cont = True
    node = 0 # ?
    trav = rtree1.leaves()
    for ix in Q.index:
        # if(ix==4):
        queryPoint = Point(Q['geometry'][ix].x*10000,Q['geometry'][ix].y*10000)
        if(Boundary.contains(queryPoint)== True):
            vrset[ix] = updateVisibility(asPoint2(queryPoint),T,arr)
            vrPoly[ix] = buildVRPolygons(vrset[ix],queryPoint)
            val = getLength(vrset[ix], queryPoint)
            VQ.put((val, Q['id'][ix]))
    print("Phase 1 OK!")
    # 1.11
    for nodeList in rtree1.leaves():
        for nodeId in nodeList[1]:
            if(end == True):
                break
            # 1.12
            # print(nodeId)
            if(nodeId not in CO):
                node = gdf.iloc[nodeId].geometry
                for ix in Q.index:
                    if(ix==4):
                        queryPoint = Point(Q['geometry'][ix].x*10000,Q['geometry'][ix].y*10000)
                        if(Boundary.contains(queryPoint)== True):
                            if(insideVisibleRegion(vrPoly[ix],node)):
                                LO[ix].put(nodeId,minDist(queryPoint,node))
                CO.add(nodeId)
            # -------------------------------------------------------
            cont = True



ans = kMaximumVisibility(targetBuilder, qpGdf, 2, gdf)