# # To add a new cell, type '# %%'
# # To add a new markdown cell, type '# %% [markdown]'
# # %%
# from IPython import get_ipython

# # %% [markdown]
# # # QD Approach to solving k-maximum visibility

# # %%
# # [TODO]: Add 10000 scaling properly


# %%
import math
import fiona
import numpy as np
import osmnx as ox
import pandas as pd
import networkx as nx
import geopandas as gp


import matplotlib.pyplot as plt
from random import random
from queue import PriorityQueue
from shapely.ops import nearest_points
from shapely.geometry import Polygon,Point, MultiPoint, box
from skgeom import Point2,Segment2,arrangement,intersection as skInter,RotationalSweepVisibility, draw

# get_ipython().run_line_magic('matplotlib', 'inline')

np.random.seed(2021)

# %% [markdown]
# ## skgeom functions

# %%
def asPoint(point2):
    return Point(point2.bbox().xmax(),point2.bbox().ymax())


# %%
def asPoint2(point):
    return Point2(point.x,point.y)


# %%
def intersectHelper(segment, node:list):
    '''
        Intersection helper to find the same between a halfedge and a node (list of segments).
        @returns-> List of segments
    '''
    ret = []
    part1 = Segment2(segment.source().point(), segment.target().point())

    for segs in node:
        isx = skInter(part1, segs)
        if type(isx) == Segment2:
            ret.append(isx)

    return ret


# %%
def updateVisibility(queryPoint:Point2, Target_:list, arrange:arrangement, visReg:list = [],node = []):
    '''
        Plain visible region calculation
        @args:
            query point:Point2  --> query point
            Target_:[segments] --> Target node
            arrangement:Arrangement --> The arrangement (see docs)
            visReg:[segments] --> list of visible segments (after first pass) 
        @returns:
            visible region:[segments]
    '''
    arrNew = arrange

    if(node != []):
        nodeVx = node.exterior.coords
        nodeBuilder = []

        # [TODO] Build a duplicate checker since skgeom throws an error
        for i in range(len(nodeVx)-1):
            nodeBuilder.append(Segment2(Point2(nodeVx[i][0],nodeVx[i][1]), Point2(nodeVx[i+1][0],nodeVx[i+1][1])))

        for bounds in nodeBuilder:
            arrNew.insert(bounds)


    vs = RotationalSweepVisibility(arrNew)
    face = arrNew.find(queryPoint)
    vx = vs.compute_visibility(queryPoint, face)
    vis = []

    # Draw visibility --optional, should be removed in prod
    # for he in arr.halfedges:
    #     draw.draw(he.curve(), visible_point=False)
    # draw.draw(queryPoint, color='magenta')

    if(visReg == []):
        for v in vx.halfedges:
            vis+=intersectHelper(v,Target_)
    else:
        # print("else!")
        for v in vx.halfedges:
            vis+=intersectHelper(v,visReg)

    ret = []
    [ret.append(x) for x in vis if x not in ret]
    # for v in ret:
    #     draw.draw(v, color='red', visible_point=False)
    return ret


# %%
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


# %%
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
            dist = round(dist,6)
    # print(dist)
    return dist


# %%
def minDist(queryPoint:Point,node:Polygon):
    nPoints = nearest_points(queryPoint,node)
    return round(nPoints[0].distance(nPoints[1]),6)


# %%
def getScaledNode(node):

    if(node.area == 0.0):
        return Point(node.x*10000,node.y*10000)

    vertices = node.exterior.coords
    scaledVx = []
    for vertex in vertices:
        scaledVx.append((vertex[0]*10000,vertex[1]*10000))
    scaledNode = Polygon(scaledVx)
    return scaledNode


# %%
def getScaledPoint(index, qdf):
    return Point(qdf['geometry'][index].x*10000,qdf['geometry'][index].y*10000)


# %%
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


# %%
def insideVisibleRegion(visibleRegionSet:list, node:Polygon):
    '''
        Test if VR and (skgeom)polygon intersect
            @args:
                visibilityRegionSet:[Polygons] --> VR from query point
                node:[segments] --> obstacle
            @returns:
                bool --> True when polygon and VR intersect, False o.w.
    '''
    for vrPolygon in visibleRegionSet:
        intersection = node.intersection(vrPolygon)
        # print(intersection.area)
        if intersection.area > 0.0:
            return True
    
    return False


# %%
def enqueue(x:{},key,value):
    '''
        Enqueue in priority queue
        @args:
            priority queue:dict
            key
            value
        @returns:
            priority queue:dict
    '''
    x[key]=-1*value
    return dict(sorted(x.items(), key=lambda item: item[1]))

def dequeue(x:{}):
    '''
        Dequeue last element from priority queue
        @args:
            priority queue:dict
        @returns:
            tuple(key,value)
    '''
    if len(x) != 0:
        element = x.popitem()
    else:
        element = (-1,-1)
    return element




# %%
# FOR QV
def getImpactRatio(vrOld:list, vrNew:list, queryPoint:Point):
    lengthOld = getLength(vrOld,queryPoint)
    lengthNew = getLength(vrNew,queryPoint)

    return (1 - round(lengthNew/lengthOld,6))


# %%
def plotFinal(scaledTarget,rtree1,gdf,qpGdf,vrPolygons,ans_vrPoly,ans):
    '''
        Plotting routine
        Plot scaled Target: Polygon
        Plot scaled Obstacles: Polygons
        Plot scaled Query Points: Points
        Plot retrieved Query Points: Points
        Plot query point visibility polygons: Polygons (already scaled)

        Append everything to a geoseries and plot it.
    '''

    Tlist = []
    Tlist.append(scaledTarget)
    Obslist = []
    for nodeList in rtree1.leaves():
        for nodeId in nodeList[1]:
            if isinstance(gdf.iloc[nodeId].geometry,Polygon):
                Obslist.append(getScaledNode(gdf.iloc[nodeId].geometry))
    QPlist = []
    for ix in qpGdf.index:
        QPlist.append(getScaledPoint(ix,qpGdf))
    VRlist = []
    for ix1 in vrPolygons:
        if ix1:
            for ix2 in ix1:
                VRlist.append(ix2)
    ansVRlist = []
    for ix1 in ans_vrPoly:
        if ix1:
            for ix2 in ix1:
                ansVRlist.append(ix2)
    Alist = []
    for ix in ans:
        Alist.append(getScaledPoint(ix,qpGdf))

    Tplot = gp.GeoSeries(Tlist)
    QPplot = gp.GeoSeries(QPlist)
    Obsplot = gp.GeoSeries(Obslist)
    VRplot = gp.GeoSeries(VRlist)
    ansVRplot = gp.GeoSeries(ansVRlist)
    Aplot = gp.GeoSeries(Alist)

    base = VRplot.plot(color='pink', edgecolor='black',alpha=0.2,figsize=(20,20))
    base2 = QPplot.plot(ax = base, color='black',markersize=3,figsize=(20,20))
    base3 = Obsplot.plot(ax = base2, color='cyan', edgecolor='black',figsize=(20,20))
    base4 = Tplot.plot(ax = base3, color='orange',figsize=(20,20))
    base5 = Aplot.plot(ax = base4, color='red',markersize=3,figsize=(20,20))
    base6 = ansVRplot.plot(ax = base5, color='lime', edgecolor='crimson',alpha=0.25,figsize=(20,20))


# %%



