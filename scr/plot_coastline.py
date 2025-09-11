from circle_from_three_points import circle_from_three_points
# Import necessary modules
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import shapely as shp
import numpy as np
from shapely.geometry import MultiLineString, Point
from shapely.geometry import Polygon

show_vtx_old = True
show_circles = False
radius_tol   = 1e3    

# Set filepath
fp = "../data/Topo-bathy/shape_format/data/Kystlinie.shp"
# fp = "../data/Topo-bathy/shape_format/data/Hav_dybde_5m_2005.shp"
# fp = "../data/Topo-bathy/shape_format/data/Land_hojde_5m_2005.shp"
# fp = "../data/coastlines-split-3857/lines.shp"
# fp = "../data/coastlines-split-4326/lines.shp"

# Read file using gpd.read_file()
gdf = gpd.read_file(fp)
tol = 100000000


multiline = shp.geometry.MultiLineString(gdf.geometry.values)
multiline = shp.line_merge(multiline)

#multiline = multiline.simplify(tolerance=5e3)

polygons = []
for line in multiline.geoms:
    polygon = Polygon(line)
    area    = polygon.area
    if area>tol:
        polygons.append(polygon)

fig,ax = plt.subplots() 


polygons_new = []
vtx_old     = []
radii_old   = []
centers_old = []


for polygon in polygons:
    x,y = polygon.exterior.xy

    vtx  = np.array([x[:-1],y[:-1]])

    #for iter in range(100):

    nvtx = vtx.shape[1]

    radii   = np.zeros(nvtx)
    centers = np.zeros((2,nvtx))

    # First point:
    centers[:,0],radii[0],tmp = circle_from_three_points(vtx[:,-1],vtx[:, 0],vtx[:, 1])
    # Middle points:
    for i in range(1,nvtx-1):
        centers[:,i],radii[i],tmp = circle_from_three_points(vtx[:,i-1],vtx[:,i  ],vtx[:,i+1])
    # Last point:
    centers[:,-1],radii[-1],tmp = circle_from_three_points(vtx[:,-2],vtx[:,-1],vtx[:, 0])

    while True:
     
        radius_min = np.min(radii)
        radius_max = np.max(radii)

        if radius_min>radius_tol:
            print(f"nvtx {nvtx} radius_min {radius_min}")
            break

        min_index = np.argmin(radii)

        radii_old.append  (radii  [min_index])
        centers_old.append(centers[:,min_index])
        vtx_old.append    (vtx    [:,min_index])

        radii     = np.delete(radii,   min_index, axis=0)
        centers   = np.delete(centers, min_index, axis=1)
        vtx       = np.delete(vtx,     min_index, axis=1)


        if min_index==0:
            centers[:,-1],radii[-1],tmp = circle_from_three_points(vtx[:,-2],vtx[:,-1],vtx[:, 0])
            centers[:, 0],radii [0],tmp = circle_from_three_points(vtx[:,-1],vtx[:, 0],vtx[:, 1])
        elif min_index==nvtx-2:
            centers[:,-2],radii[-2],tmp = circle_from_three_points(vtx[:,-3],vtx[:,-2],vtx[:,-1])
            centers[:,-1],radii[-1],tmp = circle_from_three_points(vtx[:,-2],vtx[:,-1],vtx[:, 0])
        elif min_index==nvtx-1:
            centers[:,-1],radii[-1],tmp = circle_from_three_points(vtx[:,-2],vtx[:,-1],vtx[:, 0])
            centers[:, 0],radii [0],tmp = circle_from_three_points(vtx[:,-1],vtx[:, 0],vtx[:, 1])
        else:
            i = min_index-1
            centers[:,i],radii[i],tmp = circle_from_three_points(vtx[:,i-1],vtx[:,i  ],vtx[:,i+1])
            i = min_index
            centers[:,i],radii[i],tmp = circle_from_three_points(vtx[:,i-1],vtx[:,i  ],vtx[:,i+1])
    
        nvtx = vtx.shape[1]


    #vtx_old = np.array(vtx_old)

#
#for polygon in polygons:
#    x,y = polygon.exterior.xy
#    #plt.plot(x,y)
#
#for polygon in polygons_new:
#    x,y = polygon.exterior.xy
#    plt.plot(x,y)
#
#


    plt.plot(vtx[0,:],vtx[1,:],'b-')
    plt.plot(vtx[0,:],vtx[1,:],'g.')


    if show_vtx_old:
        plt.plot(np.array(vtx_old)[:,0],np.array(vtx_old)[:,1],'r.')

    if show_circles:
        for c,r in zip(centers_old,radii_old):
            ax.add_artist(plt.Circle(c, r,fill=False))

#for line in multiline.geoms:
#    x,y = line.xy
#    plt.plot(x,y)
plt.show()