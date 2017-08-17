# Shortest-Path-Finding-on-Terrain-Surface-or-3D-Model

The main source code is main.cpp. There are three different functions in main.cpp, namely shortestpath_MMP, shortestpath_GB and shortestpath_LA. The first function is an implementation of Reference 1 and the second function is an implementation of Reference 2. The last function is a landmark-based shortest path algorithm. 

# Compiling Command

g++4 -o main main.cpp -std=c++11

# How to Use It

main "terrain_data"

The 2nd parameter is the terrain data file. 

An example: 

main small_terrain.off

In this example, samll_terrain.off is the terrain data file.

Data Format:

We used the .off format in the experiment. The content of the .off file is as follows:

OFF

Number_of_vertices Number_of_faces Number_of_edges

x_coordinate_of_1st_vertex y_coordinate_of_1st_vertex z_coordinate_of_1st_vertex

x_coordinate_of_2nd_vertex y_coordinate_of_2nd_vertex z_coordinate_of_2nd_vertex

......

x_coordinate_of_the_last_vertex y_coordinate_of_the_last_vertex z_coordinate_of_the_last_vertex

ID_of_the_1st_vertex_of_the_1st_face ID_of_the_2nd_vertex_of_the_1st_face ID_of_the_3td_vertex_of_the_1st_face

ID_of_the_1st_vertex_of_the_2nd_face ID_of_the_2nd_vertex_of_the_2nd_face ID_of_the_3td_vertex_of_the_2nd_face

......

ID_of_the_1st_vertex_of_the_last_face ID_of_the_2nd_vertex_of_the_last_face ID_of_the_3td_vertex_of_the_last_face

Each .off data could be visualized by the terrain tool (http://rwcpu1.cse.ust.hk/terrain/).

# Experimental Result

The program will save result in a file for each function. 

Each file contains "query distance" "query time (ms)" for each query performed (each row corresponds to a query and the first two columns corresponds to "query distance" and "query time (ms)" of the query). 

# References

1. Mitchell, Joseph SB, David M. Mount, and Christos H. Papadimitriou. "The discrete geodesic problem." SIAM Journal on Computing 16.4 (1987): 647-668.

2. V. Verma and J. Snoeyink: Reducing the memory required to
find a geodesic shortest path on a large mesh. ACM
SIGSPATIAL GIS 2009
