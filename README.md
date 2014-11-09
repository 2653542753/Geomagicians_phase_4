Geomagicians_phase_4
====================
Here is the Project Phase 4 submitted by

Shi Xudong A0109682Y Tang Binbin A0105559B Croze Elias A0119984M Toh Zijing A0123506R	

==================Input Format====================
IP for adding points
CE for constraint edges
DT for normal Delaunay Triangulation
CDT for constraint Delaunay Triangulation
=================Algorithm Description for DT====
1.Local Delaunay Test
Furthermore, let e be an edge in a triangulation T in the plane. If e is an edge of fewer than two triangles in T, then e is said to be locally Delaunay. If e is an edge of exactly two triangles t1 and t2 in T, then e is said to be locally Delaunay if it has an open ball containing no vertex of t1 or t2. If e is not locally Delaunay, then e needs to be flipped around.
For every edges inside the edge set, we first find its adjacent triangles. For example, for edge ab, we find the adjacent triangles abc and abd. Then we run in-circle test on abc and abc. Edge ab will be flipped if they are not local Delaunay. After flipping, the adjacent edges will be also checked by local Delaunay test.

2.Incremental Triangulation
Incremental triangulation construction mainly first construct a large enough triangle p0p1p2 to contain all the given points. For every point p, we will find out the 3 new edge to form an edge set running local Delaunay test. After finishing all edge tests, the large triangle will be removed.
=================Algorithm Description for CDT===
http://www.cescg.org/CESCG-2004/web/Domiter-Vid/index.html#algorithm


