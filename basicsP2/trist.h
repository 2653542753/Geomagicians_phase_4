#ifndef TRISTH
#define TRISTH

#include "pointSetArray.h"
#include <vector>
#include <iostream>
#include <set>
#include "..\basics\li.h"

using namespace std;

/*

  For a triangle abc, if version 0 is abc

  version 0 abc     (v:012)
  version 1 bca		(v:120)
  version 2 cab		(v:201)
  version 3 bac		(v:102)
  version 4 cba		(v:210)
  version 5 acb		(v:021)

  enext cycles   0 > 1 > 2 > 0 ,  5 > 4 > 3 > 5
  sym cycles  i <> (i + 3) % 6

  org  = ver < 3 ? v[ver] : v[(ver+1)%3]
  dest = ver < 3 ? v[(ver+1)%3] : v[ver-3] 

*/



typedef  int OrTri;  // The OrTri data structure for a triangle
typedef  int FIndex; // The index of a triangle Hint: NOT a triangle if it's negative
                     // You should be able to make all the triangle indices to be from 0 to n - 1 (n = number of triangles)

class TriangulateState {
public:
	int step;
	vector<int> linkPoints;
	std::vector<int> bigTriangle;

	bool isDone() { return step == -1; };

	TriangulateState() { this->step = -1; };
};




class Trist;

class TriRecord {
protected:
		int vi_[3];
		OrTri fnext_[6];
		int triIdx;
		bool visible;
		friend Trist;
		friend class TriangulateByEdgeCDT;

public:
	bool getVisibility(){ return visible; }
	int* getVertices(){ return vi_; }
	int getIdx(){ return triIdx; }
};

class Edge{
public:
	int vi_[2];
protected:
	friend Trist;
};

class Trist {
	private : 
		PointSetArray pointSet;
		vector<TriRecord> records;
		int maxTriIdx;
		TriRecord* findTriangle(int tIdx);
		vector<int> bigTriangle;
 
		set<int> pointsOnTri;
		set<int> lastPointsAddedIdx;
		set<int> lastTrianglesAddedIdx;
		set<int> activeEdge;
		int activePoint;
		bool flipping;

	protected:
		int en_[6];

	public:
		vector<int> addedPoints;
		vector<int> addedEdges;
		vector<Edge> constrainedEdges;
          
		Trist();
		int addPoint(LongInt x, LongInt y);
		int getPoint (int pIndex, LongInt& x1,LongInt& y1); // put the x,y values into x1,y1, and return 1 if the point pIndex exists
		int noPt();                                         // return the number of points

		int noTri(); // return the number of triangles
		int makeTri(int pIndex1,int pIndex2,int pIndex3,bool autoMerge = false); // Add a triangle into the Trist with the three point indices
		// Moreover, automatically establish the fnext pointers to its neigbhours if autoMerge = true

		void delTri(OrTri ef); // Delete a triangle, but you can assume that this is ONLY used by the IP operation
		                    // You may want to make sure all its neighbours are detached (below)
		
		void make3Tri(LongInt x, LongInt y);
		vector<int> make3Tri(int pIdx); //the point already exists

		OrTri enext(OrTri ef);
		OrTri sym(OrTri ef);
		OrTri fnext(OrTri ef);

		void getVertexIdx(OrTri ef, int& pIdx1,int& pIdx2,int& pIdx3); // return the three indices of the three vertices by OrTri
		int org(OrTri ef);  // the index of the first vertex of OrTri, e.g. org(bcd) => b
		int dest(OrTri ef); // the index of the second vertex of OrTri, e.g. org(bcd) => c
		void fmerge(OrTri abc, OrTri abd); // glue two neighbouring triangles, result abd = fnext(abc)
		void fdetach(OrTri abc); // detach triangle abc with all its neighbours (undo fmerge)

		//void incidentTriangles(int ptIndex,int& noOrTri, OrTri* otList); // A suggested function: you may want this function to return all the OrTri
		                                                                 // that are incident to this point
		                                                                 // Ignore this if you don't feel a need
		OrTri inTriangle(int ptIndex); //if returns -1, we are not in any triangle (ptIndex is not the vertex of any triangle)
		vector<int> adjacentTriangles(int pIdx1, int pIdx2); //indexes of triangles wich have pIdx1,pIdx2 as edge
		vector<int> adjacentTriangles(int pIdx); //indexes of triangles wic<h have pIdx as vertex
		
		std::vector< std::pair<OrTri, int> > findNeighbours(TriRecord tri);

		bool isLocallyDelaunay(int pIdx1, int pIdx2);
		vector<int> flipEdge(int pIdx1, int pIdx2);
		void flippingAlg(int pIdx1, int pIdx2);
		void triangulate(); //we assume there is no triangle
		TriangulateState *triangulateByPoint(int pIdx);
		bool triangulateByPointStep(TriangulateState *state);
		void hideBigTriangle();
		void showBigTriangles();
		void addPointUpdate(LongInt x, LongInt y);
		std::vector<int> getTriIdx();
		void setVisibility(int triIdx, bool visibility);
		std::vector<TriRecord> getTriangles();
		std::vector<MyPoint> getPoints();

		void clearActive() { activeEdge.clear();
						     lastPointsAddedIdx.clear();
							 lastTrianglesAddedIdx.clear(); 
							 activePoint = -1; };
		void setActiveEdge(int pIdx1, int pIdx2) { activeEdge.clear(); activeEdge.insert(pIdx1); activeEdge.insert(pIdx2); };
		void clearActiveEdge() { activeEdge.clear(); };
		bool isPointOnTri(int pIdx) { return pointsOnTri.count(pIdx) == 1; };
		bool isPointLastAdded(int pIdx) { return lastPointsAddedIdx.count(pIdx) == 1; };
		bool isTriangleLastAdded(int tIdx) { return lastTrianglesAddedIdx.count(tIdx) == 1; };
		bool isActiveEdge(int pIdx1, int pIdx2) { return activeEdge.count(pIdx1) == 1 && activeEdge.count(pIdx2) == 1; };
		bool isActivePoint(int pIdx) { return activePoint == pIdx; };

		int addConstrainedEdge(int pIdx1, int pIdx2);
		void triangulateByPointCDT(int pIdx);
		void triangulateByEdgeCDT(int edgeIdx);
		void retriangulateCDT(int a, int b, vector<int> pIdx);
		void flippingAlgCDT(int pIdx1, int pIdx2);
		bool constrained(int pIdx1, int pIdx2);
		void triangulateCDT();

		friend class TriangulateCDT;
		friend class FlippingAlgCDT;
		friend class TriangulateByPointCDT;
		friend class TriangulateByEdgeCDT;
		friend class ReTriangulate;
};

class ReTriangulate
{
protected:
	int step;
	bool done;
	int a, b, c;
	int posC;
	vector<int> pIdx;
	vector<int> leftVec, rightVec;
	Trist *triangle;

	ReTriangulate *rt;
public:
	ReTriangulate(Trist *tri, int p1, int p2, vector<int> pVec){
		a = p1;
		b = p2;
		pIdx = pVec;
		rt = NULL;
		done = false;
		triangle = tri;
		step = 0;
	}

	void next()
	{
		if (done)
			return;

		int pos;

		switch (step){
		case 0:
			//find the point c such that no point of Pidx is in circumcircle abc
			if (pIdx.size() == 0){
				done = true;
				return;
			}

			c = pIdx.front();
			pos = 0;
			posC = 0;
			for (vector<int>::iterator itV = pIdx.begin(); itV != pIdx.end(); itV++){
				if (triangle->pointSet.inCircle(a, b, c, *itV) == 1){
					c = *itV;
					posC = pos;
				}
				pos++;
			}

			triangle->makeTri(a, b, c, true);

			leftVec = vector<int>(pIdx.begin(), pIdx.begin() + posC);
			rightVec = vector<int>(pIdx.begin() + posC + 1, pIdx.end());

			if (leftVec.size() == 0 && rightVec.size() == 0){
				done = true;
				return;
			}

			step = 1;
			return;
		case 1:
			if (leftVec.size() != 0){
				if (rt == NULL){
					rt = new ReTriangulate(triangle, a, c, leftVec);
				}
				rt->next();
				if (rt->isDone()){
					delete(rt);
					rt = NULL;
					step = 2;
					if (rightVec.size() == 0){
						done = true;
					}
				}
				return;
			}
			step = 2;
		case 2:
			if (rightVec.size() != 0){
				if (rt == NULL){
					rt = new ReTriangulate(triangle, c, b, rightVec);
				}
				rt->next();
				if (rt->isDone()){
					delete(rt);
					rt = NULL;
					done = true;
				}
				return;
			}
		}

		done = true;
	}

	bool isDone()
	{
		return done;
	}
};

class TriangulateByEdgeCDT
{
protected:
	int step;
	bool done;
	int edgeIdx;
	int a, b;
	Trist *triangle;
	vector<int> linkPoints;
	vector<int> triToRemove, upperPts, lowerPts;

	ReTriangulate *rt;
public:
	TriangulateByEdgeCDT(Trist *tri, int eIdx){
		done = false;
		triangle = tri;
		edgeIdx = eIdx;
		step = 0;
		rt = NULL;
	}

	void next(){
		if (done)
			return;

		vector<int> adjTriangles;
		vector<int>::iterator it;
		TriRecord triRec;
		int c(-1), d(-1), e(-1);

		switch (step){
		case 0:
			a = triangle->constrainedEdges.at(edgeIdx).vi_[0];
			b = triangle->constrainedEdges.at(edgeIdx).vi_[1];
			//we are finding t1
			adjTriangles = triangle->adjacentTriangles(a);
			if (adjTriangles.empty()){
				done = true;
				return;
			}

			for (it = adjTriangles.begin(); it != adjTriangles.end(); ++it){
				triRec = *triangle->findTriangle(*it);
				if (a == triRec.vi_[0]){
					c = triRec.vi_[1];
					d = triRec.vi_[2];
				}
				else if (a == triRec.vi_[1]){
					c = triRec.vi_[0];
					d = triRec.vi_[2];
				}
				else{
					c = triRec.vi_[0];
					d = triRec.vi_[1];
				}
				//case where ab is already part of the triangulation: nothing to do
				if (c == b || d == b){
					done = true;
					return;
				}
				if (triangle->pointSet.intersects(a, b, c, d) == 1){
					triToRemove.push_back(*it);
					if (triangle->pointSet.turnLeft(a, b, c) == 1){
						upperPts.push_back(c);
						lowerPts.push_back(d);
					}
					else{
						upperPts.push_back(d);
						lowerPts.push_back(c);
					}
					break;
				}
			}
			//we are finding the other triangles which are crossing ab
			while (true){
				adjTriangles.clear();
				adjTriangles = triangle->adjacentTriangles(c, d);
				int tri;

				if (adjTriangles.at(0) != triToRemove.back()){
					tri = adjTriangles.at(0);
				}
				else{
					tri = adjTriangles.at(1);
				}
				triToRemove.push_back(tri);

				triRec = *triangle->findTriangle(tri);
				if ((triRec.vi_[0] == c && triRec.vi_[1] == d) || (triRec.vi_[0] == d && triRec.vi_[1] == c)){
					e = triRec.vi_[2];
				}
				else if ((triRec.vi_[1] == c && triRec.vi_[2] == d) || (triRec.vi_[1] == d && triRec.vi_[2] == c)){
					e = triRec.vi_[0];
				}
				else{
					e = triRec.vi_[1];
				}

				if (e == b){
					break;
				}

				if (triangle->pointSet.intersects(a, b, e, c) == 1){
					d = e;
				}
				else {
					c = e;
				}

				if (triangle->pointSet.turnLeft(a, b, e) == 1){
					upperPts.push_back(e);
				}
				else{
					lowerPts.push_back(e);
				}
			}

			step = 1;
			return;
		case 1:
			if (triToRemove.size() > 0){
				//we remove the triangles
				for (it = triToRemove.begin(); it != triToRemove.end(); ++it){
					triangle->delTri(*it << 3);
				}
				step = 2;
				return;
			}
			else{
				done = true;
				return;
			}
		case 2:
			//re-triangulate with respect to edge ab
			if (lowerPts.size() != 0){
				if (rt == NULL){
					rt = new ReTriangulate(triangle, a, b, lowerPts);
				}
				rt->next();
				if (rt->isDone()){
					delete(rt);
					rt = NULL;
					step = 3;
					if (upperPts.size() == 0){
						done = true;
					}
				}
				return;
			}
			step = 3;
		case 3:
			if (upperPts.size() != 0){
				if (rt == NULL){
					rt = new ReTriangulate(triangle, a, b, upperPts);
				}
				rt->next();
				if (rt->isDone()){
					delete(rt);
					rt = NULL;
					done = true;
				}
				return;
			}
		}

		done = true;
	}

	bool isDone()
	{
		return done;
	}
};

class FlippingAlgCDT
{
protected:
	bool done;
	int step;
	int pIdx1, pIdx2;
	vector<int> points;
	Trist *triangle;

	FlippingAlgCDT *facdt;
public:
	FlippingAlgCDT(Trist *tri, int p1, int p2){
		step = 0;
		done = false;
		triangle = tri;
		facdt = NULL;
		pIdx1 = p1;
		pIdx2 = p2;
	}

	void next()
	{
		if (done)
			return;

		switch (step){
		case 0:
			triangle->setActiveEdge(pIdx1, pIdx2);
			if (!triangle->constrained(pIdx1, pIdx2) && !triangle->isLocallyDelaunay(pIdx1, pIdx2)){
				points = triangle->flipEdge(pIdx1, pIdx2);
				triangle->setActiveEdge(points.at(0), points.at(1));
				step = 1;
				return;
			}
			else{
				done = true;
				return;
			}
		case 1:
			if (facdt == NULL){
				facdt = new FlippingAlgCDT(triangle, pIdx1, points.at(0));
			}
			facdt->next();
			if (facdt->isDone()){
				delete(facdt);
				facdt = NULL;
				step++;
				return;
			}
			return;
		case 2:
			if (facdt == NULL){
				facdt = new FlippingAlgCDT(triangle, pIdx1, points.at(1));
			}
			facdt->next();
			if (facdt->isDone()){
				delete(facdt);
				facdt = NULL;
				step++;
				return;
			}
			return;
		case 3:
			if (facdt == NULL){
				facdt = new FlippingAlgCDT(triangle, pIdx2, points.at(0));
			}
			facdt->next();
			if (facdt->isDone()){
				delete(facdt);
				facdt = NULL;
				step++;
				return;
			}
			return;
		case 4:
			if (facdt == NULL){
				facdt = new FlippingAlgCDT(triangle, pIdx2, points.at(1));
			}
			facdt->next();
			if (facdt->isDone()){
				delete(facdt);
				facdt = NULL;
				done = true;
			}
			return;
		}

		done = true;
	}

	bool isDone()
	{
		return done;
	}
};

class TriangulateByPointCDT
{
protected:
	int step;
	bool done;
	int pIdx;
	Trist *triangle;
	vector<int> linkPoints;

	FlippingAlgCDT *facdt;
public:
	TriangulateByPointCDT(Trist *tri, int pI){
		done = false;
		triangle = tri;
		pIdx = pI;
		step = 0;
		facdt = NULL;
	}
	void next()
	{
		if (done)
			return;

		switch (step){
		case 0:
			if (triangle->bigTriangle.empty()){
				triangle->bigTriangle = triangle->pointSet.constructCircumTri();
				triangle->makeTri(triangle->bigTriangle.at(0), triangle->bigTriangle.at(1), triangle->bigTriangle.at(2));
				step = 1;
				return;
			}
		case 1:
			linkPoints = triangle->make3Tri(pIdx);
			triangle->activePoint = pIdx;
			if (linkPoints.empty()){
				cout << "Point not in triangulation" << endl;
				done = true;
				return;
			}
			step = 2;
			return;
		case 2:
			cout << "Point in triangulation" << endl;
			triangle->setActiveEdge(linkPoints.at(0), linkPoints.at(1));
			step++;
			return;
		case 3:
			if (facdt == NULL){
				facdt = new FlippingAlgCDT(triangle, linkPoints.at(0), linkPoints.at(1));
			}
			facdt->next();
			if (facdt->isDone()){
				delete(facdt);
				facdt = NULL;
				step++;
				return;
			}
			return;
		case 4:
			triangle->setActiveEdge(linkPoints.at(0), linkPoints.at(2));
			step++;
			return;
		case 5:
			if (facdt == NULL){
				facdt = new FlippingAlgCDT(triangle, linkPoints.at(0), linkPoints.at(2));
			}
			facdt->next();
			if (facdt->isDone()){
				delete(facdt);
				facdt = NULL;
				step++;
				return;
			}
			return;
		case 6:
			triangle->setActiveEdge(linkPoints.at(1), linkPoints.at(2));
			step++;
			return;
		case 7:
			if (facdt == NULL){
				facdt = new FlippingAlgCDT(triangle, linkPoints.at(1), linkPoints.at(2));
			}
			facdt->next();
			if (facdt->isDone()){
				delete(facdt);
				facdt = NULL;
				step++;
				return;
			}
			return;
		case 8:
			triangle->clearActiveEdge();
		}
		done = true;
	}
	bool isDone()
	{
		return done;
	}
};


class TriangulateCDT
{
protected:
	bool done;
	Trist *triangle;
	vector<int>::iterator points_it;
	vector<int>::iterator edges_it;

	TriangulateByPointCDT *tpcdt;
	TriangulateByEdgeCDT *tecdt;
public:
	TriangulateCDT(Trist *tri){
		done = false;
		triangle = tri;
		tpcdt = NULL;
		tecdt = NULL;
		points_it = triangle->addedPoints.begin();
		edges_it = triangle->addedEdges.begin();
	}
	void next()
	{
		if (done)
			return;

		if (points_it != triangle->addedPoints.end()){
			if (tpcdt == NULL){
				tpcdt = new TriangulateByPointCDT(triangle, *points_it);
			}
			tpcdt->next();
			if (tpcdt->isDone()){
				delete(tpcdt);
				tpcdt = NULL;
				points_it++;
				if (points_it == triangle->addedPoints.end() && edges_it == triangle->addedEdges.end())
					done = true;
			}
			return;
		}

		triangle->clearActive();
		if (edges_it != triangle->addedEdges.end()){
			if (tecdt == NULL){
				tecdt = new TriangulateByEdgeCDT(triangle, *edges_it);
			}
			tecdt->next();
			if (tecdt->isDone()){
				delete(tecdt);
				tecdt = NULL;
				edges_it++;
			}
			else
				return;
		}

		triangle->addedPoints.clear();
		triangle->addedEdges.clear();
		done = true;
	}
	bool isDone()
	{
		return done;
	}
};


#endif
