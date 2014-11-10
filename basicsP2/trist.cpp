#include "trist.h"
#include <iostream>
#include <algorithm>
#include <gl/glut.h>

using namespace std;

Trist::Trist()
{
	en_[0] = 1;
	en_[1] = 2;
	en_[2] = 0;
	en_[3] = 5;
	en_[4] = 3;
	en_[5] = 4;
	maxTriIdx = 0;
	this->pointSet = PointSetArray();
	this->records = vector<TriRecord>();
	this->activePoint = -1;
	this->flipping = false;
}

int Trist::addPoint(LongInt x, LongInt y){
	int idx = this->pointSet.addPoint(x, y);
	lastPointsAddedIdx.insert(idx);
        addedPoints.push_back(idx);
	return idx;
}

int Trist::getPoint (int pIndex, LongInt& x1,LongInt& y1){
	return this->pointSet.getPoint(pIndex, x1, y1);
}

int Trist::noPt(){
	return this->pointSet.noPt();
}

int Trist::noTri(){
	return this->records.size();
}

int Trist::makeTri(int pIndex1,int pIndex2,int pIndex3,bool autoMerge){
	vector<pair<OrTri,int> > neigh;
	TriRecord tri;

	lastPointsAddedIdx.clear();
	pointsOnTri.insert(pIndex1);
	pointsOnTri.insert(pIndex2);
	pointsOnTri.insert(pIndex3);

	tri.visible = true;
	maxTriIdx++;
	int newTriIdx = maxTriIdx;//modifi? avant size+1 (les indices commencent ?0)
	int tLeft = this->pointSet.turnLeft(pIndex1, pIndex2, pIndex3);
	for(int i=0; i<6; i++){
		tri.fnext_[i]=-1;
	}
	
	if(tLeft==1){
		tri.vi_[0]=pIndex1;
		tri.vi_[1]=pIndex2;
		tri.vi_[2]=pIndex3;
	}else if(tLeft==-1){
		tri.vi_[0]=pIndex2;
		tri.vi_[1]=pIndex1;
		tri.vi_[2]=pIndex3;
	}else{
		cout << "Not a triangle" << endl;
		return -1;
	}
	tri.triIdx = newTriIdx;
	this->records.push_back(tri);

	if(autoMerge){
		//find the neighbours of tri
		neigh = this->findNeighbours(tri);
		for(int i=0; i<neigh.size(); i++){
			pair<OrTri,int> p = neigh.at(i);
			this->fmerge(p.first, (tri.triIdx <<3) | p.second);
		}
	}

	lastTrianglesAddedIdx.insert(tri.triIdx);
	return tri.triIdx;
}

void Trist::delTri(OrTri ef){
	this->fdetach(ef);
	int tIdx = ef >> 3;
	OrTri tri2;
	int tIdx2, ver2;
	
	for(vector<TriRecord>::iterator it = records.begin(); it != records.end(); ++it){
		if(it->triIdx == tIdx){
			records.erase(it);
			break;
		}
	}
/*
	int size = this->noTri();

	//we have to update the indexes (all triangles after tIdx have a index decreased by 1)
	for(int i=0; i<size-1; i++){
		for(int j=0; j<6; j++){
			tri2 = this->records.at(i).fnext_[j];
			if(tri2!=-1){
				tIdx2 = tri2 >> 3;
				ver2 = tri2 & 7;
				if(tIdx2 > tIdx){
					tIdx2--;
					this->records.at(i).fnext_[j] = (tIdx2 << 3) | ver2;
				}
			}
		}
	}
*/	
}
OrTri Trist::enext(OrTri ef){
	int idx = ef >> 3;
	int ver = ef & 7;
	return (idx << 3) | en_[ver];
}

void Trist::make3Tri(LongInt x, LongInt y){
	int p1Idx = -1;
	int p2Idx = -1;
	int p3Idx = -1;
	int pIdx = this->addPoint(x,y);
	OrTri tri = this->inTriangle(pIdx);

	lastTrianglesAddedIdx.clear();

	if(tri >= 0){
		this->getVertexIdx(tri, p1Idx, p2Idx, p3Idx);
		this->delTri(tri);

		this->makeTri(pIdx, p1Idx, p2Idx, true);
		this->makeTri(pIdx, p2Idx, p3Idx, true);
		this->makeTri(pIdx, p3Idx, p1Idx, true);
	}
}

vector<int> Trist::make3Tri(int pIdx){
	int p1Idx = -1;
	int p2Idx = -1;
	int p3Idx = -1;
	OrTri tri = this->inTriangle(pIdx);
	vector<int> linkPoints;

	lastTrianglesAddedIdx.clear();

	if(tri >= 0){
		this->getVertexIdx(tri, p1Idx, p2Idx, p3Idx);
		this->delTri(tri);

		this->makeTri(pIdx, p1Idx, p2Idx, true);
		this->makeTri(pIdx, p2Idx, p3Idx, true);
		this->makeTri(pIdx, p3Idx, p1Idx, true);
		linkPoints.push_back(p1Idx);
		linkPoints.push_back(p2Idx);
		linkPoints.push_back(p3Idx);
	}
	return linkPoints;
}

OrTri Trist::sym(OrTri ef){
	int idx = ef >> 3;
	int ver = ef & 7;
	return (idx << 3) | ((ver + 3) % 6);
}

OrTri Trist::fnext(OrTri ef){
	int idx = ef >> 3;
	int ver = ef & 7;
	TriRecord* tri = findTriangle(idx);
	if(tri == NULL){
		return -1;
	}
	return tri->fnext_[ver];
}

void Trist::getVertexIdx(OrTri ef, int& pIdx1,int& pIdx2,int& pIdx3){
	int tIdx = ef >> 3;
	const int ver = ef & 7;
	TriRecord *tri = findTriangle(tIdx);

	if(tri == NULL){
		pIdx1 = -1;
		pIdx2 = -1;
		pIdx3 = -1;
		return;
	}

	int p0 = tri->vi_[0];
	int p1 = tri->vi_[1];
	int p2 = tri->vi_[2];

	switch(ver){
	case 0:
		pIdx1 = p0;
		pIdx2 = p1;
		pIdx3 = p2;
		break;
	case 1:
		pIdx1 = p1;
		pIdx2 = p2;
		pIdx3 = p0;
		break;
	case 2:
		pIdx1 = p2;
		pIdx2 = p0;
		pIdx3 = p1;
		break;
	case 3:
		pIdx1 = p1;
		pIdx2 = p0;
		pIdx3 = p2;
		break;
	case 4:
		pIdx1 = p2;
		pIdx2 = p1;
		pIdx3 = p0;
		break;
	case 5:
		pIdx1 = p0;
		pIdx2 = p2;
		pIdx3 = p1;
		break;
	default:
		break;
	}
}

int Trist::org(OrTri ef){
	int tIdx = ef >> 3;
	const int ver = ef & 7;
	TriRecord *tri = findTriangle(tIdx);

	if(tri == NULL){
		return -1;
	}

	int p0 = tri->vi_[0];
	int p1 = tri->vi_[1];
	int p2 = tri->vi_[2];

	switch(ver){
	case 0:
		return p0;
	case 1:
		return p1;
	case 2:
		return p2;
	case 3:
		return p1;
	case 4:
		return p2;
	case 5:
		return p0;
	default:
		return p0;
	}
}

int Trist::dest(OrTri ef){
	int tIdx = ef >> 3;
	const int ver = ef & 7;

	TriRecord *tri = findTriangle(tIdx);

	if(tri == NULL){
		return -1;
	}

	int p0 = tri->vi_[0];
	int p1 = tri->vi_[1];
	int p2 = tri->vi_[2];
	
	switch(ver){
	case 0:
		return p1;
	case 1:
		return p2;
	case 2:
		return p0;
	case 3:
		return p0;
	case 4:
		return p1;
	case 5:
		return p2;
	default:
		return p1;
	}
}

void Trist::fdetach(OrTri abc){
	int tIdx = abc >> 3;
	OrTri tri2;
	int tIdx2, ver2;

	TriRecord *tri = findTriangle(tIdx);
	if(tri == NULL){
		return;
	}

	for(int i=0; i<6; i++){
		tri2 = tri->fnext_[i];
		if(tri2!=-1){
			tIdx2 = tri2 >> 3;
			ver2 = tri2 & 7;
			tri->fnext_[i] = -1;
			TriRecord *triRec2 = findTriangle(tIdx2);
			if(triRec2 != NULL){
				triRec2->fnext_[ver2] = -1;
			}
		}
	}
}

void Trist::fmerge(OrTri abc, OrTri abd){
	TriRecord *t1,*t2;
	
	int tIdx1 = abc >> 3;
	int ver1 = abc & 7;
	int tIdx2 = abd >> 3;
	int ver2 = abd & 7;

	t1 = findTriangle(tIdx1);
	t2 = findTriangle(tIdx2);

	if(t1 == NULL || t2 == NULL){
		return;
	}

	t1->fnext_[ver1]=abd;
	t2->fnext_[ver2]=abc;

	//fnext(bac)=bad
	abc = this->sym(abc);
	abd = this->sym(abd);
	ver1 = abc & 7;
	ver2 = abd & 7;
	t1->fnext_[ver1]=abd;
	t2->fnext_[ver2]=abc;
}

OrTri Trist::inTriangle(int ptIndex){
	int inTriangle;

	for(int i=0; i<(this->noTri()); i++){
		inTriangle = this->pointSet.inTri(this->records.at(i).vi_[0],this->records.at(i).vi_[1],this->records.at(i).vi_[2],ptIndex);
		if(inTriangle == 1){
			return records.at(i).triIdx << 3; //version 0 of triangle
		}
	}

	return -1; //we didn't find any triangle
}

//in pair<OrTri, int>, OrTri is the triangle neighbouring tri1, and int is the corresponding version of tri1
vector<pair<OrTri, int> > Trist::findNeighbours(TriRecord tri1){
	vector<pair<OrTri,int>> neigh;
	int noT = this->noTri();
	TriRecord tri2;

	//we go through every triangle
	for(int i=0; i<noT; i++){
		int count = 0;
		int commonVert1[3] = {-1,-1,-1};//indexes of common vertices in tri1
		int commonVert2[3] = {-1,-1,-1};//indexes of common vertices in tri2
		int version1 = -1;
		int version2 = -1;

		tri2 = this->records.at(i);

		//we go through every vertex
		for(int j=0; j<3; j++){
			for(int k=0; k<3; k++){
				if(tri2.vi_[j]==tri1.vi_[k]){
					commonVert1[count]=k;
					commonVert2[count]=j;
					count++;
				}
				if(count==3){
					break;
				}
			}
		}
		if(count==2){
			if(commonVert1[0]==0 && commonVert1[1]==1){
				version1=0;
			}else if(commonVert1[0]==1 && commonVert1[1]==2){
				version1=1;
			}else if(commonVert1[0]==2 && commonVert1[1]==0){
				version1=2;
			}else if(commonVert1[0]==1 && commonVert1[1]==0){
				version1=3;
			}else if(commonVert1[0]==2 && commonVert1[1]==1){
				version1=4;
			}else{
				version1=5;
			}

			if(commonVert2[0]==0 && commonVert2[1]==1){
				version2=0;
			}else if(commonVert2[0]==1 && commonVert2[1]==2){
				version2=1;
			}else if(commonVert2[0]==2 && commonVert2[1]==0){
				version2=2;
			}else if(commonVert2[0]==1 && commonVert2[1]==0){
				version2=3;
			}else if(commonVert2[0]==2 && commonVert2[1]==1){
				version2=4;
			}else{
				version2=5;
			}

			neigh.push_back(pair<OrTri,int>((tri2.triIdx<<3)|version2,version1));
		}
	}
	return neigh;
}

vector<int> Trist::adjacentTriangles(int pIdx1, int pIdx2){
	vector<int> adj;
	vector<TriRecord>::iterator it;
	int count;
	for(int i =0; i < records.size(); i++){
		count = 0;
		TriRecord r = records.at(i);
		if(r.vi_[0]==pIdx1 || r.vi_[1]==pIdx1 || r.vi_[2]==pIdx1){
			count++;
		}
		if(r.vi_[0]==pIdx2 || r.vi_[1]==pIdx2 || r.vi_[2]==pIdx2){
			count++;
		}
		if(count==2){
			adj.push_back(r.triIdx);
		}
	}
	return adj;
}

bool Trist::isLocallyDelaunay(int pIdx1, int pIdx2){
	vector<int> adj = adjacentTriangles(pIdx1, pIdx2);
	int a, b, c, d, e, f;
	vector<int> notDuplicates;
	if(adj.size()==2){
		getVertexIdx(adj.at(0) << 3, a, b, c);
		getVertexIdx(adj.at(1) << 3, d, e, f);
		notDuplicates.push_back(pIdx1);
		notDuplicates.push_back(pIdx2);
		if(find(notDuplicates.begin(), notDuplicates.end(), a)==notDuplicates.end()){
			notDuplicates.push_back(a);
		}
		if(find(notDuplicates.begin(), notDuplicates.end(), b)==notDuplicates.end()){
			notDuplicates.push_back(b);
		}
		if(find(notDuplicates.begin(), notDuplicates.end(), c)==notDuplicates.end()){
			notDuplicates.push_back(c);
		}
		if(find(notDuplicates.begin(), notDuplicates.end(), d)==notDuplicates.end()){
			notDuplicates.push_back(d);
		}
		if(find(notDuplicates.begin(), notDuplicates.end(), e)==notDuplicates.end()){
			notDuplicates.push_back(e);
		}
		if(find(notDuplicates.begin(), notDuplicates.end(), f)==notDuplicates.end()){
			notDuplicates.push_back(f);
		}
		return pointSet.inCircle(notDuplicates.at(0), notDuplicates.at(1), notDuplicates.at(2), notDuplicates.at(3))<0;
	}
	return true;
}

vector<int> Trist::flipEdge(int pIdx1, int pIdx2){
	vector<int> adj = adjacentTriangles(pIdx1, pIdx2);
	vector<int> notDuplicates;
	int a, b, c, d, e, f;
	if(adj.size()==2){
		getVertexIdx(adj.at(0) << 3, a, b, c);
		getVertexIdx(adj.at(1) << 3, d, e, f);
		notDuplicates.push_back(pIdx1);
		notDuplicates.push_back(pIdx2);
		if(find(notDuplicates.begin(), notDuplicates.end(), a)==notDuplicates.end()){
			notDuplicates.push_back(a);
		}
		if(find(notDuplicates.begin(), notDuplicates.end(), b)==notDuplicates.end()){
			notDuplicates.push_back(b);
		}
		if(find(notDuplicates.begin(), notDuplicates.end(), c)==notDuplicates.end()){
			notDuplicates.push_back(c);
		}
		if(find(notDuplicates.begin(), notDuplicates.end(), d)==notDuplicates.end()){
			notDuplicates.push_back(d);
		}
		if(find(notDuplicates.begin(), notDuplicates.end(), e)==notDuplicates.end()){
			notDuplicates.push_back(e);
		}
		if(find(notDuplicates.begin(), notDuplicates.end(), f)==notDuplicates.end()){
			notDuplicates.push_back(f);
		}
		notDuplicates.erase(std::remove(notDuplicates.begin(), notDuplicates.end(), pIdx1), notDuplicates.end());
		notDuplicates.erase(std::remove(notDuplicates.begin(), notDuplicates.end(), pIdx2), notDuplicates.end());
		delTri(adj.at(0) << 3);
		delTri(adj.at(1) << 3);
		makeTri(notDuplicates.at(0),notDuplicates.at(1), pIdx1, true);
		makeTri(notDuplicates.at(0),notDuplicates.at(1), pIdx2, true);
	}
	return notDuplicates;
}

void Trist::flippingAlg(int pIdx1, int pIdx2){
	if(!isLocallyDelaunay(pIdx1, pIdx2)){
		cout << "flip " << pIdx1 << "," << pIdx2 << endl;
		vector<int> points = flipEdge(pIdx1, pIdx2);
		setActiveEdge(points.at(0), points.at(1));

		flippingAlg(pIdx1, points.at(0));
		flippingAlg(pIdx1, points.at(1));
		flippingAlg(pIdx2, points.at(0));
		flippingAlg(pIdx2, points.at(1));
	}
}

void Trist::triangulate(){
	int maxIndx = noPt()-1;
	bigTriangle = pointSet.constructCircumTri();
	vector<int> p = bigTriangle;
	vector<int> linkPoints;
	vector<int> triToDel;
	makeTri(p.at(0), p.at(1), p.at(2));
	for(int i=0; i<= maxIndx; i++){
		linkPoints = make3Tri(i);
		flippingAlg(linkPoints.at(0), linkPoints.at(1));
		flippingAlg(linkPoints.at(0), linkPoints.at(2));
		flippingAlg(linkPoints.at(1), linkPoints.at(2));
	}

	for (int idx = 0; idx < 3; idx++) {
		triToDel = adjacentTriangles(p.at(idx));
		for (vector<int>::iterator it = triToDel.begin(); it != triToDel.end(); ++it) {
			setVisibility(*it, false);
		}
	}
}

void Trist::addPointUpdate(LongInt x, LongInt y){
	//We assume the last point 
	int pIdx = addPoint(x,y);
	vector<int> linkPoints = make3Tri(pIdx);
	if(linkPoints.empty()){
		//erase all triangles
		//triangulate
		//cout << "Point not in triangulation" << endl;
		records.clear();
		triangulate();
	}else{
		//cout << "Point in triangulation" << endl;
		flippingAlg(linkPoints.at(0), linkPoints.at(1));
		flippingAlg(linkPoints.at(0), linkPoints.at(2));
		flippingAlg(linkPoints.at(1), linkPoints.at(2));

		vector<int> triToDel;
		for (int idx = 0; idx < 3; idx++) {
			triToDel = adjacentTriangles(bigTriangle.at(idx));
			for (vector<int>::iterator it = triToDel.begin(); it != triToDel.end(); ++it) {
				setVisibility(*it, false);
			}
		}
	}
}

TriangulateState *Trist::triangulateByPoint(int pIdx){
	TriangulateState *state = new TriangulateState();

	if(bigTriangle.empty()){
		bigTriangle = pointSet.constructCircumTri();
		makeTri(bigTriangle.at(0), bigTriangle.at(1), bigTriangle.at(2));
	}
	vector<int> linkPoints = make3Tri(pIdx);
	activePoint = pIdx;
	if(linkPoints.empty()){
		//erase all triangles
		//triangulate
		//cout << "Point not in triangulation" << endl;
		records.clear();
		triangulate();
	}else{
		state->step = 0;
		state->linkPoints = linkPoints;
		state->bigTriangle = bigTriangle;
	}
	return state;
}

bool Trist::triangulateByPointStep(TriangulateState *state) {
	switch (state->step){
	case 0:
		setActiveEdge(state->linkPoints.at(0), state->linkPoints.at(1));
		break;
	case 1:
		flippingAlg(state->linkPoints.at(0), state->linkPoints.at(1));
		break;
	case 2:
		setActiveEdge(state->linkPoints.at(0), state->linkPoints.at(2));
		break;
	case 3:
		flippingAlg(state->linkPoints.at(0), state->linkPoints.at(2));
		break;
	case 4:
		setActiveEdge(state->linkPoints.at(1), state->linkPoints.at(2));
		break;
	case 5:
		flippingAlg(state->linkPoints.at(1), state->linkPoints.at(2));
		break;
	case 6:
		clearActiveEdge();
		state->step = -2;
		break;
	}
	state->step++;
	return state->isDone();
}

void Trist::hideBigTriangle(){
	vector<int> triToDel;
	for (int idx = 0; idx < 3; idx++) {
		triToDel = adjacentTriangles(bigTriangle.at(idx));
		for (vector<int>::iterator it = triToDel.begin(); it != triToDel.end(); ++it) {
			setVisibility(*it, false);
		}
	}
}

vector<int> Trist::adjacentTriangles(int pIdx){
	vector<int> tri;
	for(int i =0; i < records.size(); i++){
		TriRecord r = records.at(i);
		if(r.vi_[0]==pIdx || r.vi_[1]==pIdx || r.vi_[2]==pIdx){
			tri.push_back(r.triIdx);
		}
	}
	return tri;
}

TriRecord* Trist::findTriangle(int tIdx){
	TriRecord* tri;
	vector<TriRecord>::iterator it = records.begin();
	for(; it != records.end(); ++it){
		if(it->triIdx == tIdx){
			break;
		}
	}

	if(it == records.end()){
		tri = NULL;
	}else{
		tri = &(*it);
	}

	return tri;
}

vector<int> Trist::getTriIdx(){
	vector<int> triIdx;
	for(vector<TriRecord>::iterator it = records.begin(); it != records.end(); ++it){
		triIdx.push_back(it->triIdx);
	}
	return triIdx;
}

void Trist::setVisibility(int triIdx, bool visibility){
	TriRecord *tri = findTriangle(triIdx);
	if(tri != NULL){
		tri->visible = visibility;
	}
}

std::vector<TriRecord> Trist::getTriangles(){
	return records;
}

std::vector<MyPoint> Trist::getPoints(){
	return pointSet.getPoints();
}

void Trist::showBigTriangles(){
	vector<int> triToDel;
	for (int idx = 0; idx < 3; idx++) {
		triToDel = adjacentTriangles(bigTriangle.at(idx));
		for (vector<int>::iterator it = triToDel.begin(); it != triToDel.end(); ++it) {
			setVisibility(*it, true);
		}
	}
}

int Trist::addConstrainedEdge(int pIdx1, int pIdx2){
	Edge e;
	e.vi_[0] = pIdx1;
	e.vi_[1] = pIdx2;
	if (pIdx1 < noPt() && pIdx2 < noPt()){
		constrainedEdges.push_back(e);
		addedEdges.push_back(constrainedEdges.size() - 1);
		return constrainedEdges.size() - 1;
	}
	return -1;
}

void Trist::triangulateByPointCDT(int pIdx){
	if (bigTriangle.empty()){
		bigTriangle = pointSet.constructCircumTri();
		makeTri(bigTriangle.at(0), bigTriangle.at(1), bigTriangle.at(2));
	}
	vector<int> linkPoints = make3Tri(pIdx);
	if (linkPoints.empty()){
		cout << "Point not in triangulation" << endl;
	}
	else{
		cout << "Point in triangulation" << endl;
		flippingAlgCDT(linkPoints.at(0), linkPoints.at(1));
		flippingAlgCDT(linkPoints.at(0), linkPoints.at(2));
		flippingAlgCDT(linkPoints.at(1), linkPoints.at(2));

		vector<int> triToDel;
		triToDel = adjacentTriangles(bigTriangle.at(0));
		for (vector<int>::iterator it = triToDel.begin(); it != triToDel.end(); ++it) {
			setVisibility(*it, false);
		}
		triToDel = adjacentTriangles(bigTriangle.at(1));
		for (vector<int>::iterator it = triToDel.begin(); it != triToDel.end(); ++it) {
			setVisibility(*it, false);
		}
		triToDel = adjacentTriangles(bigTriangle.at(2));
		for (vector<int>::iterator it = triToDel.begin(); it != triToDel.end(); ++it) {
			setVisibility(*it, false);
		}
	}
}

void Trist::flippingAlgCDT(int pIdx1, int pIdx2){
	if (!constrained(pIdx1, pIdx2) && !isLocallyDelaunay(pIdx1, pIdx2)){
		vector<int> points = flipEdge(pIdx1, pIdx2);
		flippingAlgCDT(pIdx1, points.at(0));
		flippingAlgCDT(pIdx1, points.at(1));
		flippingAlgCDT(pIdx2, points.at(0));
		flippingAlgCDT(pIdx2, points.at(1));
	}
}

bool Trist::constrained(int pIdx1, int pIdx2){
	vector<Edge>::iterator it;
	for (it = constrainedEdges.begin(); it != constrainedEdges.end(); ++it){
		Edge e = *it;
		if ((e.vi_[0] == pIdx1 && e.vi_[1] == pIdx2) || (e.vi_[0] == pIdx2 && e.vi_[1] == pIdx1)){
			return true;
		}
	}
	return false;
}

void Trist::triangulateByEdgeCDT(int edgeIdx){
	vector<int> triToRemove, upperPts, lowerPts;
	int a = constrainedEdges.at(edgeIdx).vi_[0];
	int b = constrainedEdges.at(edgeIdx).vi_[1];
	int c(-1), d(-1), e(-1);
	TriRecord triRec;

	//we are finding t1
	vector<int> adjTriangles = adjacentTriangles(a);
	if (adjTriangles.empty()){
		return;
	}

	vector<int>::iterator it;
	for (it = adjTriangles.begin(); it != adjTriangles.end(); ++it){
		triRec = *findTriangle(*it);
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
			return;
		}
		if (pointSet.intersects(a, b, c, d) == 1){
			triToRemove.push_back(*it);
			if (pointSet.turnLeft(a, b, c) == 1){
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
		adjTriangles = adjacentTriangles(c, d);
		int tri;

		if (adjTriangles.at(0) != triToRemove.back()){
			tri = adjTriangles.at(0);
		}
		else{
			tri = adjTriangles.at(1);
		}
		triToRemove.push_back(tri);

		triRec = *findTriangle(tri);
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

		if (pointSet.intersects(a, b, e, c) == 1){
			d = e;
		}
		else {
			c = e;
		}

		if (pointSet.turnLeft(a, b, e) == 1){
			upperPts.push_back(e);
		}
		else{
			lowerPts.push_back(e);
		}
	}

	//we remove the triangles
	for (it = triToRemove.begin(); it != triToRemove.end(); ++it){
		delTri(*it << 3);
	}

	//re-triangulate with respect to edge ab
	retriangulateCDT(a, b, lowerPts);
	retriangulateCDT(a, b, upperPts);
}

void Trist::triangulateCDT(){
        vector<int>::iterator it;
        for (it = addedPoints.begin(); it != addedPoints.end(); it++){
                triangulateByPointCDT(*it);
        }
        for (it = addedEdges.begin(); it != addedEdges.end(); it++){
                triangulateByEdgeCDT(*it);
        }
        addedPoints.clear();
        addedEdges.clear();
}

void Trist::retriangulateCDT(int a, int b, vector<int> pIdx){
	//find the point c such that no point of Pidx is in circumcircle abc
	if (pIdx.size() == 0){
		return;
	}

	int c = pIdx.front();
	int pos = 0;
	int posC = 0;
	for (vector<int>::iterator itV = pIdx.begin(); itV != pIdx.end(); itV++){
		if (pointSet.inCircle(a, b, c, *itV) == 1){
			c = *itV;
			posC = pos;
		}
		pos++;
	}

	makeTri(a, b, c, true);
	retriangulateCDT(a, c, vector<int>(pIdx.begin(), pIdx.begin() + posC));
	retriangulateCDT(c, b, vector<int>(pIdx.begin() + posC + 1, pIdx.end()));
}
