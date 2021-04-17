#include "stdafx.h"
#include "gz.h"
#include "dda.h"

/*
-- Sort three verts by Y from small to large
*/
void Sort(GzCoord vertex[3][3]) {
	//makes A<B
	if(vertex[0][0][1] > vertex[1][0][1])
		Switch(vertex[0], vertex[1]);

	//makes B<C, C is maximum y
	if (vertex[1][0][1] > vertex[2][0][1])
		Switch(vertex[1], vertex[2]);

	//makes A<B, A is minmim y
	if (vertex[0][0][1] > vertex[1][0][1])
		Switch(vertex[0], vertex[1]);
}

/*
-- Switch two vert (position/normal/uv coordinate)
*/
void Switch(GzCoord vert1[3], GzCoord vert2[3]) {

	GzCoord tempVertex[3];

	for (int vertexInfor = 0; vertexInfor < 3; vertexInfor++) {
		for (int coord = 0; coord < 3; coord++) {
			tempVertex[vertexInfor][coord] = vert1[vertexInfor][coord];			// temp  = vert1
			vert1[vertexInfor][coord] = vert2[vertexInfor][coord];				// vert1  = vert2
			vert2[vertexInfor][coord] = tempVertex[vertexInfor][coord];			// vert2  = temp
		}
	}
}

/*
--Create an edge in DDA format with two verts coord as input
*/
Edge SetUpEdge(GzCoord vertA, GzCoord vertB) {
	Edge edge;

	edge.start[0] = vertA[0]; edge.start[1] = vertA[1]; edge.start[2] = vertA[2];
	edge.end[0] = vertB[0]; edge.end[1] = vertB[1]; edge.end[2] = vertB[2];
	edge.current[0] = vertA[0]; edge.current[1] = vertA[1]; edge.current[2] = vertA[2];
	edge.slope_x = (vertB[0] - vertA[0]) / (vertB[1] - vertA[1]);
	edge.slope_z = (vertB[2] - vertA[2]) / (vertB[1] - vertA[1]);

	return edge;
}
/*
--Create an span using  two coord from left to right
*/
Span SetUpSpan(GzCoord coordL, GzCoord coordR) {
	Span span;

	float x = coordL[0];
	span.start[0] = coordL[0]; span.start[1] = coordL[1]; span.start[2] = coordL[2];
	span.end[0] = coordR[0]; span.end[1] = coordR[1]; span.end[2] = coordR[2];
	span.current[0] = coordL[0]; span.current[1] = coordL[1]; span.current[2] = coordL[2];
	span.slope_z = (coordR[2] - coordL[2]) / (coordR[0] - coordL[0]);

	return span;
}
/*
--Update DDA edge's current position based on change on Y
*/
void UpdateDDA(Edge* edge, float deltaY) {

	(*edge).current[0] = (*edge).current[0] + (*edge).slope_x * deltaY;
	(*edge).current[1] = (*edge).current[1] + deltaY;
	(*edge).current[2] = (*edge).current[2] + (*edge).slope_z * deltaY;

}
/*
--Update  Span's current position based on change on X
*/
void UpdateSpan(Span* span, float deltaX) {
	(*span).current[0] = (*span).current[0] + deltaX;
	(*span).current[2] = (*span).current[2] + deltaX * (*span).slope_z;

}
