#ifndef DDA_H
#define DDA_H

#include "gz.h"
#include "stdafx.h"

/*struct of an edge in DDA Formate*/
#ifndef Edge_H
#define Edge_H
struct Edge {
	GzCoord start;
	GzCoord end;
	GzCoord current;
	float slope_x;
	float slope_z;
}; 
#endif

/*struct of span in DDA Formate*/
#ifndef Span_H
#define Span_H
struct Span {
	GzCoord start;
	GzCoord end;
	GzCoord current;
	float slope_z;
};

#endif
void Sort(GzCoord[3][3]);						/* sort three vertex by y coord */
void Switch(GzCoord[3], GzCoord[3]);			/* switch two vertex*/
Edge SetUpEdge(GzCoord, GzCoord);				/* return a DDA edge from two verts been input */
void UpdateDDA(Edge*, float);					/* update DDA's current value based on deltaY */
Span SetUpSpan(GzCoord, GzCoord);				/* return a span based on GzCoord of two edge's current position */
void UpdateSpan(Span*, float);					/* update span's current valuev based on deltaX */
#endif