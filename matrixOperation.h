#pragma once
#ifndef MATRIXOPERATION_H
#define MATRIXOPERATION_H

#include "gz.h"
#include "stdafx.h"

void Normalization(GzCoord);								//normalize a vector
void UnitaryMatrix(GzMatrix);								//unitary a matrix
float DotProduct(GzCoord, GzCoord);							//return dot product of two vectors
void CrossProduct(GzCoord, GzCoord, GzCoord);				//cross product of two vectors
void VectorAddition(GzCoord, GzCoord);						//add second vector to the first one
void VectorScalar(GzCoord, float);							//scale vector
void VectorSubstract(GzCoord, GzCoord);						//vector substraction
void VectorMultiplication(GzCoord, GzCoord);				//two vector multiply element at same row
GzMatrix* MatrixMultiplication(GzMatrix, GzMatrix);			//return multiplication of two matrix
int MatrixVectorMultiplication(GzMatrix, GzCoord, bool);	//mulltplication of maxtrix and vector
#endif