#include "matrixOperation.h"
#include "stdafx.h"
#include "rend.h"
#include <math.h>

/*
-- normalze an input vector
-- vec =|vec|
*/
void Normalization(GzCoord vec) {
	float magnitude = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

	if (magnitude > 0) {
		vec[0] /= magnitude;
		vec[1] /= magnitude;
		vec[2] /= magnitude;
	}
}

/*
-- change a rotation matrix into unitary matrix
-- change a sccale matrix into identity matrix
-- change a translate matrix into identity matrix 
*/
void UnitaryMatrix(GzMatrix matrix) {
	float K = sqrt(matrix[0][0] * matrix[0][0] + matrix[0][1] * matrix[0][1] + matrix[0][2] * matrix[0][2]);
	for (int row = 0; row < 4; row++) {
		for (int column = 0; column < 4; column++) {
			matrix[row][column] = matrix[row][column]/K;
		}
	}
	matrix[0][3] = 0; matrix[1][3] = 0; matrix[2][3] = 0; matrix[3][3] = 1;
}


/*
-- return a dot product of two vectors
*/
float DotProduct(GzCoord vec1, GzCoord vec2) {
	float result = 0;

	result += vec1[0] * vec2[0];
	result += vec1[1] * vec2[1];
	result += vec1[2] * vec2[2];

	return result;
}

/*
-- return a corss product of two vectors
-- third vector is the result
*/
void CrossProduct(GzCoord vec1, GzCoord vec2 ,GzCoord result) {

	result[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
	result[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
	result[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];

}

/*
-- add second vector to the first one
*/
void VectorAddition(GzCoord vec1, GzCoord vec2){

	vec1[0] += vec2[0];
	vec1[1] += vec2[1];
	vec1[2] += vec2[2];
}

/*
-- scale the vector
*/
void VectorScalar(GzCoord vec, float scalar) {
	vec[0] *= scalar;
	vec[1] *= scalar;
	vec[2] *= scalar;
}

/*
-- vector substraction  vec1 = vec1 - vec2
*/
void VectorSubstract(GzCoord vec, GzCoord vec2) {
	vec[0] = vec[0] - vec2[0];
	vec[1] = vec[1] - vec2[1];
	vec[2] = vec[2] - vec2[2];
}

/*
-- vector multiple each element
-- vec1 = [e1,e2,e3] vec2 = [E1,E2,E3]
-- vec1 = [e1*E1, e2*E2, e3*E3]
*/
void VectorMultiplication(GzCoord vec1, GzCoord vec2) {
	vec1[0] *= vec2[0];
	vec1[1] *= vec2[1];
	vec1[2] *= vec2[2];
}

/*
-- reuturn multiplication of two matrix
*/
GzMatrix* MatrixMultiplication(GzMatrix MatrixA, GzMatrix MatrixB) {
	GzMatrix result;

	for (int row = 0; row < 4; row++) {
		for (int column = 0; column < 4; column++) {
			result[row][column] = 0;
		}
	}

	for (int row = 0; row < 4; row++) {
		for (int column = 0; column < 4; column++) {
			for (int r_c = 0; r_c < 4; r_c++) {
				result[row][column] += MatrixA[row][r_c] * MatrixB[r_c][column];
			}
		}
	}

	return &result;
}

/*
--  return a vector which is the multipliciation of a matrix and vector
--  matrix is a 4x4, and vector is 1x3
--  vector_result[1x4]  = matrix [4x4] * vector[x,y,z,1]
--  then cast this vector_result into [1x3] and put into parameter:vector
*/
int MatrixVectorMultiplication(GzMatrix matrix, GzCoord vector, bool checkZ) {
	GzCoord tempVector;
	
	tempVector[0] = matrix[0][0] * vector[0] + matrix[0][1] * vector[1] + matrix[0][2] * vector[2] + matrix[0][3] * 1;	//X before projection
	tempVector[1] = matrix[1][0] * vector[0] + matrix[1][1] * vector[1] + matrix[1][2] * vector[2] + matrix[1][3] * 1;	//Y before projection
	tempVector[2] = matrix[2][0] * vector[0] + matrix[2][1] * vector[1] + matrix[2][2] * vector[2] + matrix[2][3] * 1;	//Z before projection

	float w  = matrix[3][0] * vector[0] + matrix[3][1] * vector[1] + matrix[3][2] * vector[2] + matrix[3][3] * 1;			//w=(Z/d)+1

	if (checkZ && tempVector[2] < 0) {
		return -1;
	}

	vector[0] = tempVector[0] / w;
	vector[1] = tempVector[1] / w;
	vector[2] = tempVector[2] / w;

	return 0;
}

/*
--  clamp function between min and max
*/
float Clamp(float value, float min, float  max) {
	if (value < min)
		return min;
	if (value > max)
		return max;
	return value;
}