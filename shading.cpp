#include "stdafx.h"
#include "shading.h"
#include "rend.h"
#include "matrixOperation.h"


/*
-- check if normal is on the same side of L and E
*/
float CorrectNormal(GzCoord N, GzCoord L, GzCoord E) {

	float NdotL = DotProduct(N, L);
	float NdotE = DotProduct(N, E);
	int indicator;//0: 1, 1:-1, 2, 0
	if (NdotL > 0 && NdotE > 0) {				// L,E,and N is on the same side
		indicator = 1;
	}
	else if (NdotL < 0 && NdotE < 0) {			// L,E are on the opposite side of N
		indicator = -1;
	}
	else {										// L,E are not on the same side
		indicator = 0;
	}

	return indicator;

}

/*
-- Calculate specular light
-- Specular light = Ks sumof[Ie (R.E)s] 
*/
void SpecularLight(GzCoord specular, GzColor Ks, GzLight* directLights, int numlights, GzCoord normal, GzCoord E, float spec) {
	
	GzCoord R;

	for (int lightNum = 0; lightNum < numlights; lightNum++) {

		int orientation = CorrectNormal(normal, directLights[lightNum].direction, E);		//check orientation
		if (orientation > 0) {
			R[0] = normal[0];
			R[1] = normal[1];
			R[2] = normal[2];																//R = N
		}
		else if (orientation < 0) {
			R[0] = -normal[0];
			R[1] = -normal[1];
			R[2] = -normal[2];																//R = -N
		}
		else {
			continue;
		}

		VectorScalar(R, 2 * DotProduct(normal, directLights[lightNum].direction));			//R =  2(N.L)N
		VectorSubstract(R, directLights[lightNum].direction);								//R = 2(N.L)N - L

		GzCoord lightColor = { directLights[lightNum].color[0],directLights[lightNum].color[1],directLights[lightNum].color[2] };
		float dotProductRE = DotProduct(R, E);												//R.E
		if (dotProductRE < 0) {
			dotProductRE = 0;																//clamp R.E to zero
		}
		VectorScalar(lightColor, pow(dotProductRE, spec));									//current color = le*(R.E)^SPEC
		VectorAddition(specular, lightColor);												//specular +=  current color

	}

	specular[0] *= Ks[0];																	//specular = Ks*specular 
	specular[1] *= Ks[1];
	specular[2] *= Ks[2];

	for (int i = 0; i < 3; i++) {															//clamp specular light color
		if (specular[i] < 0)
			specular[i] = 0;
		if (specular[i] > 1)
			specular[i] = 1;
	}

}


/*
-- Calculate diffuse light color
-- DiffuseLight = Kd sumof[Ie (N.L)]) 
*/
void DiffuseLight(GzCoord diffuse, GzColor Kd, GzLight* directLights, int numlights, GzCoord normal, GzCoord E) {

	for (int lightNum = 0; lightNum < numlights; lightNum++) {

		int orientation = CorrectNormal(normal, directLights[lightNum].direction, E);	//check orientation
		float NdotL;																		
		if (orientation > 0) {
			NdotL = DotProduct(normal, directLights[lightNum].direction);					//N.L
		}
		else if (orientation < 0) {
			NdotL = -DotProduct(normal, directLights[lightNum].direction);					//-N.L
		}
		else {
			continue;
		}

		if (NdotL < 0) {
			NdotL = 0;																		//clamp N.L to zero
		}

		GzCoord lightColor = { directLights[lightNum].color[0],directLights[lightNum].color[1],directLights[lightNum].color[2] };
		VectorScalar(lightColor, NdotL);													//current color = lightcolor * N.L 
		VectorAddition(diffuse, lightColor);												//diffuse +=  current color
	}

	diffuse[0] *= Kd[0];																	//diffuse = Ks*diffuse 
	diffuse[1] *= Kd[1];
	diffuse[2] *= Kd[2];

	for (int i = 0; i < 3; i++) {															//clamp diffuse color
		if (diffuse[i] < 0)
			diffuse[i] = 0;
		if (diffuse[i] > 1)
			diffuse[i] = 1;
	}

}

/*
-- Calculate Ambinet light
-- AmbientLight = Ka Ia
*/
void AmbientLight(GzCoord ambient, GzColor Ka, GzLight ambientlight) {

	ambient[0] = Ka[0] * ambientlight.color[0];												//ambient = Ka*ambient light color
	ambient[1] = Ka[1] * ambientlight.color[1];
	ambient[2] = Ka[2] * ambientlight.color[2];

	for (int i = 0; i < 3; i++) {															//clamp ambient color
		if (ambient[i] < 0)
			ambient[i] = 0;
		if (ambient[i] > 1)
			ambient[i] = 1;
	}

}



/*
-- Calcuate shaded color
*/
GzColor* LightingEq(GzColor Ks, GzColor Kd, GzColor Ka, GzLight* directLights, int numlights, GzLight ambientlight, GzCoord& normal, float spec) {
	
	GzColor result = { 0,0,0 };													//set default value	
	GzCoord E = { 0,0,-1 };
	GzCoord specular = { 0,0,0 };
	GzCoord diffuse = { 0,0,0 };
	GzCoord ambient = { 0,0,0 };
	Normalization(normal);														//normalize normal

	//Specular Light
	SpecularLight(specular, Ks, directLights, numlights, normal, E, spec);

	//diffuse Light
	DiffuseLight(diffuse, Kd,directLights, numlights,normal, E);

	//ambient Light
	AmbientLight(ambient, Ka, ambientlight);
	

	//combine 
	VectorAddition(result, specular);
	VectorAddition(result, diffuse);
	VectorAddition(result, ambient);

	//check overflow
	for (int i = 0; i < 3; i++) {
		if (result[i] < 0)
			result[i] = 0;
		if (result[i] > 1)
			result[i] = 1;
	}

	return &result;
}


/*
-- using 3 (vertex x,y and third variable) Inputs, to the interpolated result at certain location
*/
float Interpolation(GzCoord vertCoord1, GzCoord vertCoord2, GzCoord vertCoord3, float* thirdParam, float x, float y) {
	GzCoord edge1, edge2;
	GzCoord crossproduct;

	edge1[0] = vertCoord1[0] - vertCoord2[0]; edge1[1] = vertCoord1[1] - vertCoord2[1]; edge1[2] = thirdParam[0] - thirdParam[1];
	edge2[0] = vertCoord1[0] - vertCoord3[0]; edge2[1] = vertCoord1[1] - vertCoord3[1]; edge2[2] = thirdParam[0] - thirdParam[2];
	CrossProduct(edge1, edge2, crossproduct);

	float A = crossproduct[0];
	float B = crossproduct[1];
	float C = crossproduct[2];
	float D = -(vertCoord1[0] * crossproduct[0] + vertCoord1[1] * crossproduct[1] + thirdParam[0] * crossproduct[2]);

	return (A * x + B * y + D) / -C;
}

/*
-- Get interpolated normal at (i,j) using three vertices' coordinate and normal
*/
void NormalInterpolation(GzCoord interpolatedNormal, GzCoord Vertex[3][3], float i, float j) {

	float normalX[3] = { Vertex[0][1][0],Vertex[1][1][0] ,Vertex[2][1][0] };
	float normalY[3] = { Vertex[0][1][1],Vertex[1][1][1] ,Vertex[2][1][1] };
	float normalZ[3] = { Vertex[0][1][2],Vertex[1][1][2] ,Vertex[2][1][2] };


	interpolatedNormal[0] = Interpolation(Vertex[0][0], Vertex[1][0], Vertex[2][0], normalX, i, j);			//interpolated x
	interpolatedNormal[1] = Interpolation(Vertex[0][0], Vertex[1][0], Vertex[2][0], normalY, i, j);			//interpolated y
	interpolatedNormal[2] = Interpolation(Vertex[0][0], Vertex[1][0], Vertex[2][0], normalZ, i, j);			//interpolated z
}

/*
-- Get Interpolated color at (i,j) using three vertices' coordinate and color
*/
void ColorInterpolation(GzColor interpolatedColor, GzCoord Vertex[3][3], GzColor vert1Color, GzColor vert2Color, GzColor vert3Color, float i, float j) {

	float ColorR[3] = { vert1Color[0],vert2Color[0],vert3Color[0] };
	float ColorG[3] = { vert1Color[1],vert2Color[1],vert3Color[1] };
	float ColorB[3] = { vert1Color[2],vert2Color[2],vert3Color[2] };

	interpolatedColor[0] = Interpolation(Vertex[0][0], Vertex[1][0], Vertex[2][0], ColorR, i, j);			//interpolated R
	interpolatedColor[1] = Interpolation(Vertex[0][0], Vertex[1][0], Vertex[2][0], ColorG, i, j);			//interpolated G
	interpolatedColor[2] = Interpolation(Vertex[0][0], Vertex[1][0], Vertex[2][0], ColorB, i, j);			//interpolated B
}

/*
-- Get a perspective Corrected u,v texture coordinate at (x,y)
*/
void PInterpolation(float* u, float* v, GzCoord Vertex[3][3], float x, float y) {

	float v_prime_z;							// new term V'z
	float v_s_z[3];								// Vzs
	float u_p_s[3];								// U value in perspective space
	float v_p_s[3];								// V value in perspective space
	float interpolated_U, interpolated_V, interpolated_z;
	

	// three vertices u,v -> prespective space U,V
	for (int i = 0; i < 3; i++) {
		v_s_z[i] = Vertex[i][0][2];							// z of vertex at perspective space 
		v_prime_z = v_s_z[i] / (MAXINT - v_s_z[i]);			// temp value V's for interpolation
		u_p_s[i] = Vertex[i][2][0] / (v_prime_z + 1);		// U in perspective space
		v_p_s[i] = Vertex[i][2][1] / (v_prime_z + 1);		// V in perspective space
	}

	//interpolate: get U,V at x,y

	interpolated_U = Interpolation(Vertex[0][0], Vertex[1][0], Vertex[2][0], u_p_s, x, y);		// get interpolated U
	interpolated_V = Interpolation(Vertex[0][0], Vertex[1][0], Vertex[2][0], v_p_s, x, y);		// get interpolated V
	interpolated_z = Interpolation(Vertex[0][0], Vertex[1][0], Vertex[2][0], v_s_z, x, y);		// get interpolated z

	// U,V -> affine space u,v
	v_prime_z = interpolated_z / (MAXINT - interpolated_z);								// new term V'z in perspective space
	*u = interpolated_U * (v_prime_z + 1);												// new u value in affine space
	*v = interpolated_V * (v_prime_z + 1);												// new v value in affine space

}

