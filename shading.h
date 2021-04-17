#pragma once
#ifndef SHADING_H
#define SHADING_H

#include "gz.h"
#include "stdafx.h"

float CorrectNormal(GzCoord, GzCoord, GzCoord);												/*correct normoral orientation*/

void SpecularLight(GzCoord, GzColor, GzLight*, int, GzCoord, GzCoord, float);				/*calculate specular light color*/
void DiffuseLight(GzCoord, GzColor, GzLight*, int, GzCoord, GzCoord);						/*calculate diffuse light color*/
void AmbientLight(GzCoord, GzColor, GzLight);												/*calculate ambient light color*/
GzColor* LightingEq(GzColor, GzColor, GzColor, GzLight*, int, GzLight, GzCoord&, float);	/*calcuate shaded color*/

float Interpolation(GzCoord, GzCoord, GzCoord, float*, float, float);						/*interpolate third variable by x,y*/
void NormalInterpolation(GzCoord, GzCoord[3][3], float, float);								/*get interpolated normal*/
void ColorInterpolation(GzColor, GzCoord[3][3], GzColor, GzColor, GzColor, float, float);	/*get interpolated color*/
void PInterpolation(float*, float*, GzCoord[3][3], float, float);							/*get prespective corrected uv*/

#endif
