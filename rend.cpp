/* CS580 Homework 3 */
#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"
#include	"dda.h"
#include	"matrixOperation.h"
#include	"shading.h"
#include	"Integrate.h"

#ifndef PI
#define PI (float) 3.14159265358979323846
#endif

// Zero vector
GzCoord ZeroCoord = { 0.0f, 0.0f, 0.0f };

// Multiply second vector to first vector
void ProductVec3(GzCoord Vec3, GzCoord Vec3b)
{
	Vec3[0] *= Vec3b[0];
	Vec3[1] *= Vec3b[1];
	Vec3[2] *= Vec3b[2];
}

// Cross Vector and store in result variable (must be passed in)
void CrossVec3(GzCoord Vec3, GzCoord Vec3b, GzCoord result)
{
	result[0] = Vec3[1] * Vec3b[2] - Vec3[2] * Vec3b[1];
	result[1] = Vec3[2] * Vec3b[0] - Vec3[0] * Vec3b[2];
	result[2] = Vec3[0] * Vec3b[1] - Vec3[1] * Vec3b[0];
}

// Set mat A equal to B
void SetMatricesEqual(GzMatrix MatA, GzMatrix MatB)
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			MatA[i][j] = MatB[i][j];
		}
	}
}

// Set CoordA equal to CoordB
void SetCoordEqual(GzCoord CoordA, GzCoord CoordB)
{
	for (int i = 0; i < 3; i++)
	{
		CoordA[i] = CoordB[i];
	}
}

// Calculate dot product of two vectors
float DotVec3(GzCoord Vec3, GzCoord Vec3b)
{
	float DotResult = 0.0f;
	DotResult += Vec3[0] * Vec3b[0];
	DotResult += Vec3[1] * Vec3b[1];
	DotResult += Vec3[2] * Vec3b[2];
	return DotResult;
}

// Normalize by dividing by magnitude
void NormalizeVector3(GzCoord Vec3)
{
	float magnitude = (float)sqrt(Vec3[0] * Vec3[0] + Vec3[1] * Vec3[1] + Vec3[2] * Vec3[2]);
	Vec3[0] /= magnitude;
	Vec3[1] /= magnitude;
	Vec3[2] /= magnitude;
}

// Normalize by dividing by magnitude
float GetMagnitudeVector3(GzCoord Vec3)
{
	float magnitude = (float)sqrt(Vec3[0] * Vec3[0] + Vec3[1] * Vec3[1] + Vec3[2] * Vec3[2]);
	return magnitude;
}

// Multiply vector with float value
void MultiplyVector3(GzCoord Vec3, float val)
{
	Vec3[0] *= val;
	Vec3[1] *= val;
	Vec3[2] *= val;
}

// Add two vectors, store result in first vector
void AddVector3(GzCoord Vec3, GzCoord Vec3b)
{
	Vec3[0] += Vec3b[0];
	Vec3[1] += Vec3b[1];
	Vec3[2] += Vec3b[2];
}

// Subtract two vectors, store result in first vector
void SubtractVector3(GzCoord Vec3, GzCoord Vec3b)
{
	Vec3[0] -= Vec3b[0];
	Vec3[1] -= Vec3b[1];
	Vec3[2] -= Vec3b[2];
}

int GzRender::GzRotXMat(float degree, GzMatrix mat)
{
/* HW 3.1
// Create rotate matrix : rotate along x axis
// Pass back the matrix using mat value
*/
	float rad = degree * PI / 180.0;
	mat[0][0] = 1; mat[0][1] = 0; mat[0][2] = 0; mat[0][3] = 0;
	mat[1][0] = 0; mat[1][1] = cos(rad); mat[1][2] = -sin(rad); mat[1][3] = 0;
	mat[2][0] = 0; mat[2][1] = sin(rad); mat[2][2] = cos(rad); mat[2][3] = 0;
	mat[3][0] = 0; mat[3][1] = 0; mat[3][2] = 0; mat[3][3] = 1;

	return GZ_SUCCESS;
}

int GzRender::GzRotYMat(float degree, GzMatrix mat)
{
/* HW 3.2
// Create rotate matrix : rotate along y axis
// Pass back the matrix using mat value
*/
	float rad = degree * PI / 180.0;
	mat[0][0] = cos(rad); mat[0][1] = 0; mat[0][2] = sin(rad); mat[0][3] = 0;
	mat[1][0] = 0; mat[1][1] = 1; mat[1][2] = 0; mat[1][3] = 0;
	mat[2][0] = -sin(rad); mat[2][1] = 0; mat[2][2] = cos(rad); mat[2][3] = 0;
	mat[3][0] = 0; mat[3][1] = 0; mat[3][2] = 0; mat[3][3] = 1;

	return GZ_SUCCESS;
}

int GzRender::GzRotZMat(float degree, GzMatrix mat)
{
/* HW 3.3
// Create rotate matrix : rotate along z axis
// Pass back the matrix using mat value
*/
	float rad = degree * PI / 180.0;
	mat[0][0] = cos(rad); mat[0][1] = -sin(rad); mat[0][2] = 0; mat[0][3] = 0;
	mat[1][0] = sin(rad); mat[1][1] = cos(rad); mat[1][2] = 0; mat[1][3] = 0;
	mat[2][0] = 0; mat[2][1] = 0; mat[2][2] = 1; mat[2][3] = 0;
	mat[3][0] = 0; mat[3][1] = 0; mat[3][2] = 0; mat[3][3] = 1;

	return GZ_SUCCESS;
}

int GzRender::GzTrxMat(GzCoord translate, GzMatrix mat)
{
/* HW 3.4
// Create translation matrix
// Pass back the matrix using mat value
*/
	mat[0][0] = 1; mat[0][1] = 0; mat[0][2] = 0; mat[0][3] = translate[0];
	mat[1][0] = 0; mat[1][1] = 1; mat[1][2] = 0; mat[1][3] = translate[1];
	mat[2][0] = 0; mat[2][1] = 0; mat[2][2] = 1; mat[2][3] = translate[2];
	mat[3][0] = 0; mat[3][1] = 0; mat[3][2] = 0; mat[3][3] = 1;

	return GZ_SUCCESS;
}


int GzRender::GzScaleMat(GzCoord scale, GzMatrix mat)
{
/* HW 3.5
// Create scaling matrix
// Pass back the matrix using mat value
*/
	mat[0][0] = scale[0]; mat[0][1] = 0; mat[0][2] = 0; mat[0][3] = 0;
	mat[1][0] = 0; mat[1][1] = scale[1]; mat[1][2] = 0; mat[1][3] = 0;
	mat[2][0] = 0; mat[2][1] = 0; mat[2][2] = scale[2]; mat[2][3] = 0;
	mat[3][0] = 0; mat[3][1] = 0; mat[3][2] = 0; mat[3][3] = 1;

	return GZ_SUCCESS;
}


GzRender::GzRender(int xRes, int yRes)
{
/* HW1.1 create a framebuffer for MS Windows display:
 -- set display resolution
 -- allocate memory for framebuffer : 3 bytes(b, g, r) x width x height
 -- allocate memory for pixel buffer
 */
	xres = (short)xRes;													//set resolution of x
	yres = (short)yRes;													//set resolution of y
	framebuffer = (char*)malloc(3 * sizeof(char) * xRes * yRes);		//allocate framebuffer
	pixelbuffer = new GzPixel[xRes * yRes];								//allocate  pixelbuffer for each x * y
	for (int i = 0; i < xRes*yRes; i++)
	{
		pixelbuffer[i].red = 0;
		pixelbuffer[i].green = 0;
		pixelbuffer[i].blue = 0;
	}
	
	AApixelbuffers = new GzPixel*[AAKERNEL_SIZE];									//AAbuffers contains 6 pixel buffer
	for(int i = 0; i < AAKERNEL_SIZE; i++) {
		AApixelbuffers[i] = new GzPixel[xRes * yRes];					//each buffer contain xRes * yRes pixels
	}

	matlevel = 0;														//default variables used for count
	numlights = 0;
	interp_mode = GZ_FLAT;

/* HW 3.6
- setup Xsp and anything only done once 
- init default camera 
*/ 
	//Xsp
	Xsp[0][0] = xres / 2; Xsp[0][1] = 0; Xsp[0][2] = 0; Xsp[0][3] = xres / 2;
	Xsp[1][0] = 0; Xsp[1][1] = -yres / 2; Xsp[1][2] = 0; Xsp[1][3] = yres / 2;
	Xsp[2][0] = 0; Xsp[2][1] = 0; Xsp[2][2] = MAXINT; Xsp[2][3] = 0;
	Xsp[3][0] = 0; Xsp[2][1] = 0; Xsp[3][2] = 0; Xsp[3][3] = 1;

	//default camera
	m_camera.position[X] = DEFAULT_IM_X; m_camera.position[Y] = DEFAULT_IM_Y; m_camera.position[Z] = DEFAULT_IM_Z;
	m_camera.lookat[X] = 0; m_camera.lookat[Y] = 0; m_camera.lookat[Z] = 0;
	m_camera.worldup[X] = 0; m_camera.worldup[Y] = 1; m_camera.worldup[Z] = 0;
	m_camera.FOV = DEFAULT_FOV;

}

GzRender::~GzRender()
{
/* HW1.2 clean up, free buffer memory */
	free(framebuffer);
	delete pixelbuffer;
	for (int i = 0; i < AAKERNEL_SIZE; i++) {
		delete AApixelbuffers[i];					//each buffer contain xRes * yRes pixels
	}
	for (int i = 0; i < NumSamplingPlanes; i++)
	{
		for (int z = 0; z < SamplingPlaneY; z++)
		{
			delete[] samplingPlanes[i].samplingPlanePixels[z];
		}
		delete[] samplingPlanes[i].samplingPlanePixels;
	}
	delete[] AApixelbuffers;
	delete[] samplingPlanes;
}

int GzRender::GzDefault()
{
/* HW1.3 set pixel buffer to some default values - start a new frame */
		for (int j = 0; j < yres; j++) {									//go through row start from 0  to yres
			for (int i = 0; i < xres; i++) {								//go through each element on row 	
				for (int AA = 0; AA < AAKERNEL_SIZE; AA++) {				//go through each AA pixelbuffers
					AApixelbuffers[AA][ARRAY(i, j)].red = 16;				//get a color I like, a = 1, z = max
					AApixelbuffers[AA][ARRAY(i, j)].green = 607;
					AApixelbuffers[AA][ARRAY(i, j)].blue = 1279;
					AApixelbuffers[AA][ARRAY(i, j)].alpha = 1;
					AApixelbuffers[AA][ARRAY(i, j)].z = INT_MAX;
				}
			}
		}
	

	return GZ_SUCCESS;
}

int GzRender::GzBeginRender()
{
/* HW 3.7 
- setup for start of each frame - init frame buffer color,alpha,z
- compute Xiw and projection xform Xpi from camera definition 
- init Ximage - put Xsp at base of stack, push on Xpi and Xiw 
- now stack contains Xsw and app can push model Xforms when needed 
*/ 
	GzDefault();														//init frame buffer

	//Xiw
	GzCoord camera_Xaxis, camera_Yaxis, camera_Zaxis;

	camera_Zaxis[0] = m_camera.lookat[X] - m_camera.position[X];		//cl:vector from camera position to look at position
	camera_Zaxis[1] = m_camera.lookat[Y] - m_camera.position[Y];
	camera_Zaxis[2] = m_camera.lookat[Z] - m_camera.position[Z];

	Normalization(camera_Zaxis);										// camera Z-axis Z = cl/||cl||

	float scalar = DotProduct(m_camera.worldup, camera_Zaxis);			//up' = up - (up . Z)Z
	camera_Yaxis[0] = m_camera.worldup[X] - scalar * camera_Zaxis[0];
	camera_Yaxis[1] = m_camera.worldup[Y] - scalar * camera_Zaxis[1];
	camera_Yaxis[2] = m_camera.worldup[Z] - scalar * camera_Zaxis[2];

	Normalization(camera_Yaxis);										// camera Y-axis Y = up'/||up'||

	CrossProduct(camera_Yaxis, camera_Zaxis, camera_Xaxis);				// camera X-axis X = (Y x Z)

	m_camera.Xiw[0][0] = camera_Xaxis[0]; m_camera.Xiw[0][1] = camera_Xaxis[1]; m_camera.Xiw[0][2] = camera_Xaxis[2]; m_camera.Xiw[0][3] = -DotProduct(camera_Xaxis, m_camera.position);
	m_camera.Xiw[1][0] = camera_Yaxis[0]; m_camera.Xiw[1][1] = camera_Yaxis[1]; m_camera.Xiw[1][2] = camera_Yaxis[2]; m_camera.Xiw[1][3] = -DotProduct(camera_Yaxis, m_camera.position);
	m_camera.Xiw[2][0] = camera_Zaxis[0]; m_camera.Xiw[2][1] = camera_Zaxis[1]; m_camera.Xiw[2][2] = camera_Zaxis[2]; m_camera.Xiw[2][3] = -DotProduct(camera_Zaxis, m_camera.position);
	m_camera.Xiw[3][0] = 0; m_camera.Xiw[3][1] = 0; m_camera.Xiw[3][2] = 0; m_camera.Xiw[3][3] = 1;
	// store Xiw into m_camera

	//Xpi
	m_camera.Xpi[0][0] = 1; m_camera.Xpi[0][1] = 0; m_camera.Xpi[0][2] = 0; m_camera.Xpi[0][3] = 0;
	m_camera.Xpi[1][0] = 0; m_camera.Xpi[1][1] = 1; m_camera.Xpi[1][2] = 0; m_camera.Xpi[1][3] = 0;
	m_camera.Xpi[2][0] = 0; m_camera.Xpi[2][1] = 0; m_camera.Xpi[2][2] = tan((m_camera.FOV / 2)*PI / 180.0); m_camera.Xpi[2][3] = 0;
	m_camera.Xpi[3][0] = 0; m_camera.Xpi[3][1] = 0; m_camera.Xpi[3][2] = tan((m_camera.FOV / 2)*PI / 180.0); m_camera.Xpi[3][3] = 1;
	// store Xpi into m_camera

	//init Ximage
	GzPushMatrix(Xsp);													//push Xsp, Ximage has [Xsp] in stack
	GzPushMatrix(m_camera.Xpi);											//push Xpi, Ximage has [Xsp Xpi] in stack
	GzPushMatrix(m_camera.Xiw);											//push Xiw, Ximage has [Xsp Xpi Xiw] in stack

	return GZ_SUCCESS;
}

int GzRender::GzPutCamera(GzCamera camera)
{
/* HW 3.8 
/*- overwrite renderer camera structure with new camera definition
*/
	m_camera = camera;

	return GZ_SUCCESS;
}

int GzRender::GzPushMatrix(GzMatrix	matrix)
{
/* HW 3.9 
- push a matrix onto the Ximage stack
- check for stack overflow
*/
	if (matlevel > MATLEVELS || matlevel < 0) {												// stack overflow
		return GZ_FAILURE;
	}
	if (matlevel == 0) {																	// when Ximage empty put the matrix in it
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				Ximage[matlevel][i][j] = matrix[i][j];
			}
		}
	}
	else {																					// when  Ximage is not empty, matrix pushed is multiplication of last matrix on stack and new matrix 
		GzMatrix* matrix_pointer = MatrixMultiplication(Ximage[matlevel - 1], matrix);		// call multiplication function
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				Ximage[matlevel][i][j] = (*matrix_pointer)[i][j];
			}
		}
	}

	GzMatrix IdentityMatrix = { {1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1} };

	if (matlevel == 0 || matlevel == 1) {													// when the matrix pushed is Xsp or Xpi					
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				Xnorm[matlevel][i][j] = IdentityMatrix[i][j];								// push identity matrix to Xnorm
			}
		}
	}
	else {						
		UnitaryMatrix(matrix);																// unitary matrix
		GzMatrix* matrix_pointer = MatrixMultiplication(Xnorm[matlevel - 1], matrix);		// call multiplication function
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				Xnorm[matlevel][i][j] = (*matrix_pointer)[i][j];							// push new matrix to the top stack of Xnorm
			}
		}
	}
	matlevel++;																				// Ximage pointer increament

	return GZ_SUCCESS;
}

int GzRender::GzPopMatrix()
{
/* HW 3.10
- pop a matrix off the Ximage stack
- check for stack underflow
*/
	if (matlevel <= 0)																		// stack underflow
		return GZ_FAILURE;
	else {
		matlevel--;																			// stack  pointer decrement
	}

	return GZ_SUCCESS;
}

int GzRender::GzPut(int AAindex, int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z)
{
/* HW1.4 write pixel values into the buffer */
	if (r < 0)										//clamp  r within  0-4095  boundary
		r = 0;
	else if (r > 4095)
		r = 4095;

	if (g < 0)										//clamp  g within  0-4095  boundary
		g = 0;
	else if (g > 4095)
		g = 4095;

	if (b < 0)										//clamp  b within  0-4095  boundary
		b = 0;
	else if (b > 4095)
		b = 4095;

	if (i >= 0 && i < xres && j >= 0 && j < yres ) {					//check i,j value is in boudary
		if (z < AApixelbuffers[AAindex][ARRAY(i, j)].z) {				//check z  value , keep the pixel with lower z 
			AApixelbuffers[AAindex][ARRAY(i, j)] = { r, g, b, a, z };	//write proccessed pixel value in to  pixelbuffer	
		}																
	}

	return GZ_SUCCESS;
}


int GzRender::GzGet(int i, int j, GzIntensity *r, GzIntensity *g, GzIntensity *b, GzIntensity *a, GzDepth *z)
{
/* HW1.5 retrieve a pixel information from the pixel buffer */
	if (i >= 0 && i < xres && j >= 0 && j < yres) {										//check i,j value is in boudary
		GzPixel pixel = pixelbuffer[ARRAY(i, j)];										//retrieve  the pixel at locatio  (i,j)
		*r = pixel.red;
		*g = pixel.green;
		*b = pixel.blue;
		*a = pixel.alpha;
		*z = pixel.z;
	} 

	return GZ_SUCCESS;
}


int GzRender::GzFlushDisplay2File(FILE* outfile)
{
/* HW1.6 write image to ppm file -- "P6 %d %d 255\r" */
	GzIntensity r, g, b, a;																//create temp storage of r,g,b,a,z
	GzDepth z;

	fprintf(outfile, "P6 %d %d 255\r", xres, yres);										//write down the header for PPM file

	for (int j = 0; j < yres; j++) {													//go through row start from 0  to yres
		for (int i = 0; i < xres; i++) {												//go through each element on row i
			GzGet(i, j, &r, &g, &b, &a, &z);											//retrieve piixel information
			fprintf(outfile, "%c%c%c", (char)(r >> 4), (char)(g >> 4), (char)(b >> 4));	//print out RGB value to .ppm file
																						//drop  4 bits  from 16 bit signed short(12 bits valid value) to 8 bits 
																						//then cast unsigned short to unsinged char 
		}
	}

	return GZ_SUCCESS;
}

int GzRender::GzFlushDisplay2FrameBuffer()
{
/* HW1.7 write pixels to framebuffer: 
	- put the pixels into the frame buffer
	- CAUTION: when storing the pixels into the frame buffer, the order is blue, green, and red 
	- NOT red, green, and blue !!!
*/
	GzIntensity r, g, b, a;												//create temp storage of r,g,b,a,z
	GzDepth z;
	int count = 0;

	for (int j = 0; j < yres; j++) {									//go through row start from 0  to yres
		for (int i = 0; i < xres; i++) {								//go through each element on row i
			GzGet(i, j, &r, &g, &b, &a, &z);							//retrieve piixel information
			framebuffer[count++] = (char)(b >> 4);						//store b,g,r value into  framebuffer
			framebuffer[count++] = (char)(g >> 4);						//drop  4 bits  from 16 bit signed short(12 bits valid value) to 8 bits 
			framebuffer[count++] = (char)(r >> 4);						//then cast unsigned short to unsinged char 
		}
	}

	return GZ_SUCCESS;
}


/***********************************************/
/* HW2 methods: implement from here */

int GzRender::GzPutAttribute(int numAttributes, GzToken	*nameList, GzPointer *valueList) 
{
/* HW 2.1
-- Set renderer attribute states (e.g.: GZ_RGB_COLOR default color)
-- In later homeworks set shaders, interpolaters, texture maps, and lights
*/

/*
- GzPutAttribute() must accept the following tokens/values:

- GZ_RGB_COLOR					GzColor		default flat-shade color
- GZ_INTERPOLATE				int			shader interpolation mode
- GZ_DIRECTIONAL_LIGHT			GzLight
- GZ_AMBIENT_LIGHT            	GzLight		(ignore direction)
- GZ_AMBIENT_COEFFICIENT		GzColor		Ka reflectance
- GZ_DIFFUSE_COEFFICIENT		GzColor		Kd reflectance
- GZ_SPECULAR_COEFFICIENT       GzColor		Ks coef's
- GZ_DISTRIBUTION_COEFFICIENT   float		spec power
*/
	for (int attributeNum = 0; attributeNum < numAttributes; attributeNum++) {		//go through each attribute in  Tooken and value
		switch (nameList[attributeNum]) {											//find name of attriibute in name list

		case GZ_RGB_COLOR: 															//this attribute is GZ_RGB_Color	
		{
			GzColor* gzColorPointer = (GzColor*)(valueList[attributeNum]);			//conver a  GzPointer into GZColorPointer
			flatcolor[0] = ctoi((*gzColorPointer)[0]);								//conver R value from float to GzIntensity
			flatcolor[1] = ctoi((*gzColorPointer)[1]);								//conver G value from float to GzIntensity
			flatcolor[2] = ctoi((*gzColorPointer)[2]);								//conver B value from float to GzIntensity
		}
		break;																		//then pass RGB value to Render

		case GZ_INTERPOLATE:														//this attribute is GZ_INTERPOLATE
		{
			int* modePointer = (int*)(valueList[attributeNum]);						//retrieve interpolate mode from valueList
			interp_mode = *modePointer;												//store the interpolate mode
		}
		break;

		case GZ_DIRECTIONAL_LIGHT: 													//this attribute is GZ_DIRECTIONAL_LIGHT	
		{
			GzLight* gzLightPointer = (GzLight*)(valueList[attributeNum]);			//retrieve light direction and color from valuelList
			lights[attributeNum] = *gzLightPointer;									//store light infor into lights list
			numlights++;															//increament light count

			// Added a light, calculate pyramid with light position as apex and some base defined by a plane
			
			// Calculating base of pyramid
			// Base of pyramid has normal equal to the light vector pointing from origin towards light
			GzCoord planeEquation;
			SetCoordEqual(planeEquation, (*gzLightPointer).direction);

			// To define this plane we need two basis vectors of the plane
			GzCoord firstBasis;
			GzCoord secondBasis;
			GzCoord nonDirVector;

			
			SetCoordEqual(firstBasis, ZeroCoord);
			SetCoordEqual(secondBasis, ZeroCoord);


			// First take some random vector and cross it with the normal vector of the plane

			// nonDirVector is some random vector that I've tried to make sure will always be 
			// perpendicular to the vector "planeEquation"
			SetCoordEqual(nonDirVector, planeEquation);
			MultiplyVector3(nonDirVector, -1.0f);
			for (int i = 0; i < 3; i++)
			{
				nonDirVector[i] += (float)rand();
			}

			// first Basis vector is random vector cross normal
			// second basis vector is normal cross first basis vector
			CrossVec3(planeEquation, nonDirVector, firstBasis);
			CrossVec3(planeEquation, firstBasis, secondBasis);

			// Normalize these two basis vectors
			NormalizeVector3(firstBasis);
			NormalizeVector3(secondBasis);

			float BasisLength = (float)sqrt(pow(max_plane_height, 2) + pow(max_plane_width, 2));

			MultiplyVector3(firstBasis, BasisLength);
			MultiplyVector3(secondBasis, BasisLength);

			float CamToPlane[4];

			for (int i = 0; i < 4; i++)
			{
				SetCoordEqual(max_plane_pos[i], ZeroCoord);
			}

			// Now we get four points that form the base of the light pyramid

			AddVector3(max_plane_pos[0], firstBasis);
			AddVector3(max_plane_pos[1], secondBasis);

			MultiplyVector3(firstBasis, -1.0f);
			MultiplyVector3(secondBasis, -1.0f);

			AddVector3(max_plane_pos[2], firstBasis);
			AddVector3(max_plane_pos[3], secondBasis);

			//move light frustm plane to z axis for 20.0f
			for (int i = 0; i < 4; i++)
			{
				max_plane_pos[i][1] -= 10.0f;
				max_plane_pos[i][2] += 50.0f;
			}

			// We get the vector from camera position to origin

			GzCoord VectorToCamera;
			SetCoordEqual(VectorToCamera, m_camera.position);
			MultiplyVector3(VectorToCamera, -1.0f);

			// Calculate which of the points on the base are closest and furthest
			// from the camera

			float max_dist = -FLT_MAX;
			int max_i = 0;
			float min_dist = FLT_MAX;
			int min_i = 0;
			for (int i = 0; i < 4; i++)
			{

				GzCoord Dist;
				// Set vector to the vector from origin to base point
				SetCoordEqual(Dist, max_plane_pos[i]);
				// Add this vector to vector from camera to origin
				AddVector3(Dist, VectorToCamera);
				// This now gives us the vector from camera to base point
				// Calculate its magnitude and see if it is a minimum or maximum
				float newDist = GetMagnitudeVector3(Dist);
				if (newDist > max_dist)
				{
					max_dist = newDist;
					max_i = i;
				}
				if (newDist < min_dist && newDist > 0)
				{
					min_dist = newDist;
					min_i = i;
				}
			}

			// Now take the vector from the camera to the min and max and project onto camera forward
			GzCoord CameraForward;
			CameraForward[0] = m_camera.lookat[X] - m_camera.position[X];
			CameraForward[1] = m_camera.lookat[Y] - m_camera.position[Y];
			CameraForward[2] = m_camera.lookat[Z] - m_camera.position[Z];

			NormalizeVector3(CameraForward);

			GzCoord EyeToMinPoint, EyeToMaxPoint;
			SetCoordEqual(EyeToMinPoint, max_plane_pos[min_i]);
			SetCoordEqual(EyeToMaxPoint, max_plane_pos[max_i]);

			SubtractVector3(EyeToMinPoint, m_camera.position);
			SubtractVector3(EyeToMaxPoint, m_camera.position);

			//float EyeToMinDot = DotVec3(EyeToMinPoint, CameraForward);	
			//float EyeToMaxDot = DotVec3(EyeToMaxPoint, CameraForward);
			Projection(EyeToMinPoint, CameraForward);
			Projection(EyeToMaxPoint, CameraForward);
			float EyeToMin = GetMagnitudeVector3(EyeToMinPoint);
			float EyeToMax = GetMagnitudeVector3(EyeToMaxPoint);


			//SetCoordEqual(EyeToMinPoint, CameraForward);
			//SetCoordEqual(EyeToMaxPoint, CameraForward);

			// Now create sampling planes between EyeToMinPoint, EyeToMaxPoint
			NumSamplingPlanes = (int)ceil((EyeToMax - EyeToMin) / delta_t) + 1;
			samplingPlanes = new GZSAMPLINGPLANE[NumSamplingPlanes];
			for (int t = 0; t < NumSamplingPlanes; t++)
			{
				// Create a sampling plane
				samplingPlanes[t].samplingPlanePixels = new GzCoord*[SamplingPlaneY];
				for (int z = 0; z < SamplingPlaneY; z++)
				{
					samplingPlanes[t].samplingPlanePixels[z] = new GzCoord[SamplingPlaneX];
				}
				samplingPlanes[t].DistanceFromEye = EyeToMin + delta_t * t;		

				SetCoordEqual(samplingPlanes[t].midPointPosition, CameraForward);				
				MultiplyVector3(samplingPlanes[t].midPointPosition, samplingPlanes[t].DistanceFromEye);
				AddVector3(samplingPlanes[t].midPointPosition, m_camera.position);

				SetCoordEqual(samplingPlanes[t].PlaneNormal, CameraForward);
				MultiplyVector3(samplingPlanes[t].PlaneNormal, -1.0f);
				NormalizeVector3(samplingPlanes[t].PlaneNormal);

				samplingPlanes[t].samplePlaneBuffer = new GzPixel[xres * yres];

			}

			// Created sampling planes

			// How to verify? Perhaps render a triangle on each sampling plane

		}
		break;

		case GZ_AMBIENT_LIGHT: 														//this attribute is GZ_AMBIENT_LIGHT	
		{
			GzLight* gzLightPointer = (GzLight*)(valueList[attributeNum]);			//retrieve light direction and color from valuelList
			ambientlight = *gzLightPointer;											//store light as ambientlight
		}
		break;

		case GZ_AMBIENT_COEFFICIENT:												//this attribute is GZ_AMBIENT_COEFFICIENT
		{
			GzColor* gzColorPointer = (GzColor*)(valueList[attributeNum]);			//retrieve color pointer from valueList
			Ka[0] = (*gzColorPointer)[0];											//store color as Ka
			Ka[1] = (*gzColorPointer)[1];
			Ka[2] = (*gzColorPointer)[2];
		}
		break;

		case GZ_DIFFUSE_COEFFICIENT:												//this attribute is GZ_DIFFUSE_COEFFICIENT	
		{
			GzColor* gzColorPointer = (GzColor*)(valueList[attributeNum]);			//retrieve color pointer from valueList
			Kd[0] = (*gzColorPointer)[0];											//store color as Kd
			Kd[1] = (*gzColorPointer)[1];
			Kd[2] = (*gzColorPointer)[2];
		}
		break;


		case GZ_SPECULAR_COEFFICIENT:												//this attribute is GZ_SPECULAR_COEFFICIENT
		{
			GzColor* gzColorPointer = (GzColor*)(valueList[attributeNum]);			//retrieve color pointer from valueList
			Ks[0] = (*gzColorPointer)[0];											//store color as Ks
			Ks[1] = (*gzColorPointer)[1];
			Ks[2] = (*gzColorPointer)[2];
		}
		break;

		case GZ_DISTRIBUTION_COEFFICIENT:											//this attribute is GZ_DISTRIBUTION_COEFFICIENT
		{
			spec = *((float*)(valueList[attributeNum]));							//store spec value from valueList
		}
		break;

		case  GZ_TEXTURE_MAP:
		{
			tex_fun = *(GzTexture)(valueList[attributeNum]);
		}
		break;

		}
	}

	return GZ_SUCCESS;
}

int GzRender::GzDebugRenderSamplingPlanes()
{
	GzToken		nameListTriangle[3]; 	/* vertex attribute names */
	GzPointer	valueListTriangle[3]; 	/* vertex attribute pointers */
	GzCoord		vertexList[3];	/* vertex position coordinates */
	GzCoord		normalList[3];	/* vertex normals */
	GzTextureIndex  	uvList[3];		/* vertex texture map indices */

	nameListTriangle[0] = GZ_POSITION;
	nameListTriangle[1] = GZ_NORMAL;
	nameListTriangle[2] = GZ_TEXTURE_INDEX;

	GzCoord PlaneNormal;
	SetCoordEqual(PlaneNormal, samplingPlanes[0].PlaneNormal);

	// To define this plane we need two basis vectors of the plane
	//GzCoord firstBasis;
	//GzCoord secondBasis;
	//GzCoord nonDirVector;

	//SetCoordEqual(firstBasis, ZeroCoord);
	//SetCoordEqual(secondBasis, ZeroCoord);

	//// First take some random vector and cross it with the normal vector of the plane

	//// nonDirVector is some random vector that I've tried to make sure will always be 
	//// perpendicular to the vector "planeEquation"
	//SetCoordEqual(nonDirVector, PlaneNormal);
	//MultiplyVector3(nonDirVector, -1.0f);
	//for (int i = 0; i < 3; i++)
	//{
	//	nonDirVector[i] += (float)rand();
	//}

	//// first Basis vector is random vector cross normal
	//// second basis vector is normal cross first basis vector
	//CrossVec3(PlaneNormal, nonDirVector, firstBasis);
	//CrossVec3(PlaneNormal, firstBasis, secondBasis);

	//NormalizeVector3(firstBasis);
	//NormalizeVector3(secondBasis);

	//float BasisLength = (float)sqrt(pow(max_plane_height, 2) + pow(max_plane_width, 2));

	//MultiplyVector3(firstBasis, BasisLength);
	//MultiplyVector3(secondBasis, BasisLength);

	for (int t = 0; t < NumSamplingPlanes; t++)
	{

		GzCoord PlaneOrigin;
		SetCoordEqual(PlaneOrigin, samplingPlanes[t].midPointPosition);
		
		GzCoord TV1, TV2;
		SetCoordEqual(TV1, PlaneOrigin);
		SetCoordEqual(TV2, PlaneOrigin);

		// To define this plane we need two basis vectors of the plane
		GzCoord firstBasis;
		GzCoord secondBasis;
		GzCoord nonDirVector;

		SetCoordEqual(firstBasis, ZeroCoord);
		SetCoordEqual(secondBasis, ZeroCoord);

		// First take some random vector and cross it with the normal vector of the plane

		// nonDirVector is some random vector that I've tried to make sure will always be 
		// perpendicular to the vector "planeEquation"
		SetCoordEqual(nonDirVector, PlaneNormal);
		MultiplyVector3(nonDirVector, -1.0f);

		float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

		for (int i = 0; i < 3; i++)
		{
			nonDirVector[i] += (float)rand();
			nonDirVector[i] *= r;
			nonDirVector[i] -= (float)rand();
		}

		// first Basis vector is random vector cross normal
		// second basis vector is normal cross first basis vector
		CrossVec3(PlaneNormal, nonDirVector, firstBasis);
		CrossVec3(PlaneNormal, firstBasis, secondBasis);

		NormalizeVector3(firstBasis);
		NormalizeVector3(secondBasis);

		float BasisLength = (float)sqrt(pow(max_plane_height, 2) + pow(max_plane_width, 2));

		MultiplyVector3(firstBasis, BasisLength);
		MultiplyVector3(secondBasis, BasisLength);

		AddVector3(TV1, firstBasis);
		AddVector3(TV2, secondBasis);

		
		MultiplyVector3(TV1, r * 2.0f);
		MultiplyVector3(TV2, r * 2.0f);

		// Time to render triangle
		SetCoordEqual(vertexList[0], PlaneOrigin);
		SetCoordEqual(vertexList[1], TV1);
		SetCoordEqual(vertexList[2], TV2);

		SetCoordEqual(normalList[0], PlaneNormal);
		SetCoordEqual(normalList[1], PlaneNormal);
		SetCoordEqual(normalList[2], PlaneNormal);

		uvList[0][0] = 0.0f;
		uvList[0][1] = 0.0f;
		uvList[1][0] = 1.0f;
		uvList[1][1] = 1.0f;
		uvList[2][0] = 0.0f;
		uvList[2][1] = 1.0f;

		//vertexList[0][1] -= 0.25f;

		valueListTriangle[0] = (GzPointer)vertexList;
		valueListTriangle[1] = (GzPointer)normalList;
		valueListTriangle[2] = (GzPointer)uvList;
		GzPutTriangle(3, nameListTriangle, valueListTriangle);
	}



	return 0;
}

int GzRender::GzCalculateSamplingPlanes()
{
	FILE *MyFile;
	if ((MyFile = fopen("test.txt", "wb")) == NULL)
	{
		AfxMessageBox("The output file was not opened\n");
		return GZ_FAILURE;
	}


	// At this point, you have added a directional light source
	// You also have constructed the sampling planes 

	// Now calculate rendering at each point
	
	// Need to create UV coordinate system as needed for Equation 9

	// Origin of UV coordinate system is position of light source
	// light source position is defined as some distance in the direction
	// of light direction away from origin
	GzCoord lightPosition;
	SetCoordEqual(lightPosition, lights[0].direction);
	MultiplyVector3(lightPosition, -lightSourceOffset); // lightSourceOffset defined in rend.h, 10.0f
	lightPosition[2] = 20.0f;							//lightScource exact on middle of plane



	// Find Point Q from Figure 2 in UV system
	// To find Point Q take vector from camera to light position, and project it onto CameraForward
	
	GzCoord CameraToLight;
	SetCoordEqual(CameraToLight, lightPosition);
	SubtractVector3(CameraToLight, m_camera.position);

	GzCoord CameraForward;
	CameraForward[0] = m_camera.lookat[X] - m_camera.position[X];
	CameraForward[1] = m_camera.lookat[Y] - m_camera.position[Y];
	CameraForward[2] = m_camera.lookat[Z] - m_camera.position[Z];

	NormalizeVector3(CameraForward);

	GzCoord CameraLeftVector;
	SetCoordEqual(CameraLeftVector, ZeroCoord);
	CrossVec3(CameraForward, m_camera.worldup, CameraLeftVector);
	NormalizeVector3(CameraLeftVector);

	GzCoord CameraUpVector;
	SetCoordEqual(CameraUpVector, m_camera.worldup);
	NormalizeVector3(CameraUpVector);



	// Now going to calculate the textures of all the sampling planes
	for (int t = 0; t < NumSamplingPlanes; t++)
	{
		// Looking at equation 9, this is the ck constant, aka exp(-Bptk)

		float ck = exp(-extinctionCoefficient * atmosphericDensity * samplingPlanes[t].DistanceFromEye);
		float frustumheight = 2.0f * samplingPlanes[t].DistanceFromEye * tan(m_camera.FOV * 0.5f * PI / 180);


///////////////////////////////////////////////////
		//set  light frustum range is that plane

		GzCoord BaseLeft, BaseRight, BaseFar, BaseClose;
		SetCoordEqual(BaseLeft, max_plane_pos[0]);
		SetCoordEqual(BaseClose, max_plane_pos[1]);
		SetCoordEqual(BaseRight, max_plane_pos[2]);
		SetCoordEqual(BaseFar, max_plane_pos[3]);

		GzCoord BaseLeftPoint, BaseRightPoint, TopPoint;

		if (samplingPlanes[t].DistanceFromEye <= lightPosition[Z]) {		//between close palne and middle point

			BaseLeftPoint[X] = BaseLeft[X] / (BaseLeft[Z] - BaseClose[Z]) * (samplingPlanes[t].midPointPosition[Z] - BaseClose[Z]);
			BaseLeftPoint[Y] = BaseLeft[Y];
			BaseLeftPoint[Z] = samplingPlanes[t].midPointPosition[Z];

			BaseRightPoint[X] = BaseRight[X] / (BaseRight[Z] - BaseClose[Z]) * (samplingPlanes[t].midPointPosition[Z] - BaseClose[Z]);
			BaseRightPoint[Y] = BaseRight[Y];
			BaseRightPoint[Z] = samplingPlanes[t].midPointPosition[Z];

			TopPoint[X] = lightPosition[X];
			TopPoint[Y] = lightPosition[Z] / lightPosition[Y] * samplingPlanes[t].midPointPosition[Z];
			TopPoint[Z] = samplingPlanes[t].midPointPosition[Z];

		}

		else {																//middle point and far plane

			BaseLeftPoint[X] = BaseLeft[X] / (BaseLeft[Z] - BaseClose[Z]) * (BaseFar[Z] - samplingPlanes[t].midPointPosition[Z]);
			BaseLeftPoint[Y] = BaseLeft[Y];
			BaseLeftPoint[Z] = samplingPlanes[t].midPointPosition[Z];

			BaseRightPoint[X] = BaseRight[X] / (BaseRight[Z] - BaseClose[Z]) * (BaseFar[Z] - samplingPlanes[t].midPointPosition[Z]);
			BaseRightPoint[Y] = BaseRight[Y];
			BaseRightPoint[Z] = samplingPlanes[t].midPointPosition[Z];

			TopPoint[X] = lightPosition[X];
			TopPoint[Y] = lightPosition[Z] / lightPosition[Y] * samplingPlanes[t].midPointPosition[Z];
			TopPoint[Z] = samplingPlanes[t].midPointPosition[Z];

		}

////////////////////////////////////////////////////////
		// CONFUSION, if we have UV plane it is only 2D. However the sampling plane is in 3D. I don't
		// know how the sampling plane is constructed. This code will construct the sampling plane
		// along the midpoint in the direction of the V axis. It will reuse the same value at each 
		// vertical V offset in the entire horizontal row of the sampling plane. 
		for (int j = 0; j < SamplingPlaneX; j++)
		{

			// Calculate horizontal offset from midpoint position
			GzCoord CameraToSamplingPoint;
			GzCoord SamplingPointPos;
			SetCoordEqual(SamplingPointPos, CameraLeftVector);
			float vertOffsetX = -(j - SamplingPlaneX / 2) * frustumheight/ SamplingPlaneX;
			MultiplyVector3(SamplingPointPos, -vertOffsetX);
			AddVector3(SamplingPointPos, samplingPlanes[t].midPointPosition);
			
			SetCoordEqual(CameraToSamplingPoint, SamplingPointPos);
			SubtractVector3(CameraToSamplingPoint, m_camera.position);
			NormalizeVector3(CameraToSamplingPoint);

			// This is how far Q is from camera position along CameraToSamplingPoint
			double QpointDist = DotVec3(CameraToLight, CameraToSamplingPoint);

			// Actual position of Q
			GzCoord QPoint;
			SetCoordEqual(QPoint, CameraToSamplingPoint);
			MultiplyVector3(QPoint, QpointDist);
			AddVector3(QPoint, m_camera.position);

			// The U axis of the UV plane is simply the camera forward
			// The V axis is defined as a vector from light source to the Q point

			GzCoord LightToQ;
			SetCoordEqual(LightToQ, QPoint);
			SubtractVector3(LightToQ, lightPosition);


			GzCoord QToLight;
			SetCoordEqual(QToLight, LightToQ);
			MultiplyVector3(QToLight, -1.0f);
			NormalizeVector3(QToLight);

			GzCoord QToCamera;
			SetCoordEqual(QToCamera, m_camera.position);
			SubtractVector3(QToCamera, QPoint);
			NormalizeVector3(QToCamera);
		
			for (int i = 0; i < SamplingPlaneY; i++)
			{

				// Get the midpoint of the sampling plane
				GzCoord Midpoint;
				SetCoordEqual(Midpoint, SamplingPointPos);

				// Apply a vertical offset to the sampling plane along the V axis (aka QToLight)
				GzCoord MidpointVertOffset;
				SetCoordEqual(MidpointVertOffset, CameraUpVector);
				float vertOffsetY = -(i - SamplingPlaneY / 2) * frustumheight / SamplingPlaneY;
				MultiplyVector3(MidpointVertOffset, vertOffsetY);
				AddVector3(Midpoint, MidpointVertOffset);

				// The actual u value of this point is calculated by taking a vector from Q to this point
				// and projecting it onto camera forward (aka the u axis)
				GzCoord QToMidpoint;
				SetCoordEqual(QToMidpoint, Midpoint);
				SubtractVector3(QToMidpoint, QPoint);

				float u = GetMagnitudeVector3(QToMidpoint);

				if (DotVec3(QToMidpoint, QToCamera) > 0)
				{
					u *= -1.0f;
				}

				// The v value of this point is simply the distance from LightToQ
				// ASSUMPTION, LIGHT IS ABOVE THE CAMERA SO THIS VALUE WILL ALWAYS BE POSITIVE
				float v = GetMagnitudeVector3(LightToQ);

				// Now calculate the integral for the space between this sampling plane (inclusive) and the next
				// sampling plane (exclusive)

				// See Equation 9 for calculation of q
				float q = 0.0f;

////////////////////////////////////////////////////
								//test for equation 3
				GzCoord pixelPoint, LightToPoint, EyeToPoint;
				pixelPoint[0] = vertOffsetX; pixelPoint[1] = vertOffsetY; pixelPoint[2] = samplingPlanes[0].midPointPosition[2];

				SetCoordEqual(LightToPoint, pixelPoint);
				SubtractVector3(LightToPoint, lightPosition);
				float sValue = GetMagnitudeVector3(LightToPoint);

				SetCoordEqual(EyeToPoint, pixelPoint);
				SubtractVector3(EyeToPoint, m_camera.position);
				float tValue = GetMagnitudeVector3(EyeToPoint);


				float gt = atmosphericDensity * i/1000 * exp(-extinctionCoefficient * atmosphericDensity * (sValue + tValue)) / (sValue*sValue);

				GzCoord forAngle;

				float angle = DotVec3(LightToPoint, EyeToPoint) / (GetMagnitudeVector3(LightToPoint) *GetMagnitudeVector3(EyeToPoint));
				float attenu = (GetMagnitudeVector3(LightToPoint) - vertOffsetY) / GetMagnitudeVector3(LightToPoint);
				
//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
				//get H(t) 1 visible 0 unvisible
				bool visible = 0;

				float edgeEqu1 = (BaseRightPoint[1] - TopPoint[1])*(pixelPoint[0] - TopPoint[0]) - (BaseRightPoint[0] - TopPoint[0])*(pixelPoint[1] - TopPoint[1]);
				float edgeEqu2 = (BaseLeftPoint[1] - BaseRightPoint[1])*(pixelPoint[0] - BaseRightPoint[0]) - (BaseLeftPoint[0] - BaseRightPoint[0])*(pixelPoint[1] - BaseRightPoint[1]);
				float edgeEqu3 = (TopPoint[1] - BaseLeftPoint[1])*(pixelPoint[0] - BaseLeftPoint[0]) - (TopPoint[0] - BaseLeftPoint[0])*(pixelPoint[1] - BaseLeftPoint[1]);
				if (edgeEqu1 <= 0 && edgeEqu2 <= 0 && edgeEqu3 <= 0 || edgeEqu1 >= 0 && edgeEqu2 >= 0 && edgeEqu3 >= 0) visible = 1;

//////////////////////////////////////////////////////

				// This function in Integrate.h will be integrated from the sampling point's t value to the next t value
				LightFunctor func(u, v, extinctionCoefficient, atmosphericDensity, gt);
				// q = qtrap(func, u, u + delta_t);
				Trapzd<LightFunctor> s(func, u, u + delta_t);
				for (int j = 1; j <= INTEGRATION_STEPS_POW + 1; j++)
				{
					q = s.next();
				}

				// Multiply the value of q 
				q *= ck;

				samplingPlanes[t].samplingPlanePixels[i][j][0] = q * visible;
				samplingPlanes[t].samplingPlanePixels[i][j][1] = q * visible;
				samplingPlanes[t].samplingPlanePixels[i][j][2] = q * visible;
			}
		}
	}

	float accumulatedRed = 0.0f;
	float accumulatedGreen = 0.0f;
	float accumulatedBlue = 0.0f;

	float MaxValue = -FLT_MAX;
	float MinValue = FLT_MAX; 
	
	//FILE *MyFile;
	//if ((MyFile = fopen("test.txt", "wb")) == NULL)
	//{
	//	AfxMessageBox("The output file was not opened\n");
	//	return GZ_FAILURE;
	//}

	//fprintf(MyFile, "vertext1 %f%f%f\n", max_plane_pos[0][0], max_plane_pos[0][1], max_plane_pos[0][2]);
	//fprintf(MyFile, "vertext2 %f%f%f\n", max_plane_pos[1][0], max_plane_pos[1][1], max_plane_pos[1][2]);
	//fprintf(MyFile, "vertext3 %f%f%f\n", max_plane_pos[2][0], max_plane_pos[2][1], max_plane_pos[2][2]);
	//fprintf(MyFile, "vertext4 %f%f%f\n", max_plane_pos[3][0], max_plane_pos[3][1], max_plane_pos[3][2]);


	for (int j = 0; j < SamplingPlaneX; j++)
	{
		for (int i = 0; i < SamplingPlaneY ; i++)
		{
			for (int t = 0; t < NumSamplingPlanes; t++)
			{
				accumulatedRed += samplingPlanes[t].samplingPlanePixels[i][j][0];
				accumulatedGreen += samplingPlanes[t].samplingPlanePixels[i][j][1];
				accumulatedBlue += samplingPlanes[t].samplingPlanePixels[i][j][2];

			}

			if (accumulatedBlue > MaxValue) MaxValue = accumulatedBlue;
			if (accumulatedBlue < MinValue) MinValue = accumulatedBlue;

			/*double f;
			int e;

			f = frexp(accumulatedBlue, &e);
			f *= 4095;*/

			//pixelbuffer[i * SamplingPlaneX + j].red += (short)Clamp(f, 0, 4095);
			//pixelbuffer[i * SamplingPlaneX + j].green += (short)Clamp(f, 0, 4095);
			//pixelbuffer[i * SamplingPlaneX + j].blue += (short)Clamp(f, 0, 4095);
			//pixelbuffer[i * SamplingPlaneX + j].blue = (short)Clamp((float)pixelbuffer[i * SamplingPlaneX + j].blue, 0, 4095);

			accumulatedRed = 0.0f;
			accumulatedGreen = 0.0f;
			accumulatedBlue = 0.0f;

		}
	}

	for (int j = 0; j < SamplingPlaneX; j++)
	{
		for (int i = 0; i < SamplingPlaneY; i++)
		{
			for (int t = 0; t < NumSamplingPlanes; t++)
			{
				accumulatedRed += samplingPlanes[t].samplingPlanePixels[i][j][0];
				accumulatedGreen += samplingPlanes[t].samplingPlanePixels[i][j][1];
				accumulatedBlue += samplingPlanes[t].samplingPlanePixels[i][j][2];

			}

			/*double f;
			int e;

			f = frexp(accumulatedBlue, &e);
			f *= 4095;*/

			float f = (accumulatedBlue - MinValue) / (MaxValue - MinValue);
			f *= 4095;

			pixelbuffer[i * SamplingPlaneX + j].red += (short)Clamp(f, 0, 4095);
			pixelbuffer[i * SamplingPlaneX + j].green += (short)Clamp(f, 0, 4095);
			pixelbuffer[i * SamplingPlaneX + j].blue += (short)Clamp(f, 0, 4095);

			pixelbuffer[i * SamplingPlaneX + j].red = (short)Clamp((float)pixelbuffer[i * SamplingPlaneX + j].red, 0, 4095);
			pixelbuffer[i * SamplingPlaneX + j].green = (short)Clamp((float)pixelbuffer[i * SamplingPlaneX + j].green, 0, 4095);
			pixelbuffer[i * SamplingPlaneX + j].blue = (short)Clamp((float)pixelbuffer[i * SamplingPlaneX + j].blue, 0, 4095);

			accumulatedRed = 0.0f;
			accumulatedGreen = 0.0f;
			accumulatedBlue = 0.0f;

		}

	}
	



	// Now we have sampling planes... how to render?
	// Going to try creating two triangles and texturing them
	// See GzDebugRenderSamplingPlanes, where the UVs are now being set to some values
	// And see tex_samplingPlane function


	return 0;
}


int GzRender::tex_samplingPlane(float u, float v, GzColor color)
{
	unsigned char		pixel[3];
	unsigned char     dummy;
	char  		foo[8];
	int   		i, j;
	FILE			*fd;

	/* bounds-test u,v to make sure nothing will overflow image array bounds */

	u = u < 0 ? 0 : u;
	u = u > 1 ? 1 : u;
	v = v < 0 ? 0 : v;
	v = v > 1 ? 1 : v;

	/* determine texture cell corner values and perform bilinear interpolation */
	i = floor(u * (SamplingPlaneX - 1));
	j = floor(v * (SamplingPlaneY - 1));

	/*for (int k = 0; k < NumSamplingPlanes; k++)
	{
		color[0] += samplingPlanes[k].samplingPlanePixels[i][j].red;
		color[1] += samplingPlanes[k].samplingPlanePixels[i][j].green;
		color[2] += samplingPlanes[k].samplingPlanePixels[i][j].blue;
	}*/

	return GZ_SUCCESS;
}

int GzRender::TestSamplingPlaneOutput() {
	for (int j = 0; j < yres; j++) {
		for (int i = 0; i < xres; i++) {
			//pixelbuffer[ARRAY(i, j)].red = samplingPlanes[0].samplePlaneBuffer[j*yres + i].red;
			//pixelbuffer[ARRAY(i, j)].green = samplingPlanes[0].samplePlaneBuffer[j*yres + i].green;
			//pixelbuffer[ARRAY(i, j)].blue = samplingPlanes[0].samplePlaneBuffer[j*yres + i].blue;
			pixelbuffer[ARRAY(i, j)].alpha = 1;
			pixelbuffer[ARRAY(i, j)].z = 1;
		}
	}

	return 0;
}

int GzRender::TestForSamplePlanePassingInfor() {
	
	GzCoord SamplePlaneOneTopLeftPos;
	float frustumheight = 2.0f * samplingPlanes[0].DistanceFromEye * tan(m_camera.FOV * 0.5f * PI / 180);
	SamplePlaneOneTopLeftPos[0] = -frustumheight / 2;
	SamplePlaneOneTopLeftPos[1] = frustumheight / 2;
	SamplePlaneOneTopLeftPos[2] = 0;

	for (int j = 0; j < yres; j++) {								
		for (int i = 0; i < xres; i++) {
			samplingPlanes[0].samplePlaneBuffer[j*yres + i].red = 4095;
			samplingPlanes[0].samplePlaneBuffer[j*yres + i].green = 0;
			samplingPlanes[0].samplePlaneBuffer[j*yres + i].blue = 0;
			samplingPlanes[0].samplePlaneBuffer[j*yres + i].alpha =1;
			samplingPlanes[0].samplePlaneBuffer[j*yres + i].z = 1;
		}
	}

	for (int j = 0; j < yres; j++) {
		for (int i = 0; i < xres/2; i++) {
			samplingPlanes[1].samplePlaneBuffer[j*yres + i].red = 0;
			samplingPlanes[1].samplePlaneBuffer[j*yres + i].green = 4095;
			samplingPlanes[1].samplePlaneBuffer[j*yres + i].blue = 0;
			samplingPlanes[1].samplePlaneBuffer[j*yres + i].alpha = 1;
			samplingPlanes[1].samplePlaneBuffer[j*yres + i].z = 1;
		}
	}

	for (int j = 0; j < yres/2; j++) {
		for (int i = 0; i < xres; i++) {
			samplingPlanes[2].samplePlaneBuffer[j*yres + i].red = 0;
			samplingPlanes[2].samplePlaneBuffer[j*yres + i].green = 0;
			samplingPlanes[2].samplePlaneBuffer[j*yres + i].blue = 4095;
			samplingPlanes[2].samplePlaneBuffer[j*yres + i].alpha = 1;
			samplingPlanes[2].samplePlaneBuffer[j*yres + i].z = 1;
		}
	}

	for (int j = 0; j < yres; j++) {
		for (int i = 0; i < xres; i++) {
			pixelbuffer[ARRAY(i, j)].red = samplingPlanes[0].samplePlaneBuffer[j*yres + i].red + samplingPlanes[1].samplePlaneBuffer[j*yres + i].red + samplingPlanes[2].samplePlaneBuffer[j*yres + i].red;
			pixelbuffer[ARRAY(i, j)].green = samplingPlanes[0].samplePlaneBuffer[j*yres + i].green+samplingPlanes[1].samplePlaneBuffer[j*yres + i].green+ samplingPlanes[2].samplePlaneBuffer[j*yres + i].green;
			pixelbuffer[ARRAY(i, j)].blue = samplingPlanes[0].samplePlaneBuffer[j*yres + i].blue+ samplingPlanes[1].samplePlaneBuffer[j*yres + i].blue+ samplingPlanes[2].samplePlaneBuffer[j*yres + i].blue;
			pixelbuffer[ARRAY(i, j)].alpha = 1;
			pixelbuffer[ARRAY(i, j)].z = 1;
		}
	}
	
	return GZ_SUCCESS;

}

int GzRender::GzPutTriangle(int numParts, GzToken *nameList, GzPointer *valueList)
/* numParts - how many names and values */
{
/* HW 2.2
-- Pass in a triangle description with tokens and values corresponding to
      GZ_NULL_TOKEN:		do nothing - no values
      GZ_POSITION:		3 vert positions in model space
-- Return error code
*/
/*
-- Xform positions of verts using matrix on top of stack 
-- Clip - just discard any triangle with any vert(s) behind view plane 
		- optional: test for triangles with all three verts off-screen (trivial frustum cull)
-- invoke triangle rasterizer  
*/
	GzCoord	Vertex[3][3];																		//Vertex[vertex number 1/2/3][vertex infor coord/normal/texture][coord in 3D space x/y/z]
	GzColor defaultColor = { 1,1,1 };
	int blueColor = 0;

	for (int attributeNum = 0; attributeNum < numParts; attributeNum++) {							//go through each attribute in nameList and valueList

		switch (nameList[attributeNum]) {															//find the name in nameList with corresponding attribute

		case GZ_NULL_TOKEN:																			//if nameList has GZ_NULL_TOKEN
			break;																					//return with no value											

		case GZ_POSITION:																			//if nameList has GZ_POSITION
		{
			GzCoord* VertCoordPointer = (GzCoord*)(valueList[attributeNum]);						//get a GzCoorPointer with three vertex coordinate in it

			for (int vertexNum = 0; vertexNum < 3; vertexNum++) {
				for (int coord = 0; coord < 3; coord++) {
					Vertex[vertexNum][0][coord] = VertCoordPointer[vertexNum][coord];				//store three vertex coordinate in world space
				}
				if (MatrixVectorMultiplication(Ximage[matlevel - 1], Vertex[vertexNum][0],true) < 0)//transform vertex coord to screen space
					GZ_FAILURE;
			}

		}
		break;

		case GZ_NORMAL:																				//nameList is GZ_NORMAL, calculate normal
		{
			GzCoord* NormCoordPointer = (GzCoord*)(valueList[attributeNum]);						//get a GzCoorPointer with three normal coordinate in it
			
			for (int vertexNum = 0; vertexNum < 3; vertexNum++) {
				for (int coord = 0; coord < 3; coord++) {
					Vertex[vertexNum][1][coord] = NormCoordPointer[vertexNum][coord];				//store three vertex normal in world space
				}
				MatrixVectorMultiplication(Xnorm[matlevel - 1], Vertex[vertexNum][1],false);		//transform norm1 to image space	
			}
		}
		break;

		case GZ_TEXTURE_INDEX:																		//NameList is TextureIndex
		{
			GzTextureIndex* TextureIndexPointer = (GzTextureIndex*)(valueList[attributeNum]);		//get TextureIndexPointer with uv in it
			for (int vertexNum = 0; vertexNum < 3; vertexNum++) {
				for (int coord = 0; coord < 2; coord++) {
					Vertex[vertexNum][2][coord] = TextureIndexPointer[vertexNum][coord];			//store uv coordinate
				}
			}
		}
		break;

		case GZ_RGB_COLOR:
		{
			int* useBlue = (int*)(valueList[attributeNum]);
			if (*useBlue == 1)
			{
				blueColor = 1;
			}
		}
		break;

		}
	}

	GzIntensity r, g, b, a;																			//create temp storage of r,g,b,a,z
	GzDepth z;
	GzCoord	InputVertex[3][3];

	float AAFilter[AAKERNEL_SIZE][3] 		/* X-shift, Y-shift, weight */
	{
	-0.52, 0.38, 0.128, 		0.41, 0.56, 0.119,		0.27, 0.08, 0.294,
	-0.17, -0.29, 0.249,		0.58, -0.55, 0.104,		-0.31, -0.71, 0.106
	};

	memcpy(InputVertex, Vertex, sizeof(Vertex));													//make a copy of Vertex

	for (int AA = 0; AA < AAKERNEL_SIZE; AA++) {

		memcpy(Vertex, InputVertex, sizeof(Vertex));												//retrieve the Vertex infor

		for (int vertexNum = 0; vertexNum < 3; vertexNum++) {
			for (int axis = 0; axis < 2; axis++) {
				Vertex[vertexNum][0][axis] -= AAFilter[AA][axis];									//offset vertex 
			}
		}

		// 1.Sort verts by Y															
		Sort(Vertex);

		// Ground Shading: get three vertex color
		GzColor* colorPointer;
		
		if (tex_fun != 0) {
			Ks[0] = Kd[0] = Ka[0] = defaultColor[0]; Ks[0] = Kd[1] = Ka[1] = defaultColor[1]; Ks[0] = Kd[2] = Ka[2] = defaultColor[2];
		}

		colorPointer = LightingEq(Ks, Kd, Ka, lights, numlights, ambientlight, Vertex[0][1], spec);	//get vertex color when set Ks,Ka,Kd to 1
		GzColor vert1_color = { (*colorPointer)[0],(*colorPointer)[1],(*colorPointer)[2] };
		colorPointer = LightingEq(Ks, Kd, Ka, lights, numlights, ambientlight, Vertex[1][1], spec);
		GzColor vert2_color = { (*colorPointer)[0],(*colorPointer)[1],(*colorPointer)[2] };
		colorPointer = LightingEq(Ks, Kd, Ka, lights, numlights, ambientlight, Vertex[2][1], spec);
		GzColor vert3_color = { (*colorPointer)[0],(*colorPointer)[1],(*colorPointer)[2] };

		// 2.Setup  edge  DDAs for (1-2),(2-3), (1-3)
		Edge vert1_vert2 = SetUpEdge((Vertex[0][0]), (Vertex[1][0]));
		Edge vert2_vert3 = SetUpEdge((Vertex[1][0]), (Vertex[2][0]));
		Edge vert1_vert3 = SetUpEdge((Vertex[0][0]), (Vertex[2][0]));

		// 3.Sort edge by L or R
		Edge leftEdge;
		Edge rightEdge;
		bool longEdgeLeft;
		if (vert1_vert2.slope_x < vert1_vert3.slope_x) {										//L:1,2,3 R:1,3
			leftEdge = vert1_vert2;
			rightEdge = vert1_vert3;
			longEdgeLeft = false;
		}
		else {																					//L:1,3 R:1,2,3
			leftEdge = vert1_vert3;
			rightEdge = vert1_vert2;
			longEdgeLeft = true;
		}

		//5. Advance leftEdge and rightEdge DDA current position to top y-scan line (ceiling)
		float DeltaY = ceil(vert1_vert3.current[1]) - vert1_vert3.current[1];					//find deltaY as difference between Y and ceil(Y)
		UpdateDDA(&leftEdge, DeltaY);															//advance (1-2) to top y-scan line
		UpdateDDA(&rightEdge, DeltaY);															//advance (1-3) to top y-scan line
		int edgeStart = vert1_vert3.current[1];

		//12. keep running until current y > Y[3]
		while (leftEdge.current[1] < leftEdge.end[1] || rightEdge.current[1] < rightEdge.end[1]) {//end loop when both edge's current value are out of boundary	
			//11.Switch from 1-2 edge to 2-3 edge when current 1-2 pisition > Y[2]
			if (!longEdgeLeft && leftEdge.current[1] > leftEdge.end[1]) {						//when (1-3) on right and (1-2)finish rasterization
				leftEdge = vert2_vert3;															//update leftEdge as (2-3)
				float DeltaY = rightEdge.current[1] - leftEdge.current[1];						//get difference between current left and right edge
				UpdateDDA(&leftEdge, DeltaY);													//update left DDA to current y-scan line
			}
			if (longEdgeLeft && rightEdge.current[1] > rightEdge.end[1]) {						//when (1-3) on left and (1-2)finish rasterization
				rightEdge = vert2_vert3;														//update rightEdge as (2-3)
				float DeltaY = leftEdge.current[1] - rightEdge.current[1];						//get difference between current left and right edge
				UpdateDDA(&rightEdge, DeltaY);													//update right DDA to current y-scan line
			}

			//6.Setup span DDA on successuve lines based on edge DDA
			//7.Setup span's value
			Span  span = SetUpSpan(leftEdge.current, rightEdge.current);						//setup span between (1-2) and (1-3) current position
			//8.Adavance span  surrent position to left-most covered pixel
			float DeltaX = ceil(span.current[0]) - span.current[0];								//find  deltaX as difference between X and ceil(X)
			UpdateSpan(&span, DeltaX);															//advance span to left most
			//9.Interpolate span position and parameter (Z) until current position > end
			while (span.current[0] < span.end[0]) {												//go through each pixel point on span

				//Shading: deal with color here
				GzColor color = { 0,0,0 };
				if (blueColor == 1)
				{
					vert1_color[0] = 0.0f;
					vert1_color[1] = 0.0f;
					vert1_color[2] = 1.0f;

					vert2_color[0] = 0.0f;
					vert2_color[1] = 0.0f;
					vert2_color[2] = 1.0f;

					vert3_color[0] = 0.0f;
					vert3_color[1] = 0.0f;
					vert3_color[2] = 1.0f;
				}
				//texture Color
				float u, v;
				GzColor textureColor = { 1,1,1 };

				if (tex_fun != 0) {
					PInterpolation(&u, &v, Vertex, span.current[0], span.current[1]);
					// tex_fun(u, v, textureColor);
					tex_samplingPlane(u, v, textureColor);
					Kd[0] = Ka[0] = textureColor[0]; Kd[1] = Ka[1] = textureColor[1]; Kd[2] = Ka[2] = textureColor[2];
				}

				if (interp_mode == GZ_COLOR) {
					//ground shading

					// Shading: get Plane Equation for color interpolation	
					ColorInterpolation(color, Vertex, vert1_color, vert2_color, vert3_color, span.current[0], span.current[1]);		//interpolate color

					color[0] *= textureColor[0];																					//adjust color with texture color
					color[1] *= textureColor[1];
					color[2] *= textureColor[2];
					
				}

				else if (interp_mode == GZ_NORMALS) {
					//phong shading				
					GzCoord normal = { 0,0,0 };
					NormalInterpolation(normal, Vertex, span.current[0], span.current[1]);											//interpolate normal

					colorPointer = LightingEq(Ks, Kd, Ka, lights, numlights, ambientlight, normal, spec);		//get color at certain piexel with texture
					color[0] = (*colorPointer)[0]; color[1] = (*colorPointer)[1]; color[2] = (*colorPointer)[2];
				}
				r = ctoi((color)[0]);								//conver R value from float to GzIntensity
				g = ctoi((color)[1]);								//conver G value from float to GzIntensity
				b = ctoi((color)[2]);								//conver B value from float to GzIntensity

				//10. Test interpolated-Z against FB-Z for each pixel and update into buffer
				GzPut(AA, span.current[0], span.current[1], r, g, b, 1, span.current[2]);						//write  color value to FB

				UpdateSpan(&span, 1);															//move span's current position to next right pixel
			}

			UpdateDDA(&leftEdge, 1);															//move left DDA edge down one pixel
			UpdateDDA(&rightEdge, 1);															//move right DDA edge down one pixel
		}

	}

	for (int j = 0; j < yres; j++) {									//go through row start from 0  to yres
		for (int i = 0; i < xres; i++) {								//go through each element on row i
			pixelbuffer[ARRAY(i, j)].red = 0;
			pixelbuffer[ARRAY(i, j)].green = 0;
			pixelbuffer[ARRAY(i, j)].blue = 0;
			pixelbuffer[ARRAY(i, j)].alpha = 0;
			pixelbuffer[ARRAY(i, j)].z = 0;

			for (int AA = 0; AA < AAKERNEL_SIZE; AA++) {
				pixelbuffer[ARRAY(i, j)].red += AApixelbuffers[AA][ARRAY(i, j)].red * AAFilter[AA][2];
				pixelbuffer[ARRAY(i, j)].green += AApixelbuffers[AA][ARRAY(i, j)].green * AAFilter[AA][2];
				pixelbuffer[ARRAY(i, j)].blue += AApixelbuffers[AA][ARRAY(i, j)].blue * AAFilter[AA][2];
				pixelbuffer[ARRAY(i, j)].alpha += AApixelbuffers[AA][ARRAY(i, j)].alpha * AAFilter[AA][2];
				pixelbuffer[ARRAY(i, j)].z += AApixelbuffers[AA][ARRAY(i, j)].z * AAFilter[AA][2];						//write  color value to FB
			}
		}
	}

	//adding  filter as light beam

	//for (int j = 0; j < yres; j++) {									//go through row start from 0  to yres

	//	int midPoint = 120;
	//	int length = ceil ((20.0 + (80.0 / yres) * j)/2);
	//	for (int i = midPoint - length; i < midPoint + length; i++) {								//go through each element on row i
	//		
	//		float distance = sqrt((i - midPoint)*(i - midPoint) + j * j);
	//		float attenuation = (250 - distance) / 250; 

	//		float lightColorR = Clamp(pixelbuffer[ARRAY(i, j)].red * 5.0, 0, 4095);
	//		float lightColorG = Clamp(pixelbuffer[ARRAY(i, j)].green * 3.8, 0, 4095);
	//		float lightColorB = Clamp(pixelbuffer[ARRAY(i, j)].blue * 2.2, 0, 4095);

	//		float originalColorR = Clamp(pixelbuffer[ARRAY(i, j)].red, 0, 4095);
	//		float originalColorG = Clamp(pixelbuffer[ARRAY(i, j)].green, 0, 4095);
	//		float originalColorB = Clamp(pixelbuffer[ARRAY(i, j)].blue, 0, 4095);

	//		float attenuatedRed = Clamp(lightColorR * attenuation, originalColorR, lightColorR);
	//		float attenuatedGreen = Clamp(lightColorG * attenuation, originalColorG, lightColorG);
	//		float attenuatedBlue = Clamp(lightColorB * attenuation, originalColorB, lightColorB);

	//		pixelbuffer[ARRAY(i, j)].red = attenuatedRed;
	//		pixelbuffer[ARRAY(i, j)].green = attenuatedGreen;
	//		pixelbuffer[ARRAY(i, j)].blue = attenuatedBlue;
	//	
	//	}
	//}

	return GZ_SUCCESS;
}

