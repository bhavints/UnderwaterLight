/* CS580 Homework 3 */

#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"
#include	"dda.h"
#include	"matrixOperation.h"
#include	"shading.h"

#define PI (float) 3.14159265358979323846

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
	delete AApixelbuffers;
}

int GzRender::GzDefault()
{
/* HW1.3 set pixel buffer to some default values - start a new frame */
		for (int j = 0; j < yres; j++) {									//go through row start from 0  to yres
			for (int i = 0; i < xres; i++) {								//go through each element on row 	
				for (int AA = 0; AA < AAKERNEL_SIZE; AA++) {				//go through each AA pixelbuffers
					AApixelbuffers[AA][ARRAY(i, j)].red = 4095;				//get a color I like, a = 1, z = max
					AApixelbuffers[AA][ARRAY(i, j)].green = 2794;
					AApixelbuffers[AA][ARRAY(i, j)].blue = 3227;
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
		GzColor defaultColor = { 1,1,1 };
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

				//texture Color
				float u, v;
				GzColor textureColor = { 1,1,1 };

				if (tex_fun != 0) {
					PInterpolation(&u, &v, Vertex, span.current[0], span.current[1]);
					tex_fun(u, v, textureColor);
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

	return GZ_SUCCESS;
}

