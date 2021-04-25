#include	"gz.h"
#ifndef GZRENDER_
#define GZRENDER_


/* Camera defaults */
#define	DEFAULT_FOV		35.0
#define	DEFAULT_IM_Z	(-10.0)  /* world coords for image plane origin (-10,5,-10) */
#define	DEFAULT_IM_Y	(0.0)    /* default look-at point = 0,0,0 */
#define	DEFAULT_IM_X	(0.0)

#define	DEFAULT_AMBIENT	{0.1, 0.1, 0.1}
#define	DEFAULT_DIFFUSE	{0.7, 0.6, 0.5}
#define	DEFAULT_SPECULAR	{0.2, 0.3, 0.4}
#define	DEFAULT_SPEC		32

#define	MATLEVELS	100		/* how many matrix pushes allowed */
#define	MAX_LIGHTS	10		/* how many lights allowed */

#define INTEGRATION_STEPS_POW 5

class GzRender{			/* define a renderer */
  

public:
	unsigned short	xres;
	unsigned short	yres;
	GzPixel		*pixelbuffer;		/* frame buffer array */
	GzPixel		**AApixelbuffers;
	char* framebuffer;
	int x, y, z;

	GzCamera		m_camera;
	short		    matlevel;	        /* top of stack - current xform */
	GzMatrix		Ximage[MATLEVELS];	/* stack of xforms (Xsm) */
	GzMatrix		Xnorm[MATLEVELS];	/* xforms for norms (Xim) */
	GzMatrix		Xsp;		        /* NDC to screen (pers-to-screen) */
	GzColor		flatcolor;          /* color state for flat shaded triangles */
	int			interp_mode;
	int			numlights;
	GzLight		lights[MAX_LIGHTS];
	GzLight		ambientlight;
	GzColor		Ka, Kd, Ks;
	float		    spec;		/* specular power */
	GzTexture		tex_fun;    /* tex_fun(float u, float v, GzColor color) */

	const int 	rWidth = 256;		// frame buffer and display width
	const int	rHeight = 256;

	// DEFINE VARIABLES FOR SAMPLING PLANES
	const float delta_t = 0.2f;; // Distance between sampling planes
	const float tmin = 0.5f;;
	const float tmax = 50.0f;
	const float d0 = 100.0; // Distance of the furthest sampling plane from the light source
	const float max_plane_width = 10.0f;
	const float max_plane_height = 10.0f;
	GZSAMPLINGPLANE		*samplingPlanes;		/* sampling plane buffer */
	int NumSamplingPlanes;
	const float lightSourceOffset = 10.0f;
	const float atmosphericDensity = 0.01f; // rho value
	const float extinctionCoefficient = 0.01f; // Beta value

	// Not sure what "units" sampling planes are in
	// Going to make up a random value between each sampling plane "pixel"

	const float samplingPixelDist = max_plane_width / (float)rWidth;

	const int SamplingPlaneX = 256;
	const int SamplingPlaneY = 256;
	
  	// Constructors
	GzRender(int xRes, int yRes);
	~GzRender();

	// Function declaration

	// HW1: Display methods
	int GzDefault();
	int GzBeginRender();
	int GzPut(int AA, int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z);
	int GzGet(int i, int j, GzIntensity *r, GzIntensity *g, GzIntensity *b, GzIntensity *a, GzDepth	*z);

	int GzFlushDisplay2File(FILE* outfile);
	int GzFlushDisplay2FrameBuffer();

	// HW2: Render methods
	int GzPutAttribute(int numAttributes, GzToken *nameList, GzPointer *valueList);
	int GzPutTriangle(int numParts, GzToken *nameList, GzPointer *valueList);

	// HW3
	int GzPutCamera(GzCamera camera);
	int GzPushMatrix(GzMatrix	matrix);
	int GzPopMatrix();

	// Extra methods: NOT part of API - just for general assistance */
	inline int ARRAY(int x, int y){return (x+y*xres);}	/* simplify fbuf indexing */
	inline short	ctoi(float color) {return(short)((int)(color * ((1 << 12) - 1)));}		/* convert float color to GzIntensity short */


	// Object Translation
	int GzRotXMat(float degree, GzMatrix mat);
	int GzRotYMat(float degree, GzMatrix mat);
	int GzRotZMat(float degree, GzMatrix mat);
	int GzTrxMat(GzCoord translate, GzMatrix mat);
	int GzScaleMat(GzCoord scale, GzMatrix mat);

	int GzDebugRenderSamplingPlanes();
	int GzCalculateSamplingPlanes();
	int TestForSamplePlanePassingInfor();
	int TestSamplingPlaneOutput();
	int tex_samplingPlane(float u, float v, GzColor color);
};
#endif