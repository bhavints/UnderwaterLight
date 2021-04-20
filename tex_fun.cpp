/* Texture functions for cs580 GzLib	*/
#include    "stdafx.h" 
#include	"stdio.h"
#include	"Gz.h"

GzColor	*image=NULL;
GzColor	*imageCat, *imageDog;
GzColor	*imageWater;
int xs, ys;
int reset = 1;

/* Image texture function */
int tex_fun(float u, float v, GzColor color)
{
  unsigned char		pixel[3];
  unsigned char     dummy;
  char  		foo[8];
  int   		i, j;
  FILE			*fd;

  if (reset) {          /* open and load texture file */
    fd = fopen ("cat.ppm", "rb");
    if (fd == NULL) {
      fprintf (stderr, "texture file not found\n");
      exit(-1);
    }
    fscanf (fd, "%s %d %d %c", foo, &xs, &ys, &dummy);
    image = (GzColor*)malloc(sizeof(GzColor)*(xs+1)*(ys+1));
    if (image == NULL) {
      fprintf (stderr, "malloc for texture image failed\n");
      exit(-1);
    }

    for (i = 0; i < xs*ys; i++) {	/* create array of GzColor values */
      fread(pixel, sizeof(pixel), 1, fd);
      image[i][RED] = (float)((int)pixel[RED]) * (1.0 / 255.0);
      image[i][GREEN] = (float)((int)pixel[GREEN]) * (1.0 / 255.0);
      image[i][BLUE] = (float)((int)pixel[BLUE]) * (1.0 / 255.0);
      }

    reset = 0;          /* init is done */
	fclose(fd);
  }

/* bounds-test u,v to make sure nothing will overflow image array bounds */

  u = u < 0 ? 0 : u;
  u = u > 1 ? 1 : u;
  v = v < 0 ? 0 : v;
  v = v > 1 ? 1 : v;

/* determine texture cell corner values and perform bilinear interpolation */
  i = floor(u * (xs - 1));
  j = floor(v * (ys - 1));

  GzColor A_Color = { image[j * xs + i][RED],image[j  *xs + i][GREEN],image[j * xs + i][BLUE] };
  GzColor B_Color = { image[j * xs + i + 1][RED],image[j * xs + i + 1][GREEN],image[j * xs + i + 1][BLUE] };
  GzColor C_Color = { image[(j + 1) * xs + i + 1][RED],image[(j + 1) * xs + i + 1][GREEN],image[(j + 1) * xs + i + 1][BLUE] };
  GzColor D_Color = { image[(j + 1) * xs + i][RED],image[(j + 1) * xs + i][GREEN],image[(j + 1) * xs + i][BLUE] };
  
  float s = u * (xs - 1) - i;
  float t = v * (ys - 1) - j;

  /* set color to interpolated GzColor value and return */
  for (int k = 0; k < 3; k++) {
	  color[k] = s * t * C_Color[k] + (1 - s) * t * D_Color[k] + s * (1 - t) * B_Color[k] + (1 - s) * (1 - t) * A_Color[k];
  }

	return GZ_SUCCESS;
}

/* Procedural texture function */
int ptex_water(float u, float v, GzColor color)
{
	// checkerboard
	float N = 6;
	int u_interval = ceil(u * 6);
	int v_interval = ceil(v * 6);

	u = u * 6 - floor(u * 6);					//offset u,v value
	v = v * 6 - floor(v * 6);


	unsigned char		pixelc[3], pixeld[3];
	unsigned char     dummy;
	char  		foo[8];
	int   		i, j;
	FILE		*fdc;

	if (reset) {          /* open and load texture file */
		fdc = fopen("water_transparent.ppm", "rb");
		if (fdc == NULL) {
			fprintf(stderr, "texture file not found\n");
			exit(-1);
		}
		fscanf(fdc, "%s %d %d %c", foo, &xs, &ys, &dummy);
		imageWater = (GzColor*)malloc(sizeof(GzColor)*(xs + 1)*(ys + 1));
		if (imageWater == NULL) {
			fprintf(stderr, "malloc for texture image failed\n");
			exit(-1);
		}

		for (i = 0; i < xs*ys; i++) {	/* create array of GzColor values */
			fread(pixelc, sizeof(pixelc), 1, fdc);
			imageWater[i][RED] = (float)((int)pixelc[RED]) * (1.0 / 255.0);
			imageWater[i][GREEN] = (float)((int)pixelc[GREEN]) * (1.0 / 255.0);
			imageWater[i][BLUE] = (float)((int)pixelc[BLUE]) * (1.0 / 255.0);
		}

		reset = 0;          /* init is done */
		fclose(fdc);
	}

	/* bounds-test u,v to make sure nothing will overflow image array bounds */

	u = u < 0 ? 0 : u;
	u = u > 1 ? 1 : u;
	v = v < 0 ? 0 : v;
	v = v > 1 ? 1 : v;

	/* determine texture cell corner values and perform bilinear interpolation */
	i = floor(u * (xs - 1));
	j = floor(v * (ys - 1));

	image = imageWater;

	GzColor A_Color = { image[j * xs + i][RED],image[j  *xs + i][GREEN],image[j * xs + i][BLUE] };
	GzColor B_Color = { image[j * xs + i + 1][RED],image[j * xs + i + 1][GREEN],image[j * xs + i + 1][BLUE] };
	GzColor C_Color = { image[(j + 1) * xs + i + 1][RED],image[(j + 1) * xs + i + 1][GREEN],image[(j + 1) * xs + i + 1][BLUE] };
	GzColor D_Color = { image[(j + 1) * xs + i][RED],image[(j + 1) * xs + i][GREEN],image[(j + 1) * xs + i][BLUE] };

	float s = u * (xs - 1) - i;
	float t = v * (ys - 1) - j;

	/* set color to interpolated GzColor value and return */
	for (int k = 0; k < 3; k++) {
		color[k] = s * t * C_Color[k] + (1 - s) * t * D_Color[k] + s * (1 - t) * B_Color[k] + (1 - s) * (1 - t) * A_Color[k];
	}

	return GZ_SUCCESS;
}

/* Procedural texture function */
int ptex_fun(float u, float v, GzColor color)
{
	// checkerboard
	float N = 6;
	int u_interval = ceil(u * 6);
	int v_interval = ceil(v * 6);
	int imageType = 0;							//dog image

	if (u_interval % 2 == v_interval % 2) {
		imageType = 1;							// cat image
	}
	
	u = u * 6 - floor(u * 6);					//offset u,v value
	v = v * 6 - floor(v * 6);


	unsigned char		pixelc[3], pixeld[3];
	unsigned char     dummy;
	char  		foo[8];
	int   		i, j;
	FILE		*fdc, *fdd;

	if (reset) {          /* open and load texture file */
		fdc = fopen("cat.ppm", "rb");
		fdd = fopen("dog.ppm", "rb");
		if (fdc == NULL || fdd == NULL) {
			fprintf(stderr, "texture file not found\n");
			exit(-1);
		}
		fscanf(fdc, "%s %d %d %c", foo, &xs, &ys, &dummy);
		fscanf(fdd, "%s %d %d %c", foo, &xs, &ys, &dummy);
		imageCat = (GzColor*)malloc(sizeof(GzColor)*(xs + 1)*(ys + 1));
		imageDog = (GzColor*)malloc(sizeof(GzColor)*(xs + 1)*(ys + 1));
		if (imageCat == NULL || imageDog == NULL) {
			fprintf(stderr, "malloc for texture image failed\n");
			exit(-1);
		}

		for (i = 0; i < xs*ys; i++) {	/* create array of GzColor values */
			fread(pixelc, sizeof(pixelc), 1, fdc);
			imageCat[i][RED] = (float)((int)pixelc[RED]) * (1.0 / 255.0);
			imageCat[i][GREEN] = (float)((int)pixelc[GREEN]) * (1.0 / 255.0);
			imageCat[i][BLUE] = (float)((int)pixelc[BLUE]) * (1.0 / 255.0);

			fread(pixeld, sizeof(pixeld), 1, fdd);
			imageDog[i][RED] = (float)((int)pixeld[RED]) * (1.0 / 255.0);
			imageDog[i][GREEN] = (float)((int)pixeld[GREEN]) * (1.0 / 255.0);
			imageDog[i][BLUE] = (float)((int)pixeld[BLUE]) * (1.0 / 255.0);
		}

		reset = 0;          /* init is done */
		fclose(fdc); 
		fclose(fdd);
	}

	/* bounds-test u,v to make sure nothing will overflow image array bounds */

	u = u < 0 ? 0 : u;
	u = u > 1 ? 1 : u;
	v = v < 0 ? 0 : v;
	v = v > 1 ? 1 : v;

	/* determine texture cell corner values and perform bilinear interpolation */
	i = floor(u * (xs - 1));
	j = floor(v * (ys - 1));

	if (imageType == 0)
		image = imageCat;
	else
		image = imageDog;

	GzColor A_Color = { image[j * xs + i][RED],image[j  *xs + i][GREEN],image[j * xs + i][BLUE] };
	GzColor B_Color = { image[j * xs + i + 1][RED],image[j * xs + i + 1][GREEN],image[j * xs + i + 1][BLUE] };
	GzColor C_Color = { image[(j + 1) * xs + i + 1][RED],image[(j + 1) * xs + i + 1][GREEN],image[(j + 1) * xs + i + 1][BLUE] };
	GzColor D_Color = { image[(j + 1) * xs + i][RED],image[(j + 1) * xs + i][GREEN],image[(j + 1) * xs + i][BLUE] };

	float s = u * (xs - 1) - i;
	float t = v * (ys - 1) - j;

	/* set color to interpolated GzColor value and return */
	for (int k = 0; k < 3; k++) {
		color[k] = s * t * C_Color[k] + (1 - s) * t * D_Color[k] + s * (1 - t) * B_Color[k] + (1 - s) * (1 - t) * A_Color[k];
	}

	return GZ_SUCCESS;
}

/* Free texture memory */
int GzFreeTexture()
{
	if(image!=NULL)
		free(image);
	return GZ_SUCCESS;
}