// from https://github.com/mnhrdt/imscript
#ifndef _BICUBIC_C
#define _BICUBIC_C

#include <math.h>

// extrapolate by nearest value
inline
static float getsample_1(float *x, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (l < 0) l = 0;
	if (i >= w) i = w-1;
	if (j >= h) j = h-1;
	if (l >= pd) l = pd-1;
	return x[(i+j*w)*pd + l];
}


static float cubic_interpolation(float v[4], float x)
{
	return v[1] + 0.5 * x*(v[2] - v[0]
			+ x*(2.0*v[0] - 5.0*v[1] + 4.0*v[2] - v[3]
			+ x*(3.0*(v[1] - v[2]) + v[3] - v[0])));
}

static float bicubic_interpolation_cell(float p[4][4], float x, float y)
{
	float v[4];
	v[0] = cubic_interpolation(p[0], y);
	v[1] = cubic_interpolation(p[1], y);
	v[2] = cubic_interpolation(p[2], y);
	v[3] = cubic_interpolation(p[3], y);
	return cubic_interpolation(v, x);
}

void bicubic_interpolation(float *result, float *img, int w, int h, int pd, float x, float y)
{
	x -= 1;
	y -= 1;

	int ix = floor(x);
	int iy = floor(y);
	for (int l = 0; l < pd; l++) {
		float c[4][4];
		for (int j = 0; j < 4; j++)
			for (int i = 0; i < 4; i++)
				c[i][j] = getsample_1(img, w, h, pd, ix + i, iy + j, l);
		float r = bicubic_interpolation_cell(c, x - ix, y - iy);
		result[l] = r;
	}
}

#endif//_BICUBIC_C

