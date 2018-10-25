// from https://github.com/mnhrdt/imscript
#ifndef _EROSION_C
#define _EROSION_C

#include <stdlib.h>
#include <math.h>
#include <assert.h>

static float getpixel_nan(float *x, int w, int h, int i, int j)
{
	if (i < 0 || i >= w || j < 0 || j >= h)
		return NAN;
	return x[i + j*w];
}

void morsi_erosion(float *y, float *x, int w, int h, int *e)
{
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		float a = INFINITY;
		for (int k = 0; k < e[0]; k++)
			a = fminf(a, getpixel_nan(x,w,h, i-e[2]+e[2*k+4], j-e[3]+e[2*k+5]));
		y[j*w+i] = a;
	}
}

int *build_disk(float radius)
{
	if (!(radius > 1)) return NULL;
	int side = 2*radius+4, elen = 2*side*side+4;
	int *e = malloc(elen*sizeof*e), cx = 0;
	for (int i = -radius-1; i <= radius+1; i++)
	for (int j = -radius-1; j <= radius+1; j++)
		if (hypot(i,j) < radius) {
			e[2*cx+4] = i;
			e[2*cx+5] = j;
			cx += 1;
		}
	assert(cx < side*side);
	e[0] = cx;
	e[1] = e[2] = e[3] = 0;
	return e;
}

#endif// _EROSION_C

