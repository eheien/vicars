#include "Matrix.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

VCDenseStd::VCDenseStd(const unsigned int &ncols, const unsigned int &nrows) : VCDenseMatrix(ncols, nrows) {
	data = (double*)valloc(sizeof(double)*ncols*nrows);
	assert(data);
	for (unsigned int i=0;i<ncols*nrows;++i) data[i] = 0;
}

VCDenseStd::~VCDenseStd(void) {
	if (data) free(data);
	data = NULL;
}

unsigned long VCDenseStd::mem_bytes(void) const {
	return sizeof(double)*width*height;
}

double *VCDenseStdStraight::getRow(double *buf, const unsigned int &row) const {
	return &data[row*width];
}

double *VCDenseStdStraight::getCol(double *buf, const unsigned int &col) const {
	printf("not implemented");
	assert(false);
}
