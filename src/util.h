///////////////////////////////////////////////////////////////////////// 
//
// util.h 
//
// copyright (c) 2016-2016, Renato Rodrigues Silva
// last modified Mar, 2016
// first written Jan, 2016
//
//     This program is free software; you can redistribute it and/or
//     modify it under the terms of the GNU General Public License,
//     version 3, as published by the Free Software Foundation.
//
//     This program is distributed in the hope that it will be useful,
//     but without any warranty; without even the implied warranty of
//     merchantability or fitness for a particular purpose.  See the GNU
//     General Public License, version 3, for more details.
//
//     A copy of the GNU General Public License, version 3, is available
//     at http://www.r-project.org/Licenses/GPL-3
//
// Part of the BayesPSE package
// Contains: addlog 
/////////////////////////////////////////////////////////////////////////// 


#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Applic.h>
#include <stdio.h>
#include <R_ext/Applic.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>


double addlog(double a, double b);

/* Ax = b */

void prodXy( double *X, int *dimX, double *y, double *ans);


