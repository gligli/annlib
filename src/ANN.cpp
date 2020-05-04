//----------------------------------------------------------------------
// File:			ANN.cpp
// Programmer:		Sunil Arya and David Mount
// Description:		Methods for ANN.h and ANNx.h
// Last modified:	01/27/10 (Version 1.1.2)
//----------------------------------------------------------------------
// Copyright (c) 1997-2010 University of Maryland and Sunil Arya and
// David Mount.  All Rights Reserved.
// 
// This software and related documentation is part of the Approximate
// Nearest Neighbor Library (ANN).  This software is provided under
// the provisions of the Lesser GNU Public License (LGPL).  See the
// file ../ReadMe.txt for further information.
// 
// The University of Maryland (U.M.) and the authors make no
// representations about the suitability or fitness of this software for
// any purpose.  It is provided "as is" without express or implied
// warranty.
//----------------------------------------------------------------------
// History:
//	Revision 0.1  03/04/98
//		Initial release
//	Revision 1.0  04/01/05
//		Added performance counting to annDist()
//	Revision 1.1.2  01/27/10
//		Fixed minor compilation bugs for new versions of gcc
//----------------------------------------------------------------------

#include <cstdlib>						// C standard lib defs
#include <ANN/ANN.h>
#include <ANN/ANNx.h>					// all ANN includes
#include <ANN/ANNperf.h>				// ANN performance 

using namespace std;					// make std:: accessible

//----------------------------------------------------------------------
//	Point methods
//----------------------------------------------------------------------

//----------------------------------------------------------------------
//	Distance utility.
//		(Note: In the nearest neighbor search, most distances are
//		computed using partial distance calculations, not this
//		procedure.)
//----------------------------------------------------------------------

ANNdist annDist(						// interpoint squared distance
	int					dim,
	ANNpoint			p,
	ANNpoint			q)
{
	register int d;
	register ANNcoord diff;
	register ANNcoord dist;

	dist = 0;
	for (d = 0; d < dim; d++) {
		diff = p[d] - q[d];
		dist = ANN_SUM(dist, ANN_POW(diff));
	}
	ANN_FLOP(3*dim)					// performance counts
	ANN_PTS(1)
	ANN_COORD(dim)
	return dist;
}

//----------------------------------------------------------------------
//	annPrintPoint() prints a point to a given output stream.
//----------------------------------------------------------------------

void annPrintPt(						// print a point
	ANNpoint			pt,				// the point
	int					dim,			// the dimension
	std::ostream		&out)			// output stream
{
	for (int j = 0; j < dim; j++) {
		out << pt[j];
		if (j < dim-1) out << " ";
	}
}

//----------------------------------------------------------------------
//	Point allocation/deallocation:
//
//		Because points (somewhat like strings in C) are stored
//		as pointers.  Consequently, creating and destroying
//		copies of points may require storage allocation.  These
//		procedures do this.
//
//		annAllocPt() and annDeallocPt() allocate a deallocate
//		storage for a single point, and return a pointer to it.
//
//		annAllocPts() allocates an array of points as well a place
//		to store their coordinates, and initializes the points to
//		point to their respective coordinates.  It allocates point
//		storage in a contiguous block large enough to store all the
//		points.  It performs no initialization.
//
//		annDeallocPts() should only be used on point arrays allocated
//		by annAllocPts since it assumes that points are allocated in
//		a block.
//
//		annCopyPt() copies a point taking care to allocate storage
//		for the new point.
//
//		annAssignRect() assigns the coordinates of one rectangle to
//		another.  The two rectangles must have the same dimension
//		(and it is not possible to test this here).
//----------------------------------------------------------------------

ANNpoint annAllocPt(int dim, ANNcoord c)		// allocate 1 point
{
	ANNpoint p = new ANNcoord[dim];
	for (int i = 0; i < dim; i++) p[i] = c;
	return p;
}
   
ANNpointArray annAllocPts(int n, int dim)		// allocate n pts in dim
{
	ANNpointArray pa = new ANNpoint[n];			// allocate points
	ANNpoint	  p  = new ANNcoord[n*dim];		// allocate space for coords
	for (int i = 0; i < n; i++) {
		pa[i] = &(p[i*dim]);
	}
	return pa;
}

void annDeallocPt(ANNpoint &p)					// deallocate 1 point
{
	delete [] p;
	p = NULL;
}
   
void annDeallocPts(ANNpointArray &pa)			// deallocate points
{
	delete [] pa[0];							// dealloc coordinate storage
	delete [] pa;								// dealloc points
	pa = NULL;
}
   
ANNpoint annCopyPt(int dim, ANNpoint source)	// copy point
{
	ANNpoint p = new ANNcoord[dim];
	for (int i = 0; i < dim; i++) p[i] = source[i];
	return p;
}
   
												// assign one rect to another
void annAssignRect(int dim, ANNorthRect &dest, const ANNorthRect &source)
{
	for (int i = 0; i < dim; i++) {
		dest.lo[i] = source.lo[i];
		dest.hi[i] = source.hi[i];
	}
}

												// is point inside rectangle?
ANNbool ANNorthRect::inside(int dim, ANNpoint p)
{
	for (int i = 0; i < dim; i++) {
		if (p[i] < lo[i] || p[i] > hi[i]) return ANNfalse;
	}
	return ANNtrue;
}

//----------------------------------------------------------------------
//	Error handler
//----------------------------------------------------------------------

void annError(const char* msg, ANNerr level)
{
	if (level == ANNabort) {
		cerr << "ANN: ERROR------->" << msg << "<-------------ERROR\n";
		exit(1);
	}
	else {
		cerr << "ANN: WARNING----->" << msg << "<-------------WARNING\n";
	}
}

//----------------------------------------------------------------------
//	Gli: 'basic' procedural API
//----------------------------------------------------------------------

extern "C"
{
	DLL_API void * ann_kdtree_create(double ** pa, int n, int dd, int bs, ANNsplitRule split)
	{
		return new ANNkd_tree(pa, n, dd, bs, split);
	}

	DLL_API void ann_kdtree_destroy(void * akd)
	{
		delete (ANNkd_tree *)akd;
	}

	DLL_API int ann_kdtree_search(void * akd, double * q, double eps, double * err)
	{
		int idx = -1;
		double dd = 0.0;

		((ANNkd_tree *)akd)->annkSearch(q, 1, &idx, &dd, eps);

		if (err) *err = dd;
		return idx;
	}

	DLL_API int ann_kdtree_pri_search(void * akd, double * q, double eps, double * err)
	{
		int idx = -1;
		double dd = 0.0;

		((ANNkd_tree *)akd)->annkPriSearch(q, 1, &idx, &dd, eps);

		if (err) *err = dd;
		return idx;
	}

	DLL_API void ann_kdtree_search_multi(void * akd, int * idxs, double * errs, int cnt, double * q, double eps)
	{
		((ANNkd_tree *)akd)->annkSearch(q, cnt, idxs, errs, eps);
	}
}
