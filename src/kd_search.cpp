//----------------------------------------------------------------------
// File:			kd_search.cpp
// Programmer:		Sunil Arya and David Mount
// Description:		Standard kd-tree search
// Last modified:	01/04/05 (Version 1.0)
//----------------------------------------------------------------------
// Copyright (c) 1997-2005 University of Maryland and Sunil Arya and
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
//		Changed names LO, HI to ANN_LO, ANN_HI
//----------------------------------------------------------------------

#include "kd_search.h"					// kd-search declarations

#include "tmmintrin.h"


//----------------------------------------------------------------------
//	Approximate nearest neighbor searching by kd-tree search
//		The kd-tree is searched for an approximate nearest neighbor.
//		The point is returned through one of the arguments, and the
//		distance returned is the squared distance to this point.
//
//		The method used for searching the kd-tree is an approximate
//		adaptation of the search algorithm described by Friedman,
//		Bentley, and Finkel, ``An algorithm for finding best matches
//		in logarithmic expected time,'' ACM Transactions on Mathematical
//		Software, 3(3):209-226, 1977).
//
//		The algorithm operates recursively.  When first encountering a
//		node of the kd-tree we first visit the child which is closest to
//		the query point.  On return, we decide whether we want to visit
//		the other child.  If the box containing the other child exceeds
//		1/(1+eps) times the current best distance, then we skip it (since
//		any point found in this child cannot be closer to the query point
//		by more than this factor.)  Otherwise, we visit it recursively.
//		The distance between a box and the query point is computed exactly
//		(not approximated as is often done in kd-tree), using incremental
//		distance updates, as described by Arya and Mount in ``Algorithms
//		for fast vector quantization,'' Proc.  of DCC '93: Data Compression
//		Conference, eds. J. A. Storer and M. Cohn, IEEE Press, 1993,
//		381-390.
//
//		The main entry points is annkSearch() which sets things up and
//		then call the recursive routine ann_search().  This is a recursive
//		routine which performs the processing for one node in the kd-tree.
//		There are two versions of this virtual procedure, one for splitting
//		nodes and one for leaves.  When a splitting node is visited, we
//		determine which child to visit first (the closer one), and visit
//		the other child on return.  When a leaf is visited, we compute
//		the distances to the points in the buckets, and update information
//		on the closest points.
//
//		Some trickery is used to incrementally update the distance from
//		a kd-tree rectangle to the query point.  This comes about from
//		the fact that which each successive split, only one component
//		(along the dimension that is split) of the squared distance to
//		the child rectangle is different from the squared distance to
//		the parent rectangle.
//----------------------------------------------------------------------

//----------------------------------------------------------------------
//	annkSearch - search for the k nearest neighbors
//----------------------------------------------------------------------

void ANNkd_tree::annkSearch(
	ANNpoint			q,				// the query point
	int					k,				// number of near neighbors to return
	ANNidxArray			nn_idx,			// nearest neighbor indices (returned)
	ANNdistArray		dd,				// the approximate nearest neighbor
	ANNfp				eps)			// the error bound
{

	if (k > n_pts) {					// too many near neighbors?
		annError("Requesting more near neighbors than data points", ANNabort);
	}

	ANNfp max_err = ANN_POW(1.0f + eps);
	ANN_FLOP(2)							// increment floating op count

	ANNmin_k pointMK(k);
										// search starting at the root
	root->ann_search(annBoxDistance(q, bnd_box_lo, bnd_box_hi, dim), pts, q, max_err, dim, &pointMK);

	for (int i = 0; i < k; i++) {		// extract the k-th closest points
		dd[i] = pointMK.ith_smallest_key(i);
		nn_idx[i] = pointMK.ith_smallest_info(i);
	}
}

//----------------------------------------------------------------------
//	kd_split::ann_search - search a splitting node
//----------------------------------------------------------------------

void ANNkd_split::ann_search(ANNdist box_dist, ANNpointArray pts, ANNpoint q, ANNfp max_err, int dim, ANNmin_k * pointMK)
{
										// distance to cutting plane
	ANNcoord cut_diff = q[cut_dim] - cut_val;

	if (cut_diff < 0) {					// left of cutting plane
		child[ANN_LO]->ann_search(box_dist, pts, q, max_err, dim, pointMK);// visit closer child first

		ANNcoord box_diff = cd_bnds[ANN_LO] - q[cut_dim];
		if (box_diff < 0)				// within bounds - ignore
			box_diff = 0;
										// distance to further box
		box_dist = (ANNdist) ANN_SUM(box_dist,
				ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));

										// visit further child if close enough
		if (box_dist * max_err < pointMK->max_key())
			child[ANN_HI]->ann_search(box_dist, pts, q, max_err, dim, pointMK);

	}
	else {								// right of cutting plane
		child[ANN_HI]->ann_search(box_dist, pts, q, max_err, dim, pointMK);// visit closer child first

		ANNcoord box_diff = q[cut_dim] - cd_bnds[ANN_HI];
		if (box_diff < 0)				// within bounds - ignore
			box_diff = 0;
										// distance to further box
		box_dist = (ANNdist) ANN_SUM(box_dist,
				ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));

										// visit further child if close enough
		if (box_dist * max_err < pointMK->max_key())
			child[ANN_LO]->ann_search(box_dist, pts, q, max_err, dim, pointMK);

	}
	ANN_FLOP(10)						// increment floating ops
	ANN_SPL(1)							// one more splitting node visited
}

//----------------------------------------------------------------------
//	kd_leaf::ann_search - search points in a leaf node
//		Note: The unreadability of this code is the result of
//		some fine tuning to replace indexing by pointer operations.
//----------------------------------------------------------------------

void ANNkd_leaf::ann_search(ANNdist box_dist, ANNpointArray pts, ANNpoint q, ANNfp max_err, int dim, ANNmin_k * pointMK)
{
	register ANNdist dist;				// distance to data point
	register ANNcoord* pp;				// data coordinate pointer
	register ANNcoord* qq;				// query coordinate pointer
	register ANNdist min_dist;			// distance to k-th closest point
	register ANNcoord t;
	register int d;

	min_dist = pointMK->max_key(); // k-th smallest distance so far

	for (int i = 0; i < n_pts; i++) {	// check points in bucket

		pp = pts[bkt[i]];			// first coord of next data point
		qq = q;					// first coord of query point
		dist = 0;

#if ANNcoordPrec > 7
		for (d = 0; d < dim; d++) {
			ANN_COORD(1)				// one more coordinate hit
			ANN_FLOP(4)					// increment floating ops

			t = *(qq++) - *(pp++);		// compute length and adv coordinate
										// exceeds dist to k-th smallest?
			if( (dist = ANN_SUM(dist, ANN_POW(t))) > min_dist) {
				break;
			}
		}
#else
		_CRT_UNUSED(t);

		for (d = 0; d < dim; d += 8) {
			const __m128i qqq = _mm_loadu_si128((const __m128i *)qq);
			const __m128i ppp = _mm_loadu_si128((const __m128i *)pp);
			qq += 8;
			pp += 8;

			const __m128i qq_minus_pp = _mm_subs_epi16(qqq, ppp);
			const __m128i qq_minus_pp_sq = _mm_madd_epi16(qq_minus_pp, qq_minus_pp);

			const __m128i sumt = _mm_hadd_epi32(qq_minus_pp_sq, qq_minus_pp_sq);
			const __m128i sum = _mm_hadd_epi32(sumt, sumt);

			dist += sum.m128i_u32[0];

			if (dist > min_dist) {
				break;
			}
		}
#endif

		if (d >= dim &&					// among the k best?
		   (ANN_ALLOW_SELF_MATCH || dist!=0)) { // and no self-match problem
												// add it to the list
			pointMK->insert(dist, bkt[i]);
			min_dist = pointMK->max_key();
		}
	}
	ANN_LEAF(1)							// one more leaf node visited
	ANN_PTS(n_pts)						// increment points visited
}
