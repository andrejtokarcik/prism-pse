//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <dxp@cs.bham.uc.uk> (University of Birmingham)
//	* Rashid Mehmood <rxm@cs.bham.uc.uk> (University of Birmingham)
//	
//------------------------------------------------------------------------------
//	
//	This file is part of PRISM.
//	
//	PRISM is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//	
//	PRISM is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//	
//	You should have received a copy of the GNU General Public License
//	along with PRISM; if not, write to the Free Software Foundation,
//	Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//	
//==============================================================================

// includes
#include "PrismHybrid.h"
#include <math.h>
#include <util.h>
#include <cudd.h>
#include <dd.h>
#include <odd.h>
#include <dv.h>
#include "sparse.h"
#include "hybrid.h"
#include "PrismHybridGlob.h"
#include "jnipointer.h"

// local prototypes
static void sor_rec(HDDNode *hdd, int level, int row_offset, int col_offset, int r, int c, bool transpose);
static void sor_rm(RMSparseMatrix *rmsm, int row_offset, int col_offset, int r, int c, bool is_diag);
static void sor_cmsr(CMSRSparseMatrix *cmsrsm, int row_offset, int col_offset, int r, int c, bool is_diag);

// globals (used by local functions)
static HDDNode *zero;
static int num_levels;
static bool compact_d, compact_sm;
static double *sm_dist;
static int sm_dist_shift;
static int sm_dist_mask;
static double *diags_vec;
static DistVector *diags_dist;
static double *soln, *soln2;
static double x, sup_norm, omega;
static bool forwards; 

//------------------------------------------------------------------------------

// solve the linear equation system Ax=b with Gauss-Seidel/SOR

JNIEXPORT jlong __pointer JNICALL Java_hybrid_PrismHybrid_PH_1SOR
(
JNIEnv *env,
jclass cls,
jlong __pointer _odd,	// odd
jlong __pointer rv,	// row vars
jint num_rvars,
jlong __pointer cv,	// col vars
jint num_cvars,
jlong __pointer _a,	// matrix A
jlong __pointer _b,	// vector b (if null, assume all zero)
jlong __pointer _init,	// init soln
jboolean transpose,	// transpose A? (i.e. solve xA=b not Ax=b?)
jboolean row_sums,	// use row sums for diags instead? (strictly speaking: negative sum of non-diagonal row elements)
jdouble om,			// omega (over-relaxation parameter)
jboolean fwds		// forwards or backwards?
)
{
	// cast function parameters
	ODDNode *odd = jlong_to_ODDNode(_odd);		// odd
	DdNode **rvars = jlong_to_DdNode_array(rv); 	// row vars
	DdNode **cvars = jlong_to_DdNode_array(cv); 	// col vars
	DdNode *a = jlong_to_DdNode(_a);		// matrix A
	DdNode *b = jlong_to_DdNode(_b);		// vector b
	DdNode *init = jlong_to_DdNode(_init);		// init soln

	omega = om;
	forwards = fwds;

	// mtbdds
	DdNode *reach, *diags, *id, *tmp;
	// model stats
	int n;
	long nnz;
	// flags
	bool compact_b, l_b_max;
	// matrix mtbdd
	HDDMatrix *hddm;
	HDDNode *hdd;
	// vectors
	double *b_vec;
	DistVector *b_dist;
	// timing stuff
	long start1, start2, start3, stop;
	double time_taken, time_for_setup, time_for_iters;
	// misc
	int i, j, fb, l, h, i2, j2, fb2, l2, h2, iters;
	double d, kb, kbt;
	bool done, diag_done;
	
	// start clocks
	start1 = start2 = util_cpu_time();
	
	// get number of states
	n = odd->eoff + odd->toff;
	
	// get reachable states
	reach = odd->dd;
	
	// make local copy of a
	Cudd_Ref(a);
	
	// remove and keep diagonal entries of matrix A
	id =  DD_Identity(ddman, rvars, cvars, num_rvars);
	Cudd_Ref(reach);
	id = DD_And(ddman, id, reach);
	Cudd_Ref(id);
	Cudd_Ref(a);
	diags = DD_Apply(ddman, APPLY_TIMES, id, a);
	Cudd_Ref(id);
	a = DD_ITE(ddman, id, DD_Constant(ddman, 0), a);
	
	// build hdd for matrix
	PH_PrintToMainLog(env, "\nBuilding hybrid MTBDD matrix... ");
	hddm = build_hdd_matrix(a, rvars, cvars, num_rvars, odd, true, transpose);
	hdd = hddm->top;
	zero = hddm->zero;
	num_levels = hddm->num_levels;
	kb = hddm->mem_nodes;
	kbt = kb;
	PH_PrintToMainLog(env, "[levels=%d, nodes=%d] [%.1f KB]\n", hddm->num_levels, hddm->num_nodes, kb);
	
	// split hdd matrix into blocks
	// nb: in terms of memory, this gets precedence over sparse matrices
	PH_PrintToMainLog(env, "Splitting into blocks... ");
	split_hdd_matrix(hddm, compact, false, transpose);
	compact_b = hddm->compact_b;
	rearrange_hdd_blocks(hddm, false);
	kb = hddm->mem_b;
	kbt += kb;
	PH_PrintToMainLog(env, "[levels=%d, n=%d, nnz=%d%s] [%.1f KB]\n", hddm->l_b, hddm->blocks->n, hddm->blocks->nnz, compact_b?", compact":"", kb);
	
	// add sparse matrices
	PH_PrintToMainLog(env, "Adding explicit sparse matrices... ");
	add_sparse_matrices(hddm, compact, true, transpose);
	compact_sm = hddm->compact_sm;
	if (compact_sm) {
		sm_dist = hddm->dist;
		sm_dist_shift = hddm->dist_shift;
		sm_dist_mask = hddm->dist_mask;
	}
	l_b_max = (hddm->l_b == hddm->num_levels);
	kb = hddm->mem_sm;
	kbt += kb;
	PH_PrintToMainLog(env, "[levels=%d, num=%d%s] [%.1f KB]\n", hddm->l_sm, hddm->num_sm, compact_sm?", compact":"", kb);
	
	// get vector of diags, either by extracting from mtbdd or
	// by doing (negative, non-diagonal) row sums of original A matrix (and then setting to 1 if sum is 0)
	PH_PrintToMainLog(env, "Creating vector for diagonals... ");
	if (!row_sums) {
		diags = DD_MaxAbstract(ddman, diags, cvars, num_cvars);
		diags_vec = mtbdd_to_double_vector(ddman, diags, rvars, num_rvars, odd);
	} else {
		diags_vec = hdd_negative_row_sums(hddm, n, transpose);
	}
	// if any of the diagonals are zero, set them to one - avoids division by zero errors later
	// strictly speaking, such matrices shouldn't work for this iterative method
	// but they do occur, e.g. for steady-state computation of a bscc, this fixes it
	for (i = 0; i < n; i++) diags_vec[i] = (diags_vec[i] == 0) ? 1.0 : diags_vec[i];
	// try and convert to compact form if required
	compact_d = false;
	if (compact) {
		if (diags_dist = double_vector_to_dist(diags_vec, n)) {
			compact_d = true;
			free(diags_vec);
		}
	}
	kb = (!compact_d) ? n*8.0/1024.0 : (diags_dist->num_dist*8.0+n*2.0)/1024.0;
	kbt += kb;
	if (!compact_d) PH_PrintToMainLog(env, "[%.1f KB]\n", kb);
	else PH_PrintToMainLog(env, "[dist=%d, compact] [%.1f KB]\n", diags_dist->num_dist, kb);
	
	// invert diagonal
	if (!compact_d) {
		for (i = 0; i < n; i++) diags_vec[i] = 1.0 / diags_vec[i];
	} else {
		for (i = 0; i < diags_dist->num_dist; i++) diags_dist->dist[i] = 1.0 / diags_dist->dist[i];
	}
	
	// build b vector (if present)
	if (b != NULL) {
		PH_PrintToMainLog(env, "Creating vector for RHS... ");
		b_vec = mtbdd_to_double_vector(ddman, b, rvars, num_rvars, odd);
		// try and convert to compact form if required
		compact_b = false;
		if (compact) {
			if (b_dist = double_vector_to_dist(b_vec, n)) {
				compact_b = true;
				free(b_vec);
			}
		}
		kb = (!compact_b) ? n*8.0/1024.0 : (b_dist->num_dist*8.0+n*2.0)/1024.0;
		kbt += kb;
		if (!compact_b) PH_PrintToMainLog(env, "[%.1f KB]\n", kb);
		else PH_PrintToMainLog(env, "[dist=%d, compact] [%.1f KB]\n", b_dist->num_dist, kb);
	}
	
	// create solution/iteration vectors
	PH_PrintToMainLog(env, "Allocating iteration vectors... ");
	soln = mtbdd_to_double_vector(ddman, init, rvars, num_rvars, odd);
	soln2 = (double*)calloc(hddm->blocks->max, sizeof(double));
	if (!soln2) fatal(" soln2 buffer allocation problem");
	kb = (n*8.0/1024.0)+(hddm->blocks->max*8.0/1024.0);
	kbt += kb;
	PH_PrintToMainLog(env, "[%.1f + %.1f = %.1f KB]\n", (n*8.0/1024.0), (hddm->blocks->max*8.0/1024.0), kb);
	
	// print total memory usage
	PH_PrintToMainLog(env, "TOTAL: [%.1f KB]\n", kbt);
	
	// get setup time
	stop = util_cpu_time();
	time_for_setup = (double)(stop - start2)/1000;
	start2 = stop;
	
	// start iterations
	iters = 0;
	done = false;
	PH_PrintToMainLog(env, "\nStarting iterations...\n");
	
	while (!done && iters < max_iters) {
		
		iters++;
		
//		PH_PrintToMainLog(env, "Iteration %d: ", iters);
//		start3 = util_cpu_time();
		
		sup_norm = 0.0;
		
		// stuff for block storage
		int b_n = hddm->blocks->n;
		int b_nnz = hddm->blocks->nnz;
		HDDNode **b_blocks = hddm->blocks->blocks;
		unsigned int *b_rowscols = hddm->blocks->rowscols;
		unsigned char *b_counts = hddm->blocks->counts;
		int *b_starts = (int *)hddm->blocks->counts;
		bool b_use_counts = hddm->blocks->use_counts;
		int *b_offsets = hddm->blocks->offsets;
		HDDNode **b_nodes = hddm->row_tables[hddm->l_b];
		int b_dist_shift = hddm->blocks->dist_shift;
		int b_dist_mask = hddm->blocks->dist_mask;
		int row_offset, col_offset;
		HDDNode *node;
		
		// loop through rows of blocks
		l = b_nnz; h = 0;
		for(fb = 0; fb < b_n; fb++)
		{
			// loop actually over i (can do forwards or backwards sor/gs)
			i = (forwards) ? fb : b_n-1-fb;
			
			// store block row offset
			row_offset = b_offsets[i];
			
			// initialise (partial) solution vector
			h2 = b_offsets[i+1] - b_offsets[i];
			// initialise vector
			if (b == NULL) {
				for (i2 = 0; i2 < h2; i2++) { soln2[i2] = 0.0; }
			} else if (!compact_b) {
				for (i2 = 0; i2 < h2; i2++) { soln2[i2] = b_vec[row_offset + i2]; }
			} else {
				for (i2 = 0; i2 < h2; i2++) { soln2[i2] = b_dist->dist[b_dist->ptrs[row_offset + i2]]; }
			}
			
			// loop through blocks in this row of blocks
			if (!b_use_counts) { l = b_starts[i]; h = b_starts[i+1]; }
			else if (forwards) { l = h; h += b_counts[i]; }
			else { h = l; l -= b_counts[i]; }
			diag_done = false;
			for(j = l; j < h; j++)
			{
				// get node for block and its col offset
				if (!compact_b) {
					node = b_blocks[j];
					col_offset = b_offsets[b_rowscols[j]];
				} else {
					node = b_nodes[(int)(b_rowscols[j] & b_dist_mask)];
					col_offset = b_offsets[(int)(b_rowscols[j] >> b_dist_shift)];
				}
				
				// trivial case where we are the bottom of the mtbdd already
				if (l_b_max) {
					soln2[0] -= soln[col_offset] * node->type.val;
					//printf("(%d,%d)=%f\n", row_offset, col_offset, node->type.val);
					continue;
				}
				
				// non-diagonal blocks treated normally
				// (diagonal should be the last block, unless it is absent because empty)
				if ((j != h-1) || (j == h-1 && row_offset !=col_offset)) {
					sor_rec(node, hddm->l_b, row_offset, col_offset, 0, 0, transpose);
				}
				// diagonal blocks (last blocks in row/col) are different
				// call sparse matrix traversal directly with "is_diag" flag = true
				else {
					diag_done = true;
					if (!compact_sm) {
						sor_rm((RMSparseMatrix *)node->sm.ptr, row_offset, col_offset, 0, 0, true);
					} else {
						sor_cmsr((CMSRSparseMatrix *)node->sm.ptr, row_offset, col_offset, 0, 0, true);
					}
				}
			}
			
			// if we never found a diagonal block (because it is empty and so not there),
			// then we do the stuff that should have been sone after the processing of the diagonal block
			for (i2 = 0; i2 < h2; i2++) {
				// divide by diagonal
				if (!compact_d) {
					soln2[i2] *= diags_vec[row_offset + i2];
				} else {
					soln2[i2] *= (diags_dist->dist[(int)diags_dist->ptrs[row_offset + i2]]);
				}
				// do over-relaxation if necessary
				if (omega != 1) {
					soln2[i2] = ((1-omega) * soln[row_offset + i2]) + (omega * soln2[i2]);
				}       
				// compute norm for convergence
				x = fabs(soln2[i2] - soln[row_offset + i2]);
				if (term_crit == TERM_CRIT_RELATIVE) {
					x /= soln2[i2];
				}
				if (x > sup_norm) sup_norm = x;
				// set vector element
				soln[row_offset + i2] = soln2[i2];
			}

			// trivial case where we are the bottom of the mtbdd already
			if (l_b_max) {
				soln2[0] *= ((!compact_d)?(diags_vec[row_offset]):(diags_dist->dist[(int)diags_dist->ptrs[row_offset]]));
				if (omega != 1) soln2[0] = ((1-omega) * soln[row_offset]) + (omega * soln2[0]);
				x = fabs(soln2[0] - soln[row_offset]);
				if (term_crit == TERM_CRIT_RELATIVE) x /= soln2[0];
				if (x > sup_norm) sup_norm = x;
				soln[row_offset] = soln2[0];
			}
		}
		
		// check convergence
		if (sup_norm < term_crit_param) {
			done = true;
		}
		
//		PH_PrintToMainLog(env, "%.2f %.2f sec\n", ((double)(util_cpu_time() - start3)/1000), ((double)(util_cpu_time() - start2)/1000)/iters);
	}
	
	// stop clocks
	stop = util_cpu_time();
	time_for_iters = (double)(stop - start2)/1000;
	time_taken = (double)(stop - start1)/1000;
	
	// print iters/timing info
	PH_PrintToMainLog(env, "\n%s%s: %d iterations in %.2f seconds (average %.6f, setup %.2f)\n", forwards?"":"Backwards ", (omega == 1.0)?"Gauss-Seidel":"SOR", iters, time_taken, time_for_iters/iters, time_for_setup);
	
	// free memory
	Cudd_RecursiveDeref(ddman, a);
	Cudd_RecursiveDeref(ddman, id);
	Cudd_RecursiveDeref(ddman, diags);
	free_hdd_matrix(hddm);
	if (compact_d) free_dist_vector(diags_dist); else free(diags_vec);
	if (b != NULL) if (compact_b) free_dist_vector(b_dist); else free(b_vec);
	free(soln2);
	
	// if the iterative method didn't terminate, this is an error
	if (!done) { delete soln; PH_SetErrorMessage("Iterative method did not converge within %d iterations.\nConsider using a different numerical method or increasing the maximum number of iterations", iters); return 0; }
	
	return ptr_to_jlong(soln);
}

//------------------------------------------------------------------------------

static void sor_rec(HDDNode *hdd, int level, int row_offset, int col_offset, int r, int c, bool transpose)
{
	HDDNode *e, *t;
	
	// if it's the zero node
	if (hdd == zero) {
		return;
	}
	// or if we've reached a submatrix
	// (check for non-null ptr but, equivalently, we could just check if level==l_sm)
	else if (hdd->sm.ptr) {
		if (!compact_sm) {
			sor_rm((RMSparseMatrix *)hdd->sm.ptr, row_offset, col_offset, r, c, false);
		} else {
			sor_cmsr((CMSRSparseMatrix *)hdd->sm.ptr, row_offset, col_offset, r, c, false);
		}
		return;
	}
	// or if we've reached the bottom
	else if (level == num_levels) {
		//printf("(%d,%d)=%f\n", row_offset, col_offset, hdd->type.val);
		soln2[r] -= soln[col_offset+c] * hdd->type.val;
		return;
	}
	// otherwise recurse
	e = hdd->type.kids.e;
	if (e != zero) {
		if (!transpose) {
			sor_rec(e->type.kids.e, level+1, row_offset, col_offset, r, c, transpose);
			sor_rec(e->type.kids.t, level+1, row_offset, col_offset, r, c+e->off.val, transpose);
		} else {
			sor_rec(e->type.kids.e, level+1, row_offset, col_offset, r, c, transpose);
			sor_rec(e->type.kids.t, level+1, row_offset, col_offset, r+e->off.val, c, transpose);
		}
	}
	t = hdd->type.kids.t;
	if (t != zero) {
		if (!transpose) {
			sor_rec(t->type.kids.e, level+1, row_offset, col_offset, r+hdd->off.val, c, transpose);
			sor_rec(t->type.kids.t, level+1, row_offset, col_offset, r+hdd->off.val, c+t->off.val, transpose);
		} else {
			sor_rec(t->type.kids.e, level+1, row_offset, col_offset, r, c+hdd->off.val, transpose);
			sor_rec(t->type.kids.t, level+1, row_offset, col_offset, r+t->off.val, c+hdd->off.val, transpose);
		}
	}
}

//-----------------------------------------------------------------------------------

static void sor_rm(RMSparseMatrix *rmsm, int row_offset, int col_offset, int r, int c, bool is_diag)
{
	int fb2, i2, j2, l2, h2;
	int sm_n = rmsm->n;
	int sm_nnz = rmsm->nnz;
	double *sm_non_zeros = rmsm->non_zeros;
	unsigned char *sm_row_counts = rmsm->row_counts;
	int *sm_row_starts = (int *)rmsm->row_counts;
	bool sm_use_counts = rmsm->use_counts;
	unsigned int *sm_cols = rmsm->cols;
	
	// loop through rows of submatrix
	l2 = sm_nnz; h2 = 0;
	for (fb2 = 0; fb2 < sm_n; fb2++) {
		
		// loop actually over i2 (can do forwards or backwards sor/gs)
		i2 = (forwards) ? fb2 : sm_n-1-fb2;
		
		// loop through entries in this row
		if (!sm_use_counts) { l2 = sm_row_starts[i2]; h2 = sm_row_starts[i2+1]; }
		else if (forwards) { l2 = h2; h2 += sm_row_counts[i2]; }
		else { h2 = l2; l2 -= sm_row_counts[i2]; }
		for (j2 = l2; j2 < h2; j2++) {
			soln2[r + i2] -= soln[col_offset + c + sm_cols[j2]] * sm_non_zeros[j2];
			//printf("(%d,%d)=%f\n", r + i2, col_offset + c + sm_cols[j2], sm_non_zeros[j2]);
		}
		
		if (is_diag) {
			// divide by diagonal
			if (!compact_d) {
				soln2[r + i2] *= diags_vec[row_offset + r + i2];
			} else {
				soln2[r + i2] *= (diags_dist->dist[(int)diags_dist->ptrs[row_offset + r + i2]]);
			}
			// do over-relaxation if necessary
			if (omega != 1) {
				soln2[r + i2] = ((1-omega) * soln[row_offset + r + i2]) + (omega * soln2[r + i2]);
			}	
			// compute norm for convergence
			x = fabs(soln2[r + i2] - soln[row_offset + r + i2]);
			if (term_crit == TERM_CRIT_RELATIVE) {
				x /= soln2[r + i2];
			}
			if (x > sup_norm) sup_norm = x;
			// set vector element
			soln[row_offset + r + i2] = soln2[r + i2];
		}
	}
}

//-----------------------------------------------------------------------------------

static void sor_cmsr(CMSRSparseMatrix *cmsrsm, int row_offset, int col_offset, int r, int c, bool is_diag)
{
	int fb2, i2, j2, l2, h2;
	int sm_n = cmsrsm->n;
	int sm_nnz = cmsrsm->nnz;
	unsigned char *sm_row_counts = cmsrsm->row_counts;
	int *sm_row_starts = (int *)cmsrsm->row_counts;
	bool sm_use_counts = cmsrsm->use_counts;
	unsigned int *sm_cols = cmsrsm->cols;
	
	// loop through rows of submatrix
	l2 = sm_nnz; h2 = 0;
	for (fb2 = 0; fb2 < sm_n; fb2++) {
		
		// loop actually over i2 (can do forwards or backwards sor/gs)
		i2 = (forwards) ? fb2 : sm_n-1-fb2;
		
		// loop through entries in this row
		if (!sm_use_counts) { l2 = sm_row_starts[i2]; h2 = sm_row_starts[i2+1]; }
		else if (forwards) { l2 = h2; h2 += sm_row_counts[i2]; }
		else { h2 = l2; l2 -= sm_row_counts[i2]; }
		for (j2 = l2; j2 < h2; j2++) {
			soln2[r + i2] -= soln[col_offset + c + (int)(sm_cols[j2] >> sm_dist_shift)] * sm_dist[(int)(sm_cols[j2] & sm_dist_mask)];
			//printf("(%d,%d)=%f\n", row_offset + r + i2, col_offset + c + (int)(sm_cols[j2] >> sm_dist_shift), sm_dist[(int)(sm_cols[j2] & sm_dist_mask)]);
		}
		
		if (is_diag) {
			// divide by diagonal
			if (!compact_d) {
				soln2[r + i2] *= diags_vec[row_offset + r + i2];
			} else {
				soln2[r + i2] *= (diags_dist->dist[(int)diags_dist->ptrs[row_offset + r + i2]]);
			}
			// do over-relaxation if necessary
			if (omega != 1) {
				soln2[r + i2] = ((1-omega) * soln[row_offset + r + i2]) + (omega * soln2[r + i2]);
			}	
			// compute norm for convergence
			x = fabs(soln2[r + i2] - soln[row_offset + r + i2]);
			if (term_crit == TERM_CRIT_RELATIVE) {
				x /= soln2[r + i2];
			}
			if (x > sup_norm) sup_norm = x;
			// set vector element
			soln[row_offset + r + i2] = soln2[r + i2];
		}
	}
}

//------------------------------------------------------------------------------
