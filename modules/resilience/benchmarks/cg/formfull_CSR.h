#include <cassert>

/*! \brief
 * Form a full sparse matrix, in CSR format, from a half (upper triangular) 
 * matrix also read in CSR format 
 * <pre>
 * On input, nonz/nzval/colind/rowptr represents upper part of a symmetric
 * matrix. On exit, it represents the full matrix with lower and upper parts.
 * </pre>
 */
static void
FormFullA(int n, int *nonz, double **nzval, int **colind, int **rowptr)
{
    register int i, j, k, row, new_nnz;
    int *t_colind, *t_rowptr, *al_colind, *al_rowptr, *a_colind, *a_rowptr;
    int *marker;
    double *t_val, *al_val, *a_val;

    al_colind = *colind;
    al_rowptr = *rowptr;
    al_val = *nzval;

    assert ( (marker =(int *) malloc( (n+1) * sizeof(int)) ) );

    assert ( (t_rowptr = (int *) malloc( (n+1) * sizeof(int)) ) );

    assert ( (t_colind = (int *) malloc( *nonz * sizeof(int)) ) );

    assert ( (t_val = (double*) malloc( *nonz * sizeof(double)) ) );


    /* Get counts of each row of T, and set up row pointers */
    for (i = 0; i < n; ++i) marker[i] = 0;
    for (j = 0; j < n; ++j) {
	for (i = al_rowptr[j]; i < al_rowptr[j+1]; ++i)
	    ++marker[al_colind[i]];
    }
    t_rowptr[0] = 0;
    for (i = 0; i < n; ++i) {
	t_rowptr[i+1] = t_rowptr[i] + marker[i];
	marker[i] = t_rowptr[i];
    }

    /* Transpose matrix A to T */
    for (j = 0; j < n; ++j)
	for (i = al_rowptr[j]; i < al_rowptr[j+1]; ++i) {
	    row = al_colind[i];
	    t_colind[marker[row]] = j;
	    t_val[marker[row]] = al_val[i];
	    ++marker[row];
	}

    new_nnz = *nonz * 2 - n;
    assert ( (a_rowptr = (int *) malloc( (n+1) * sizeof(int)) ) );

    assert ( (a_colind = (int *) malloc( new_nnz * sizeof(int)) ) );

    assert ( (a_val = (double*) malloc( new_nnz * sizeof(double)) ) );

    
    a_rowptr[0] = 0;
    k = 0;
    for (j = 0; j < n; ++j) {
      for (i = t_rowptr[j]; i < t_rowptr[j+1]; ++i) {
	if ( t_colind[i] != j ) { /* not diagonal */
	  a_colind[k] = t_colind[i];
	  a_val[k] = t_val[i];
#ifdef DEBUG
	  if ( fabs(a_val[k]) < 4.047e-300 )
	      printf("%5d: %e\n", k, a_val[k]);
#endif
	  ++k;
	}
      }

      for (i = al_rowptr[j]; i < al_rowptr[j+1]; ++i) {
	a_colind[k] = al_colind[i];
	a_val[k] = al_val[i];
#ifdef DEBUG
	if ( fabs(a_val[k]) < 4.047e-300 )
	    printf("%5d: %e\n", k, a_val[k]);
#endif
	++k;
      }
      
      a_rowptr[j+1] = k;
    }

    printf("FormFullA: new_nnz = %d, k = %d\n", new_nnz, k);

    free(al_val);
    free(al_colind);
    free(al_rowptr);
    free(marker);
    free(t_val);
    free(t_colind);
    free(t_rowptr);

    *nzval = a_val;
    *colind = a_colind;
    *rowptr = a_rowptr;
    *nonz = new_nnz;
}
