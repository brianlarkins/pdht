#include "csr_matmatmult.h"
#include <smvm_malloc.h>
#include <stdio.h>
#include <assert.h>


int
csr_matmatmult_double (int** pCptr, int** pCind, double** pCval, int* pCnnz,
		       double alpha,
		       int* Bptr, int* Bind, double* Bval,
		       int* Aptr, int* Aind, double* Aval,
		       const int m, const int p, const int n)
{
  int i, j, k, l, mm;
  int* Cptr, *Cind;
  double* Cval;

  int* marker;
  int* previous_index;
  int Cnnz = 0;
  
  marker = (int*) smvm_malloc (n * sizeof (int));
  Cptr = (int*) smvm_malloc ((m+1) * sizeof (int));

  for (i = 0; i < n; i++)
    marker[i] = -1;
  
  Cptr[0] = Cnnz; /* = 0 */

  // for each row i of C
  for (i = 0; i < m; i++)
    {
      /* Go through the i-th row of B, using the index k.  This row
       * gives us the coefficients in the linear combination of the
       * rows of A.  Elements in the i-th row of B that are
       * (structurally) zero are skipped. */
      for (k = Bptr[i]; k < Bptr[i+1]; k++)
	{
	  /* Get the column index of the current coefficient
	   * (by which we will multiply the j-th row of A). */
	  j = Bind[k];
	  if (j < 0 || j >= p)
	    {
	      /* Column index out of range: clean up */
	      smvm_free (marker);
	      smvm_free (Cptr);
	      return -1; /* error flag */
	    }

	  /* Get the j-th row of A.  Any time we encounter a nonzero
	   * in this scan, there will be another nonzero in C. */
	  for (l = Aptr[j]; l < Aptr[j+1]; l++)
	    {
	      mm = Aind[l]; 
	      if (mm < 0 || mm >= n)
		{
		  /* Column index out of range: clean up */
		  smvm_free (marker);
		  smvm_free (Cptr);
		  return -2; /* error flag */
		}
	      /* Do this to avoid overcounting the number of nonzeros
	       * in C.  If we are currently computing row i of C, and
	       * we've already added C(i,mm) to C, then don't add it 
	       * again.  marker[mm] == i is true (on this i iteration)
	       * iff we have already added C(i,mm) to C. */
	      if (marker[mm] != i)
		{
		  marker[mm] = i; 
		  Cnnz++;
		}
	    }
	}
      Cptr[i+1] = Cnnz;
    }

  /* 
   * Now that we know how many nonzeros there will be in C, we can
   * allocate space for it. 
   */
  Cind = smvm_malloc (Cnnz * sizeof (int));
  Cval = smvm_calloc (Cnnz, sizeof (double));

  for (i = 0; i < n; i++)
    marker[i] = -1;

  previous_index = smvm_malloc (n * sizeof (int));
  for (i = 0; i < n; i++)
    previous_index[i] = -1;

  Cnnz = 0;
  for (i = 0; i < m; i++)
    {
      for (k = Bptr[i]; k < Bptr[i+1]; k++)
	{
	  const double b_coeff = Bval[k];
	  j = Bind[k];
	  for (l = Aptr[j]; l < Aptr[j+1]; l++)
	    {
	      mm = Aind[l];
	      if (marker[mm] != i)
		{
		  /* Remember that we've added a nonzero to C at (i,
		   * mm).  If while computing this i-th row of C, we
		   * encounter this column index mm again, then we
		   * don't need to increment Cnnz or store to Cind,
		   * because there is a structural nonzero already
		   * there.  (However, in that case we still need to
		   * increment the value of C(i,mm).  That's what the
		   * previous_index array is for.) */
		  marker[mm] = i;
		  /* Remember the location in the Cval array at which
		   * we stored this nonzero.  We'll need that location
		   * later, if while we are computing the i-th row of
		   * C, we encounter the same column index in a row of
		   * A. */
		  previous_index[mm] = Cnnz;
		  Cind[Cnnz] = mm;
		  Cval[Cnnz] = Cval[Cnnz] + alpha * b_coeff * Aval[l];
		  Cnnz++;
		}
	      else
		{
		  /* There's already a nonzero value in C at that
		   * location.  This means we need to add in the new
		   * nonzero value at that location, getting the
		   * location from the previous_index array. */
		  const int loc = previous_index[mm];
		  assert (loc != -1);
		  Cval[loc] = Cval[loc] + alpha * b_coeff * Aval[l];
		}
	    }
	}
    }

  /* Free up the temporary workspace */
  smvm_free (marker);
  smvm_free (previous_index);

  /* Return the matrix C through the output variables */
  *pCptr = Cptr;
  *pCind = Cind;
  *pCval = Cval;
  *pCnnz = Cnnz;

  /* No errors, so return zero */
  return 0;
}


int
csr_matmatmult_complex (int** pCptr, int** pCind, 
			double_Complex** pCval, int* pCnnz,
			double_Complex alpha,
			int* Bptr, int* Bind, double_Complex* Bval,
			int* Aptr, int* Aind, double_Complex* Aval,
			const int m, const int p, const int n)
{
  int i, j, k, l, mm;
  int* Cptr, *Cind;
  double_Complex* Cval;

  int* marker;
  int* previous_index;
  int Cnnz = 0;
  
  marker = (int*) smvm_malloc (n * sizeof (int));
  Cptr = (int*) smvm_malloc ((m+1) * sizeof (int));

  for (i = 0; i < n; i++)
    marker[i] = -1;
  
  Cptr[0] = Cnnz; /* = 0 */

  // for each row i of C
  for (i = 0; i < m; i++)
    {
      /* Go through the i-th row of B, using the index k.  This row
       * gives us the coefficients in the linear combination of the
       * rows of A.  Elements in the i-th row of B that are
       * (structurally) zero are skipped. */
      for (k = Bptr[i]; k < Bptr[i+1]; k++)
	{
	  /* Get the column index of the current coefficient
	   * (by which we will multiply the j-th row of A). */
	  j = Bind[k];
	  if (j < 0 || j >= p)
	    {
	      /* Column index out of range: clean up */
	      smvm_free (marker);
	      smvm_free (Cptr);
	      return -1; /* error flag */
	    }

	  /* Get the j-th row of A.  Any time we encounter a nonzero
	   * in this scan, there will be another nonzero in C. */
	  for (l = Aptr[j]; l < Aptr[j+1]; l++)
	    {
	      mm = Aind[l]; 
	      if (mm < 0 || mm >= n)
		{
		  /* Column index out of range: clean up */
		  smvm_free (marker);
		  smvm_free (Cptr);
		  return -2; /* error flag */
		}
	      /* Do this to avoid overcounting the number of nonzeros
	       * in C.  If we are currently computing row i of C, and
	       * we've already added C(i,mm) to C, then don't add it 
	       * again.  marker[mm] == i is true (on this i iteration)
	       * iff we have already added C(i,mm) to C. */
	      if (marker[mm] != i)
		{
		  marker[mm] = i; 
		  Cnnz++;
		}
	    }
	}
      Cptr[i+1] = Cnnz;
    }

  /* 
   * Now that we know how many nonzeros there will be in C, we can
   * allocate space for it. 
   */
  Cind = smvm_malloc (Cnnz * sizeof (int));
  Cval = smvm_calloc (Cnnz, sizeof (double_Complex));

  for (i = 0; i < n; i++)
    marker[i] = -1;

  previous_index = smvm_malloc (n * sizeof (int));
  for (i = 0; i < n; i++)
    previous_index[i] = -1;

  Cnnz = 0;
  for (i = 0; i < m; i++)
    {
      for (k = Bptr[i]; k < Bptr[i+1]; k++)
	{
	  const double_Complex b_coeff = Bval[k];
	  j = Bind[k];
	  for (l = Aptr[j]; l < Aptr[j+1]; l++)
	    {
	      mm = Aind[l];
	      if (marker[mm] != i)
		{
		  /* Remember that we've added a nonzero to C at (i,
		   * mm).  If while computing this i-th row of C, we
		   * encounter this column index mm again, then we
		   * don't need to increment Cnnz or store to Cind,
		   * because there is a structural nonzero already
		   * there.  (However, in that case we still need to
		   * increment the value of C(i,mm).  That's what the
		   * previous_index array is for.) */
		  marker[mm] = i;
		  /* Remember the location in the Cval array at which
		   * we stored this nonzero.  We'll need that location
		   * later, if while we are computing the i-th row of
		   * C, we encounter the same column index in a row of
		   * A. */
		  previous_index[mm] = Cnnz;
		  Cind[Cnnz] = mm;
		  Cval[Cnnz] = double_Complex_add (Cval[Cnnz], double_Complex_multiply (alpha, double_Complex_multiply (b_coeff, Aval[l])));
		  Cnnz++;
		}
	      else
		{
		  /* There's already a nonzero value in C at that
		   * location.  This means we need to add in the new
		   * nonzero value at that location, getting the
		   * location from the previous_index array. */
		  const int loc = previous_index[mm];
		  assert (loc != -1);
		  Cval[loc] = double_Complex_add (Cval[loc], double_Complex_multiply (alpha, double_Complex_multiply (b_coeff, Aval[l])));
		}
	    }
	}
    }

  /* Free up the temporary workspace */
  smvm_free (marker);
  smvm_free (previous_index);

  /* Return the matrix C through the output variables */
  *pCptr = Cptr;
  *pCind = Cind;
  *pCval = Cval;
  *pCnnz = Cnnz;

  /* No errors, so return zero */
  return 0;
}

int
csr_matmatmult_pattern (int** pCptr, int** pCind, int* pCnnz,
			int* Bptr, int* Bind, 
			int* Aptr, int* Aind, 
			const int m, const int p, const int n)
{
  int i, j, k, l, mm;
  int* Cptr, *Cind;

  int* marker;
  int Cnnz = 0;
  
  marker = (int*) smvm_malloc (n * sizeof (int));
  Cptr = (int*) smvm_malloc ((m+1) * sizeof (int));

  for (i = 0; i < n; i++)
    marker[i] = -1;
  
  Cptr[0] = Cnnz; /* = 0 */

  // for each row i of C
  for (i = 0; i < m; i++)
    {
      /* Go through the i-th row of B, using the index k.  This row
       * gives us the coefficients in the linear combination of the
       * rows of A.  Elements in the i-th row of B that are
       * (structurally) zero are skipped. */
      for (k = Bptr[i]; k < Bptr[i+1]; k++)
	{
	  /* Get the column index of the current coefficient
	   * (by which we will multiply the j-th row of A). */
	  j = Bind[k];
	  if (j < 0 || j >= p)
	    {
	      /* Column index out of range: clean up */
	      smvm_free (marker);
	      smvm_free (Cptr);
	      return -1; /* error flag */
	    }

	  /* Get the j-th row of A.  Any time we encounter a nonzero
	   * in this scan, there will be another nonzero in C. */
	  for (l = Aptr[j]; l < Aptr[j+1]; l++)
	    {
	      mm = Aind[l]; 
	      if (mm < 0 || mm >= n)
		{
		  /* Column index out of range: clean up */
		  smvm_free (marker);
		  smvm_free (Cptr);
		  return -2; /* error flag */
		}
	      /* Do this to avoid overcounting the number of nonzeros
	       * in C.  If we are currently computing row i of C, and
	       * we've already added C(i,mm) to C, then don't add it 
	       * again.  marker[mm] == i is true (on this i iteration)
	       * iff we have already added C(i,mm) to C. */
	      if (marker[mm] != i)
		{
		  marker[mm] = i; 
		  Cnnz++;
		}
	    }
	}
      Cptr[i+1] = Cnnz;
    }

  /* 
   * Now that we know how many nonzeros there will be in C, we can
   * allocate space for it. 
   */
  Cind = smvm_malloc (Cnnz * sizeof (int));

  for (i = 0; i < n; i++)
    marker[i] = -1;

  Cnnz = 0;
  for (i = 0; i < m; i++)
    {
      for (k = Bptr[i]; k < Bptr[i+1]; k++)
	{
	  j = Bind[k];
	  for (l = Aptr[j]; l < Aptr[j+1]; l++)
	    {
	      mm = Aind[l];
	      if (marker[mm] != i)
		{
		  /* Remember that we've added a nonzero to C at (i,
		   * mm).  If while computing this i-th row of C, we
		   * encounter this column index mm again, then we
		   * don't need to increment Cnnz or store to Cind,
		   * because there is a structural nonzero already
		   * there. */
		  marker[mm] = i;
		  Cind[Cnnz] = mm;
		  Cnnz++;
		}
	      /* No need for the else: this is a pattern matrix, so 
	       * we don't have to add to the nonzero value. */
	    }
	}
    }

  /* Free up the temporary workspace */
  smvm_free (marker);

  /* Return the matrix C through the output variables */
  *pCptr = Cptr;
  *pCind = Cind;
  *pCnnz = Cnnz;

  /* No errors, so return zero */
  return 0;
}

