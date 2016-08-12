#ifndef _csr_matmatmult_h
#define _csr_matmatmult_h

#include <__complex.h>

/**
 * Computes C := alpha * B * A, in which A, B and C are CSR format
 * sparse matrices storing double precision real values.  B is
 * dimension m x p and A is dimension p x n.
 *
 * @param pCptr [OUT] pointer to the ptr array in the CSR structure of C
 * @param pCind [OUT] pointer to the ind array in the CSR structure of C
 * @param pCval [OUT] pointer to the val array in the CSR structure of C
 * @param pCnnz [OUT] pointer to the number of (structural) nonzeros in C
 * @param alpha [IN]  scalar multiplier
 * @param Bptr [IN]   ptr array in the CSR structure of B
 * @param Bind [IN]   ind array in the CSR structure of B
 * @param Bval [IN]   val array in the CSR structure of B
 * @param Aptr [IN]   ptr array in the CSR structure of A
 * @param Aind [IN]   ind array in the CSR structure of A
 * @param Aval [IN]   val array in the CSR structure of A
 * @param m [IN]      number of rows in C and B
 * @param p [IN]      number of columns in B, and number of rows in A
 * @param n [IN]      number of columns in C, and number of columns in A
 *
 * @return Error code (zero if no error)
 */
int
csr_matmatmult_double (int** pCptr, int** pCind, 
		       double** pCval, int* pCnnz,
		       double alpha,
		       int* Bptr, int* Bind, double* Bval,
		       int* Aptr, int* Aind, double* Aval,
		       const int m, const int p, const int n);

/**
 * Computes C := alpha * B * A, in which A, B and C are CSR format
 * sparse matrices storing double precision complex values.  B is
 * dimension m x p and A is dimension p x n.
 *
 * @param pCptr [OUT] pointer to the ptr array in the CSR structure of C
 * @param pCind [OUT] pointer to the ind array in the CSR structure of C
 * @param pCval [OUT] pointer to the val array in the CSR structure of C
 * @param pCnnz [OUT] pointer to the number of (structural) nonzeros in C
 * @param alpha [IN]  scalar multiplier
 * @param Bptr [IN]   ptr array in the CSR structure of B
 * @param Bind [IN]   ind array in the CSR structure of B
 * @param Bval [IN]   val array in the CSR structure of B
 * @param Aptr [IN]   ptr array in the CSR structure of A
 * @param Aind [IN]   ind array in the CSR structure of A
 * @param Aval [IN]   val array in the CSR structure of A
 * @param m [IN]      number of rows in C and B
 * @param p [IN]      number of columns in B, and number of rows in A
 * @param n [IN]      number of columns in C, and number of columns in A
 *
 * @return Error code (zero if no error)
 */
int
csr_matmatmult_complex (int** pCptr, int** pCind, 
			double_Complex** pCval, int* pCnnz,
			double_Complex alpha,
			int* Bptr, int* Bind, double_Complex* Bval,
			int* Aptr, int* Aind, double_Complex* Aval,
			const int m, const int p, const int n);
/**
 * Computes C := alpha * B * A, in which A, B and C are CSR format
 * sparse matrices storing a pattern (no values, just indices).  B is
 * dimension m x p and A is dimension p x n.
 *
 * @param pCptr [OUT] pointer to the ptr array in the CSR structure of C
 * @param pCind [OUT] pointer to the ind array in the CSR structure of C
 * @param pCnnz [OUT] pointer to the number of (structural) nonzeros in C
 * @param Bptr [IN]   ptr array in the CSR structure of B
 * @param Bind [IN]   ind array in the CSR structure of B
 * @param Aptr [IN]   ptr array in the CSR structure of A
 * @param Aind [IN]   ind array in the CSR structure of A
 * @param m [IN]      number of rows in C and B
 * @param p [IN]      number of columns in B, and number of rows in A
 * @param n [IN]      number of columns in C, and number of columns in A
 *
 * @return Error code (zero if no error)
 */
int
csr_matmatmult_pattern (int** pCptr, int** pCind, int* pCnnz,
			int* Bptr, int* Bind, 
			int* Aptr, int* Aind,
			const int m, const int p, const int n);


#endif /* _csr_matmatmult_h */
