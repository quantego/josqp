package josqp;

import com.quantego.josqp.AMD;

public class AMDDemo {
	
	public static void main(String... args) {
		/* The symmetric can_24 Harwell/Boeing matrix, including upper and lower
	     * triangular parts, and the diagonal entries.  Note that this matrix is
	     * 0-based, with row and column indices in the range 0 to n-1. */
	    int n = 24, nz;
	    int Ap [ ] = { 0, 9, 15, 21, 27, 33, 39, 48, 57, 61, 70, 76, 82, 88, 94, 100,
		106, 110, 119, 128, 137, 143, 152, 156, 160 };
	    int Ai [ ] = {
		/* column  0: */    0, 5, 6, 12, 13, 17, 18, 19, 21,
		/* column  1: */    1, 8, 9, 13, 14, 17,
		/* column  2: */    2, 6, 11, 20, 21, 22,
		/* column  3: */    3, 7, 10, 15, 18, 19,
		/* column  4: */    4, 7, 9, 14, 15, 16,
		/* column  5: */    0, 5, 6, 12, 13, 17,
		/* column  6: */    0, 2, 5, 6, 11, 12, 19, 21, 23,
		/* column  7: */    3, 4, 7, 9, 14, 15, 16, 17, 18,
		/* column  8: */    1, 8, 9, 14,
		/* column  9: */    1, 4, 7, 8, 9, 13, 14, 17, 18,
		/* column 10: */    3, 10, 18, 19, 20, 21,
		/* column 11: */    2, 6, 11, 12, 21, 23,
		/* column 12: */    0, 5, 6, 11, 12, 23,
		/* column 13: */    0, 1, 5, 9, 13, 17,
		/* column 14: */    1, 4, 7, 8, 9, 14,
		/* column 15: */    3, 4, 7, 15, 16, 18,
		/* column 16: */    4, 7, 15, 16,
		/* column 17: */    0, 1, 5, 7, 9, 13, 17, 18, 19,
		/* column 18: */    0, 3, 7, 9, 10, 15, 17, 18, 19,
		/* column 19: */    0, 3, 6, 10, 17, 18, 19, 20, 21,
		/* column 20: */    2, 10, 19, 20, 21, 22,
		/* column 21: */    0, 2, 6, 10, 11, 19, 20, 21, 22,
		/* column 22: */    2, 20, 21, 22,
		/* column 23: */    6, 11, 12, 23 } ;

	    int[] P = new int[24], Pinv = new int[24];
	    int i, j, k, jnew, p, inew;
	    AMD.Status result ;
	    char[][] A = new char[24][24] ;



	    System.out.printf("AMD demo, with the 24-by-24 Harwell/Boeing matrix, can_24:\n") ;


	    /* print the input matrix */
	    nz = Ap [n] ;
	    System.out.printf("\nInput matrix:  %d-by-%d, with %d entries.\n"
		  + "   Note that for a symmetric matrix such as this one, only the\n"
		  + "   strictly lower or upper triangular parts would need to be\n"
		  + "   passed to AMD, since AMD computes the ordering of A+A'.  The\n"
		  + "   diagonal entries are also not needed, since AMD ignores them.\n"
		   , n, n, nz) ;
	    for (j = 0 ; j < n ; j++)
	    {
		System.out.printf ("\nColumn: %d, number of entries: %d, with row indices in"
			+" Ai [%d ... %d]:\n    row indices:",
			j, Ap [j+1] - Ap [j], Ap [j], Ap [j+1]-1) ;
		for (p = Ap [j] ; p < Ap [j+1] ; p++)
		{
		    i = Ai [p] ;
		    System.out.printf (" %d", i) ;
		}
		System.out.printf ("\n") ;
	    }

	    /* print a character plot of the input matrix.  This is only reasonable
	     * because the matrix is small. */
	    System.out.printf ("\nPlot of input matrix pattern:\n") ;
	    for (j = 0 ; j < n ; j++)
	    {
		for (i = 0 ; i < n ; i++) A [i][j] = '.' ;
		for (p = Ap [j] ; p < Ap [j+1] ; p++)
		{
		    i = Ai [p] ;
		    A [i][j] = 'X' ;
		}
	    }
	    System.out.printf ("    ") ;
	    for (j = 0 ; j < n ; j++) System.out.printf (" %1d", j % 10) ;
	    System.out.printf ("\n") ;
	    for (i = 0 ; i < n ; i++)
	    {
		System.out.printf ("%2d: ", i) ;
		for (j = 0 ; j < n ; j++)
		{
		    System.out.printf (" %c", A [i][j]) ;
		}
		System.out.printf ("\n") ;
	    }

	    /* order the matrix */
	    result = AMD.amd_order (n, Ap, Ai, P) ;
	    System.out.printf ("return value from amd_order: %s\n",
		result) ;

//	    /* print the statistics */
//	    amd_info (Info) ;

//	    if (result != AMD_OK)
//	    {
//		System.out.printf ("AMD failed\n") ;
//		exit (1) ;
//	    }

	    /* print the permutation vector, P, and compute the inverse permutation */
	    System.out.printf ("Permutation vector:\n") ;
	    for (k = 0 ; k < n ; k++)
	    {
		/* row/column j is the kth row/column in the permuted matrix */
		j = P [k] ;
		Pinv [j] = k ;
		System.out.printf (" %2d", j) ;
	    }
	    System.out.printf ("\n\n") ;

	    System.out.printf ("Inverse permutation vector:\n") ;
	    for (j = 0 ; j < n ; j++)
	    {
		k = Pinv [j] ;
		System.out.printf (" %2d", k) ;
	    }
	    System.out.printf ("\n\n") ;

	    /* print a character plot of the permuted matrix. */
	    System.out.printf ("\nPlot of permuted matrix pattern:\n") ;
	    for (jnew = 0 ; jnew < n ; jnew++)
	    {
		j = P [jnew] ;
		for (inew = 0 ; inew < n ; inew++) A [inew][jnew] = '.' ;
		for (p = Ap [j] ; p < Ap [j+1] ; p++)
		{
		    inew = Pinv [Ai [p]] ;
		    A [inew][jnew] = 'X' ;
		}
	    }
	    System.out.printf ("    ") ;
	    for (j = 0 ; j < n ; j++) System.out.printf (" %1d", j % 10) ;
	    System.out.printf ("\n") ;
	    for (i = 0 ; i < n ; i++)
	    {
		System.out.printf ("%2d: ", i) ;
		for (j = 0 ; j < n ; j++)
		{
		    System.out.printf (" %c", A [i][j]) ;
		}
		System.out.printf ("\n") ;
	    }

	}

}
