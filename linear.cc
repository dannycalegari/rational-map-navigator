/* linear.cc linear algebra operations */

/* 	Given a matrix A of full rank, and a column vector v with a single nonzero entry,
	we want to find A^{-1}(v). We do this by column operations. 
	
	First, after a permutation, we can assume that v is the vector (0,0, . . . , 1)^T.
	Start with matrix A and B=Id. 
	
	For j=0 to n-2,
		For i=0 to n-2-j,
			Find L_{j,i} such that A_{j,i} + L_{j,i}*A_{j,i+1} = 0.
			A_{*,i} -> A_{*,i}+L_{*,i}*A_{*,i+1}
			B_{*,i} -> B_{*,i}+L_{*,i}*B_{*,i+1}
			
	After doing this, the resulting matrix B is the matrix we want.	*/
 
	
	
	