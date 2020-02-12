#pragma once

void invert(float matrix[][MAX_ROW], int n)
{
	float m[(2*MAX_COLUMN)][MAX_ROW];		/*augmented matrix*/
	float i[MAX_COLUMN][MAX_ROW];				/*identity*/
	int cur_row;
	int cur_column;
	
	/*identity*/
	for (cur_row = 0; cur_row < n; cur_row++)
	{
		for (cur_column = 0; cur_column < n; cur_column++)
		{
			if (cur_row == cur_column)
			{
				i[cur_column][cur_row] = 1;
			}
			else
			{
				i[cur_column][cur_row] = 0;
			}
		}
	}
	
	/*augmented matrix, 1st half*/
	for (cur_row = 0; cur_row < n; cur_row++)
	{
		for (cur_column = 0; cur_column < n; cur_column++)
		{
			m[cur_column][cur_row] = matrix[cur_column][cur_row];
		}
	}
	/*augmented matrix, 2nd half*/
	for (cur_row = 0; cur_row < n; cur_row++)
	{
		for (cur_column = n; cur_column < (2 * n); cur_column++)
		{
			m[cur_column][cur_row] = i[(cur_column - n)][cur_row];
		}
	}
	
	echelon(m, n);
}

/*
 *	This is a subfunction of inverse where the augmented matrix is being reduced to echelon
 * form. If through the process of row operations a row zeros is obtained at the left side of
 * the augmented matrix or if any of the diagonal values is zero, it prints out an error 
 * message that the matrix has no inverse.
 *	Arguments:	
 *		m: the augmented matrix
 *		n: the size of the matrix being inverted
 *		inverse: an array to store the values of the inverted matrix
 *		temp: a temporary storage
 *
 */
 
void echelon(float m[(2*MAX_COLUMN)][MAX_ROW], int n)
{
	int i;										/*counter*/
	int cur_row;
	int cur_column;
	float inverse[MAX_COLUMN][MAX_ROW];
   float temp;  
												 
	
	/*pivoting and swapping of rows*/
	for (i = 0; i < (n-1); i++)
	{
		for (cur_row = i; cur_row < n; cur_row++)
		{
			int row_index_max = find_max_row(m, i, cur_row, n);
	
			for (cur_column = 0; cur_column < (2*n); cur_column++)
			{
				temp = m[cur_column][cur_row];
				m[cur_column][cur_row] = m[cur_column][row_index_max];
				m[cur_column][row_index_max] = temp;
			}

		}
		reduce_echelon(m, i, n);
		
		if (m[i][i] == 0)
		{
			printf("\nMatrix has no inverse!\n\n");
		}
																																																								
	}
	
	/*reverse row operations*/
	for (i = (n-1); i > 0;i--)
	{

		reverse_echelon(m, i, n);

	}
	
	/*reduces m to canonical form*/
	for (i = 0; i < n; i++)
	{
		canonical_form(m, i, n);
	}
	
	/*gets the inverse*/
	for (cur_row = 0; cur_row < n; cur_row++)
	{
		for (cur_column = 0; cur_column < n; cur_column++)
		{
			inverse[cur_column][cur_row] = m[(cur_column + n)][cur_row];
		}
	}
		print_matrix("\n----- INVERSE -----\n", inverse, n, n);

}

/*
 * This function is used to obtain the row which has the leading non-zero entry, an important
 * parameter in doing elementary row operations in obtaining the echelon form of an augmented
 * matrix.
 *	Arguments:
 * 	m: the augmented matrix
 *		cur_column, start_row, end: indices and counters
 * 	index_max_row: the index of the row which has the leading non-zero entry
 * Returns a value.
 *
 */

int find_max_row(float m[(2*MAX_COLUMN)][MAX_ROW], int cur_column, int start_row, int end)
{
	int index_max_row = start_row;
	int cur_row;
	
	for (cur_row = start_row; cur_row < end; cur_row++)
	{
		if (m[cur_column][cur_row] > m[cur_column][index_max_row])
		{
			index_max_row = cur_row;
		}
	}
	return index_max_row;
}

/*
 *	This function is used to get the echelon form of the augmented matrix.
 * Arguments:
 *		m: the augmented matrix
 *		index: the current index of the values of the augmented matrix being reduced to echelon
 *		n: the size of the matrix nxn being inverted
 *		temp: temporary storage
 *		factor: the factor used to reduce the values of the matrix into echelon
 *
 */
 
void reduce_echelon(float m[(2*MAX_COLUMN)][MAX_ROW], int index, int n)
{
	int cur_row;
	int cur_column;
	float temp;
	float factor;
	
	for (cur_row = index; cur_row < (n-1); cur_row++)
	{
		factor = m[index][(cur_row + 1)]/m[index][index];
		for (cur_column = index; cur_column < (n * 2); cur_column++)
		{
			temp = m[cur_column][(cur_row + 1)] - (factor * m[cur_column][index]);
			m[cur_column][(cur_row + 1)] = temp;
		}
	}
}

/*
 * This function is used for the reverse elementary row operations.
 *
 */
 
void reverse_echelon(float m[(2*MAX_COLUMN)][MAX_ROW], int index, int n)
{
	int cur_row;
	int cur_column;
	float temp;
	float factor;
		
	for (cur_row = index; cur_row > 0; cur_row--)
	{
		factor = m[index][(cur_row - 1)]/m[index][index];
		for (cur_column = index; cur_column < (2 * n); cur_column++)
		{
			temp = m[cur_column][(cur_row - 1)] - (factor * m[cur_column][index]);
			m[cur_column][(cur_row -1)] = temp;
		}
	}	
}

/*
 * This function is used to obtain the canonical form of the reduced echolon form of the 
 * augmented matrix.
 *
 */
 
void canonical_form(float m[(2*MAX_COLUMN)][MAX_ROW], int index, int n)
{
	int cur_column;
	float temp;
	float factor;
	
	factor = m[index][index];
	
	for (cur_column = index; cur_column < (n*2); cur_column++)
	{
		temp = m[cur_column][index]/factor; 
		m[cur_column][index] = temp;
	}
}
