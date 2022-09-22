package com.quantego.josqp;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * Sparse matrix implementation based on {@link RealMatrix} which stores the matrix in CRS (Yale) internally. 
 * @author Nils Loehndorf
 *
 */
public class CSCMatrixBuilder  {
	
	int _numRows, _numCols, _numElements;
	ArrayList<Double> _elements;
	ArrayList<Integer> _indices;
	int[] _starts;
	
	/**
	 * Creates a new empty matrix with a given number of columns. The matrix can start empty, but adding non-zero elements
	 * on a matrix of unknown dimension requires internal resize operation for each additional column added which incurs
	 * a minor overhead of copying the respective arrays. 
	 * @param numCols total number of (known) columns (e.g. variables)
	 */
	public CSCMatrixBuilder(int numCols) {
		_numCols = numCols;
		_elements = new ArrayList<>();
		_indices = new ArrayList<>();
		_starts = new int[numCols+1];
	};
	
	/**
	 * Set the element of the matrix. If the column index is equal or greater than the predefined number of columns, 
	 * the matrix will be automatically resized.
	 * @param row row index
	 * @param column column index
	 * @param value real value
	 */
	public CSCMatrixBuilder set(int row, int col, double value) {
		if (col >= _numCols)
			resize(col);
		if (row >= _numRows)
			_numRows = row+1;
		if (_starts[col] == _starts[col+1]) {
			insert(_starts[col], col, row, value);
			return this;
		}
		for (int i=_starts[col]; i<_starts[col+1]; i++) {
			int j = _indices.get(i);
			if (j>row) { //column gets inserted in row
				insert(i, col, row, value);
				return this;
			}
			if (j==row) { //column exists in row
				 _elements.set(i, value);
				 return this;
			}
		}
		//colum is last column in row
		insert(_starts[col+1], col, row, value);
		return this;
	}
	
	private void insert(int index, int col, int row, double value) {
		_indices.add(index, row);
		_elements.add(index, value);
		_numElements++;
		for (int k=col+1; k<_numCols+1; k++)
			_starts[k]++;
	}
	
	private void resize(int col) {
		_starts = Arrays.copyOfRange(_starts, 0, col+2);
		for (int i=_numCols+1; i<col+2; i++)
			_starts[i] = _starts[_numCols];
		_numCols = col+1;
	}

	public int getNumRows() {
		return _numRows;
	}
	
	public CSCMatrix build(boolean negate) {
		return new CSCMatrix(
				this._numRows,
				this._numCols,
				this._numElements,
    			this._starts,
				Utils.toIntArray(this._indices),
				negate ? Utils.toDoubleArrayNeg(this._elements) : Utils.toDoubleArray(this._elements)
		);
    }
	 

}
