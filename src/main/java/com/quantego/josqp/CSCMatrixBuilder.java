package com.quantego.josqp;

import java.util.*;

/**
 * Sparse matrix implementation which stores the matrix in CRS (Yale) internally.
 * @author Nils Loehndorf
 *
 */
public class CSCMatrixBuilder  {

	int _numRows, _numCols, _numElements;
	ArrayList<Map<Integer,Double>> _values;
	private double[] _elements;
	private int[] _indices;
	private int[] _starts;
	private boolean _hasUpdate = false;
	
	/**
	 * Creates a new empty matrix with a given number of columns. The matrix can start empty, but adding non-zero elements
	 * on a matrix of unknown dimension requires internal resize operation for each additional column added which incurs
	 * a minor overhead of copying the respective arrays. 
	 * @param numCols total number of (known) columns (e.g. variables)
	 */
	public CSCMatrixBuilder() {
		_numCols = 0;
		_values = new ArrayList<>();
		_starts = new int[1];
	};
	
	/**
	 * Set the element of the matrix. If the column index is equal or greater than the predefined number of columns, 
	 * the matrix will be automatically resized.
	 * @param row row index
	 * @param col column index
	 * @param value real value
	 */
	/*public CSCMatrixBuilder set(Integer row, int col, Double value) {
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
	}*/
	public CSCMatrixBuilder set(int row, int col, double value) {
		if (row>=_numRows)
			_numRows = row+1;
		if (col>=_numCols) {
			for (int i=_numCols; i<=col; i++)
				_values.add(new TreeMap<>());
			_numCols = col+1;
		}
		if(_values.get(col).put(row,value)==null)
			_numElements++;
		_hasUpdate = true;
		return this;
	}

	public void update(boolean negate) {
		_elements = new double[_numElements];
		_indices = new int[_numElements];
		_starts = new int[_numCols+1];
		int n=0;
		for (int j=0; j<_numCols; j++) {
			Map<Integer,Double> map = _values.get(j);
			_starts[j+1] = _starts[j] + map.size();
			for (Map.Entry<Integer,Double> e : map.entrySet()) {
				_elements[n] = negate?-e.getValue():e.getValue();
				_indices[n] = e.getKey();
				n++;
			}
		}
		_hasUpdate = false;
	}

	public double[] getElements() {
		return _elements;
	}

	public int[] getStarts() {
		return _starts;
	}

	public int[] getIndices() {
		return _indices;
	}
	
	/*private void insert(int index, int col, Integer row, Double value) {
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
	}*/

	public int getNumRows() {
		return _numRows;
	}
	
	public CSCMatrix build(boolean negate) {
		update(negate);
		return new CSCMatrix(
				this._numRows,
				this._numCols,
				this._numElements,
    			this._starts,
				this._indices,
				this._elements
		);
    }
	 

}
