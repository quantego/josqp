package com.quantego.josqp;

/**
 * A Java style implementation of a sparse matrix. Each row is represented by two arrays, an integer array of column indices and a double array of elements. 
 * Queries of elements are supported by a bisection search of column indices in each row.
 * 
 * @author Nils Loehndorf
 *
 */
public class SparseMatrix {
	final int _n;
	final Row[] _rows;
	
	public static class Row {
		
		final int _size;
		final int[] _index;
		final double[] _values;
		
		public Row(int[] index, double[] values) {
			if (index.length != values.length)
				throw new IllegalArgumentException("Values and index must have the same size.");
			_index = index;
			_values = values;
			_size = _index.length;
		}
		
		int getIndex(final int i) {
			int from = 0;
			int to = _size-1;
			while (from <= to) {
				int mid = (from + to) >>> 1;
				int midPoint = _index[mid];
				if (midPoint < i) {
					from = mid + 1;
				} else if (midPoint > i) {
					to = mid - 1;
				} else {
					return mid;
				}
			}
			return -(from + 1);
		}
		
		double getValue(int i) {
			int index = getIndex(i);
			return _values[index];
		}
	}
	
	public SparseMatrix(int rows) {
		_rows = new Row[rows];
		_n = _rows.length;
	}
	
	public void addRow(int i, int[] index, double[] values) {
		_rows[i] = new Row(index,values);
	}
	
	public double getValue(int i, int j) {
		return _rows[i].getValue(j);
	}
	
	public int[] getIndex(int i) {
		return _rows[i]._index;
	}
	
	public double[] getValues(int i) {
		return _rows[i]._values;
	}
	
	public SparseMatrix clone() {
		SparseMatrix clone = new SparseMatrix(_n);
		for (int i=0; i<_n; i++) {
			Row row = _rows[i];
			clone.addRow(i, row._index.clone(), row._values.clone());
		}
		return clone;
	}
	
	/*public void cholesky(double[] b) { 
		 SparseMatrix R = new SparseMatrix(_n); 
		 boolean isspd = (b.length == _n);
		 for (int j = 0; j < _n; j++) { 
			 double d = 0.0;
			 Row row = _rows[j];
			 int[] index = row._index;
			 double[] values = row._values;
			 int[] index2;
			 double[] values2;
			 for (int l = 0; l < row._size; l++) {
				 int k = index[l];
				 double s = getValue(k,j); 
				 for (int i = 0; i < k; i++) { 
					 s -= R.getValue(i, k)*R.getValue(i, j);
				 }
				 R[k][j] = s = s/R[k][k];
				 index2[l] = j
				 d += s*s; 
//				 isspd = isspd & (A[k][j] == A[j][k]); 
			 } 
			d = A[j][j] - d; 
			isspd = isspd & (d > 0.0); 
			R[j][j] = Math.sqrt(Math.max(d,0.0)); 
			for (int k = j+1; k < n; k++) { 
				  R[k][j] = 0.0; 
			} 
		} 
	}*/
	
}
