package com.quantego.josqp;

import java.util.Arrays;

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
	
	public SparseMatrix(int n) {
		_rows = new Row[n];
		_n = n;
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
	
	public class CholeskyDecomposition {
		final SparseMatrix L;
		final double[] d;
		final double[] dinv;
		public CholeskyDecomposition(SparseMatrix l, double[] d, double[] dinv) {
			L = l;
			this.d = d;
			this.dinv = dinv;
		}
	}
	
	public CholeskyDecomposition cholesky() { 
		 SparseMatrix L = new SparseMatrix(_n);
		 double[] D = new double[_n];
		 double[] Dinv = new double[_n];
		 for (int j = 0; j < _n; j++) { 
			 double d = 0.0;
			 Row ArowJ = _rows[j];
			 int[] ArowJindex = ArowJ._index;
			 double[] ArowJvalues = ArowJ._values;
			 int[] LrowJindex = ArowJ._index;
			 double[] LrowJvalues = new double[ArowJ._size];
			 double Ajj = 0.0;
			 for (int k = 0; k < ArowJindex.length; k++) {
				 if (ArowJindex[k]<j) {
					 double s = ArowJvalues[k]; 
					 Row LrowK = L._rows[k];
					 double[] LrowKvalues = LrowK._values;				 
					 for (int i = 0; i < LrowJvalues.length; i++) {
						 if (LrowJindex[i]<LrowJindex[k])
							 s -= LrowKvalues[k]*LrowJvalues[i];
						 else
							 break;
					 }
					 LrowJvalues[k] = s = s*Dinv[k];
					 d += s*s; 
				 }
				 else {
					 if (ArowJindex[k]==j) Ajj = ArowJvalues[k];
					 break;
				 } 
			 }
			 L.addRow(j, LrowJindex, LrowJvalues);
			 D[j] = Math.sqrt(Ajj - d);
			 Dinv[j] = 1./D[j];
		} 
		return new CholeskyDecomposition(L,D,Dinv);
	}
	
	public static void main(String[] args){
		SparseMatrix A = new SparseMatrix(3);
		A.addRow(0, new int[] {0, 1, 2},new double[]{25, 15, -5});
		A.addRow(1, new int[] {0, 1}, new double[] {15,18});
		A.addRow(2, new int[] {0, 2}, new double[] {-5,11});
		CholeskyDecomposition C = A.cholesky();
		System.out.println(Arrays.toString(C.L._rows[0]._index));
		System.out.println(Arrays.toString(C.L._rows[0]._values));
		System.out.println(Arrays.toString(C.L._rows[1]._index));
		System.out.println(Arrays.toString(C.L._rows[1]._values));
		System.out.println(Arrays.toString(C.L._rows[2]._index));
		System.out.println(Arrays.toString(C.L._rows[2]._values));
	}
	
}
