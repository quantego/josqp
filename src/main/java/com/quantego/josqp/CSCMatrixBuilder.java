package com.quantego.josqp;

import java.util.Arrays;

/**
 * Sparse matrix builder that collects entries as primitive (row, column, value) triplets and
 * assembles a CSC matrix in one linear counting sort.
 *
 * <p>The builder never allocates a boxed object per entry. Assembly is stable by insertion order,
 * so setting the same (row, column) twice keeps the last value — the same semantics as the
 * previous map-based builder. Repeated {@link #update(boolean)} calls with unchanged entries and
 * the same sign reuse the assembled arrays.</p>
 *
 * @author Nils Loehndorf
 */
public class CSCMatrixBuilder {

    private int[] _row;
    private int[] _col;
    private double[] _value;
    private int _entries;
    private int _numRows, _numCols;

    private double[] _elements;
    private int[] _indices;
    private int[] _starts;
    private boolean _hasUpdate = true;
    private boolean _builtNegate;

    /**
     * Creates a new empty matrix. Rows and columns are grown on demand.
     */
    public CSCMatrixBuilder() {
        this(16);
    }

    /**
     * Creates a new empty matrix with room for the given number of entries.
     *
     * @param capacity expected number of entries
     */
    public CSCMatrixBuilder(int capacity) {
        int size = Math.max(capacity, 1);
        _row = new int[size];
        _col = new int[size];
        _value = new double[size];
    }

    private void grow() {
        int capacity = _row.length * 2;
        _row = Arrays.copyOf(_row, capacity);
        _col = Arrays.copyOf(_col, capacity);
        _value = Arrays.copyOf(_value, capacity);
    }

    /**
     * Set the element of the matrix. If the row or column index is equal to or greater than the
     * current dimension, the matrix is resized.
     *
     * @param row row index
     * @param col column index
     * @param value real value
     */
    public CSCMatrixBuilder set(int row, int col, double value) {
        if (row < 0 || col < 0)
            throw new IllegalArgumentException("Negative index in set(" + row + "," + col + ")");
        if (row >= _numRows)
            _numRows = row + 1;
        if (col >= _numCols)
            _numCols = col + 1;
        if (_entries == _row.length)
            grow();
        _row[_entries] = row;
        _col[_entries] = col;
        _value[_entries] = value;
        _entries++;
        _hasUpdate = true;
        return this;
    }

    /**
     * Assemblies the CSC arrays. Called implicitly by {@link #build(boolean)}.
     *
     * @param negate store the negated values (used for maximization models)
     */
    public void update(boolean negate) {
        if (!_hasUpdate && _elements != null && _builtNegate == negate)
            return;
        final int n = _entries;
        final int cols = _numCols;
        final int rows = _numRows;

        // Stable counting sort by row, then by column: result is ordered by (column, row),
        // ties keep insertion order so the last write to a cell wins.
        int[] order = new int[n];
        int[] sorted = new int[n];
        int[] rowCount = new int[Math.max(rows, 1)];
        int[] colStart = new int[cols + 1];
        for (int i = 0; i < n; i++) {
            rowCount[_row[i]]++;
            colStart[_col[i] + 1]++;
        }
        int acc = 0;
        for (int r = 0; r < rowCount.length; r++) {
            int count = rowCount[r];
            rowCount[r] = acc;
            acc += count;
        }
        for (int i = 0; i < n; i++)
            order[rowCount[_row[i]]++] = i;
        for (int c = 0; c < cols; c++)
            colStart[c + 1] += colStart[c];
        int[] colPos = Arrays.copyOf(colStart, cols + 1);
        for (int k = 0; k < n; k++) {
            int i = order[k];
            sorted[colPos[_col[i]]++] = i;
        }

        // Drop all but the last entry of every (column, row) run.
        int unique = 0;
        for (int k = 0; k < n; k++) {
            int i = sorted[k];
            if (k + 1 < n) {
                int next = sorted[k + 1];
                if (_col[next] == _col[i] && _row[next] == _row[i])
                    continue;
            }
            sorted[unique++] = i;
        }

        _elements = new double[unique];
        _indices = new int[unique];
        _starts = new int[cols + 1];
        for (int k = 0; k < unique; k++)
            _starts[_col[sorted[k]] + 1]++;
        for (int c = 0; c < cols; c++)
            _starts[c + 1] += _starts[c];
        for (int k = 0; k < unique; k++) {
            int i = sorted[k];
            _indices[k] = _row[i];
            _elements[k] = negate ? -_value[i] : _value[i];
        }

        _hasUpdate = false;
        _builtNegate = negate;
    }

    /** @return values of the assembled CSC matrix, sorted by column then row */
    public double[] getElements() {
        return _elements;
    }

    /** @return column pointers of the assembled CSC matrix (length columns + 1) */
    public int[] getStarts() {
        return _starts;
    }

    /** @return row indices of the assembled CSC matrix */
    public int[] getIndices() {
        return _indices;
    }

    /** @return number of rows */
    public int getNumRows() {
        return _numRows;
    }

    /** @return number of columns */
    public int getNumCols() {
        return _numCols;
    }

    /** @return number of entries currently collected (including duplicates) */
    public int getNumEntries() {
        return _entries;
    }

    /**
     * Builds the CSC matrix. The returned matrix owns the assembled arrays.
     *
     * @param negate store the negated values
     * @return CSC matrix
     */
    public CSCMatrix build(boolean negate) {
        update(negate);
        return new CSCMatrix(_numRows, _numCols, _elements.length, _starts, _indices, _elements);
    }
}
