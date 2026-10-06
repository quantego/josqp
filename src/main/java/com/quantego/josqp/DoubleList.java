package com.quantego.josqp;

import java.util.Arrays;

/** Growable primitive double list without boxing, used while a model is being built. */
final class DoubleList {
    private double[] _values;
    private int _size;

    DoubleList() {
        _values = new double[16];
    }

    DoubleList(int capacity) {
        _values = new double[Math.max(capacity, 1)];
    }

    void add(double value) {
        if (_size == _values.length)
            _values = Arrays.copyOf(_values, _size * 2);
        _values[_size++] = value;
    }

    void set(int index, double value) {
        if (index < 0 || index >= _size)
            throw new IndexOutOfBoundsException("Index " + index + " outside [0," + _size + ")");
        _values[index] = value;
    }

    double get(int index) {
        if (index < 0 || index >= _size)
            throw new IndexOutOfBoundsException("Index " + index + " outside [0," + _size + ")");
        return _values[index];
    }

    int size() {
        return _size;
    }

    double[] toArray() {
        return Arrays.copyOf(_values, _size);
    }

    double[] toArrayNegated() {
        double[] result = new double[_size];
        for (int i = 0; i < _size; i++)
            result[i] = -_values[i];
        return result;
    }
}
