package com.quantego.josqp;

import java.util.Arrays;

/** Growable primitive int list without boxing, used while a model is being built. */
final class IntList {
    private int[] _values;
    private int _size;

    IntList() {
        _values = new int[16];
    }

    IntList(int capacity) {
        _values = new int[Math.max(capacity, 1)];
    }

    void add(int value) {
        if (_size == _values.length)
            _values = Arrays.copyOf(_values, _size * 2);
        _values[_size++] = value;
    }

    void set(int index, int value) {
        if (index < 0 || index >= _size)
            throw new IndexOutOfBoundsException("Index " + index + " outside [0," + _size + ")");
        _values[index] = value;
    }

    int get(int index) {
        if (index < 0 || index >= _size)
            throw new IndexOutOfBoundsException("Index " + index + " outside [0," + _size + ")");
        return _values[index];
    }

    int size() {
        return _size;
    }

    int[] toArray() {
        return Arrays.copyOf(_values, _size);
    }
}
