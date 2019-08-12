#include "Array3D.h"

template <typename T, size_t width, size_t height, size_t depth>
Array3D<T, width, height, depth>::Array3D() {
	arr = new T[width * height * depth];
	this->_width = width;
	this->_height = height;
	this->_depth = depth;
	this->length = width * height * depth;
}

template <typename T, size_t width, size_t height, size_t depth>
Array3D<T, width, height, depth>::Array3D(T init) : Array3D() {
	for (size_t i = 0; i < length; i++)
		arr[i] = init;
}

template <typename T, size_t width, size_t height, size_t depth>
Array3D<T, width, height, depth>::~Array3D() {
	delete[] arr;
}

template <typename T, size_t width, size_t height, size_t depth>
T& Array3D<T, width, height, depth>::operator ()(size_t i, size_t j, size_t k) {
	return arr[i * _width * _depth + j * _depth + k]; 
}
