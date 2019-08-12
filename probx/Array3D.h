#pragma once

template <typename T, size_t width, size_t height, size_t depth>
class Array3D {

	T* arr;

public:
	size_t _width;
	size_t _height;
	size_t _depth;
	size_t length;

	Array3D();
	Array3D(T init) : Array3D();
	~Array3D();

	T& operator ()(size_t i, size_t j, size_t k);
};

