//#include <iostream>
//#include <cmath>
//
//using namespace std;
//
//enum AlgorithmType {
//	FirstFit,
//	NextFit
//};
//AlgorithmType bin_algorithm;
//
//// The maximum input weight U~[1, M]
//double M = 6;
//
//template <typename T>
//class Array3D {
//
//	/*
//	This class implements a 3D Array in a 1D Array
//	Where the access to a cell is ARRAY_NAME(X, Y, Z)
//	*/
//
//	T* arr;
//
//public:
//	size_t width;
//	size_t height;
//	size_t depth;
//	size_t length;
//
//	// Default CTOR, will create a 1 cell array
//	Array3D() : arr(nullptr), width(0), height(0), depth(0), length(0) {}
//
//	// Creates a 3D array with size width = w, height = h, depth = d (w*h*d)
//	Array3D(size_t w, size_t h, size_t d) {
//		this->width = w;
//		this->height = h;
//		this->depth = d;
//		this->length = w * h * d;
//		arr = new T[this->length];
//	}
//
//	// Creates a 3D array and initialize it to init_val
//	Array3D(size_t w, size_t h, size_t d, T init_val) : Array3D(w, h, d) {
//		for (size_t i = 0; i < length; i++)
//			arr[i] = init_val;
//	}
//
//	// Copy CTOR
//	Array3D(Array3D& _arr) : Array3D(_arr.width, _arr.height, _arr.depth) {
//		for (size_t i = 0; i < length; i++)
//			arr[i] = _arr.arr[i];
//	}
//
//	// Move CTOR
//	Array3D(Array3D&& _arr) {
//		width = _arr.width;
//		height = _arr.height;
//		depth = _arr.depth;
//		length = _arr.length;
//
//		arr = _arr.arr;
//		_arr.arr = nullptr;
//	}
//
//	// Destructor
//	~Array3D() {
//		if (arr != nullptr) delete[] arr;
//	}
//
//	// Access to the array arr_name(x, y, z)
//	T& operator () (size_t x, size_t y, size_t z) {
//		return arr[height * depth * x + depth * y + z];
//	}
//
//	// Copy Assignment
//	Array3D<T>& operator = (const Array3D<T>& other) {
//		if (this == &other) return *this;
//
//		if (length != other.length) {
//			delete[] arr;
//			width = other.width;
//			height = other.height;
//			depth = other.depth;
//			length = other.length;
//			arr = new T[length];
//		}
//
//		for (size_t i = 0; i < length; i++)
//			arr[i] = other.arr[i];
//
//		return *this;
//	}
//
//	// Move Assignment
//	Array3D<T>& operator = (Array3D<T>&& other) {
//		if (this == &other) return *this;
//
//		width = other.width;
//		height = other.height;
//		depth = other.depth;
//		length = other.length;
//
//		arr = other.arr;
//		other.arr = nullptr;
//
//		return *this;
//	}
//};
//
//// Returns the probability that in the first bin will be exactly i items
//// Using dynamic programming the function fills the iter_x array with all the values
//// so all the probabilities are in there to use later
//// Where:
////			- input_length = How many items are in the input
////			- max_bin_size = The max size the bin can hold
////			- max_items_in_bin = The max items that can be in the bin
//double prob_X(Array3D<double> &iter_x, size_t input_length, size_t max_bin_size, size_t max_items_in_bin);
//// Returns the probability that k bins will be open
//// Using dynamic programming the function fills the iter_y array with all the values
//// Where:
////			- input_length = How many items are in the input
////			- open_bins = The max count of open bins that can be
////			- max_bin_size = The max items that can be in the bin
//double prob_Y(Array3D<double> &iter_y, Array3D<double> x, size_t input_length, size_t open_bins, size_t max_bin_size);
//// Returns the Expected Value (Mean) of open bins
//// Where:
////			- input_length = How many items are in the input
//double mean_E(Array3D<double> &iter_e, Array3D<double> y, size_t input_length);
//
//void main() {
//	size_t n, k;
//
//	bin_algorithm = FirstFit;
//
//	cout << "Enter n: " << endl;
//	cin >> n;
//
//	Array3D<double> x_prob_3d(n + 1, M + 1, M + 1, 0);
//	prob_X(x_prob_3d, n, M, M);
//
//	//for (size_t i = 0; i < n+1; i++)
//	//{
//	//	for (size_t j = 0; j < M+1; j++)
//	//	{
//	//		for (size_t k = 0; k < M+1; k++)
//	//		{
//	//			printf("x(%d,%d,%d) = %lf\n", i, j, k, x_prob_3d(i, j, k));
//	//		}
//	//	}
//	//}
//
//	Array3D<double> y_prob_3d(1, n + 1, n + 1, 0);
//	prob_Y(y_prob_3d, x_prob_3d, n, n, M);
//
//	Array3D<double> e_mean_3d(1, 1, n + 1, 0);
//	double expected_e = mean_E(e_mean_3d, y_prob_3d, n);
//
//	printf("Expected E(%d) = %lf\n", n, expected_e);
//
//	// this shit it used so the console dont exit, dont know why but other shit dont work
//	cin >> n;
//}
//
//double prob_X(Array3D<double> &iter_x, size_t input_length, size_t max_bin_size, size_t max_items_in_bin) {
//	// Calculates and insert the start conditions to the array
//	// X(1, m, 1) = (m/M), m=[0,M]
//	// X(n, m, 0) = (1 - (m/M))^(n), n=[0, N], m=[0,M]
//	for (size_t n = 0; n <= input_length; n++) {
//		for (size_t m = 0; m <= max_bin_size; m++) {
//			// Different calculation for different algorithms.. duhhh...
//			if (bin_algorithm == FirstFit)
//				iter_x(n, m, 0) = pow(1.0 - (((double)m) / ((double)M)), n);
//			else if (bin_algorithm == NextFit)
//				iter_x(n, m, 0) = 1.0 - (((double)m) / ((double)M));
//			if (n == 1)
//				iter_x(1, m, 1) = (((double)m) / ((double)M));
//		}
//	}
//
//	// Filling the array with all the dynamic programming calculations
//	// X(n, m, i) = 
//	//				SIGMA(j=1, m){ (1/M) * X(n-1, m-j, i-1) } + 
//	//				SIGMA(j=m+1, M){ (1/M) * X(n-1, m, i) }
//	for (size_t n = 1; n <= input_length; n++) {
//		for (size_t m = 1; m <= max_bin_size; m++) {
//			for (size_t i = 1; i <= max_items_in_bin; i++) {
//				// We did this calculations in the start conditions
//				// (i > n) ther cannot be i items in the bin if there only left n items in the input
//				if (i == 1 && n == 1 || i > n) continue;
//
//				printf("X(%d, %d, %d) =\n", n, m, i);
//
//				// The first SIGMA
//				for (size_t j = 1; j <= m; j++)
//					if (m - j >= 0) {
//						iter_x(n, m, i) += (iter_x(n - 1, m - j, i - 1) / ((double)M));
//						printf("\tx(%d,%d,%d) = %lf\n", n-1, m-j, i-1, iter_x(n-1, m-j, i-1));
//					}
//
//				printf("\tSecond Sigma:\n");
//				// The second SIGMA
//				if (bin_algorithm == FirstFit)
//					for (size_t j = m + 1; j <= M; j++) {
//						iter_x(n, m, i) += (iter_x(n - 1, m, i) / ((double)M));
//						printf("\tx(%d,%d,%d) = %lf\n", n - 1, m, i, iter_x(n - 1, m, i));
//					}
//				printf("\t\t Final x(%d,%d,%d) = %lf\n", n, m, i, iter_x(n, m, i));
//			}
//		}
//	}
//
//	return iter_x(input_length, max_bin_size, max_items_in_bin);
//}
//
//double prob_Y(Array3D<double> &iter_y, Array3D<double> x, size_t input_length, size_t open_bins, size_t max_bin_size) {
//	// Y(n, k) = SIGMA(j=1, M){ X(n, M,j) * Y(n-j, k-1) }
//	// With Starting condiotions:
//	//	- Y(0,0) = Y(1,1) = 1
//	//	- Y(1,k) = 0, k != 1
//	//	- Y(n,0) = 0, n > 0
//
//	// Setting the start conditions
//	iter_y(0, 0, 0) = 1;
//	iter_y(0, 1, 1) = 1;
//
//	// Y(1,k) for every k != 1
//	for (size_t k = 2; k <= open_bins; k++)
//		iter_y(0, 1, k) = 0;
//
//	// Y(n,0) for every n > 0
//	for (size_t n = 1; n <= input_length; n++)
//		iter_y(0, n, 0) = 0;
//
//	// Calculating the sigma itself
//	for (size_t n = 2; n <= input_length; n++)
//		for (size_t k = 0; k <= open_bins; k++)
//			for (size_t j = 1; j <= M; j++) {
//				if (n < j || k == 0) continue;
//				iter_y(0, n, k) += x(n, max_bin_size, j) * iter_y(0, n - j, k - 1);
//			}
//
//	return iter_y(0, input_length, open_bins);
//}
//
//double mean_E(Array3D<double> &iter_e, Array3D<double> y, size_t input_length) {
//
//	// Calculating the expected value (mean)
//	// E(n) = SIGMA(j=1, n){ Y(n, j) * j }
//	for (size_t j = 1; j <= input_length; j++)
//		iter_e(0, 0, j) = y(0, input_length, j) * j + iter_e(0, 0, j - 1);
//
//	return iter_e(0, 0, input_length);
//}