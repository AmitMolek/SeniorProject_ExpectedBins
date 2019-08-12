#pragma once

#include <fstream>

using namespace std;

class CSVWriter {

public:
	static void create_csv_from_3d_array(string file_name, double arr[], size_t length);

};

