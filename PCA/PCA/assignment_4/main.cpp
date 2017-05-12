#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include "svdcmp.c"

using namespace std;

double** memory_alloc(int height, int width){
	double** m;

	// memory allocation
	m = new double*[height];
	for (int i = 0; i < height; i++)
		m[i] = new double[width];

	return m;
}
int memory_free(double** ppMemAllocated){
	if (ppMemAllocated == 0)
		return -1;

	// free memory
	free(ppMemAllocated[0]);
	free(ppMemAllocated);

	return 0;
}
void bmp_to_raw(char* f_name, double** image, int k){
	ifstream file(f_name, ios::binary);		// file open
	double* pixel, **f_image;				// pixel, file image
	int i, index = 0;					// index

	// memory allocation
	f_image = memory_alloc(10304, 3);
	pixel = new double[10304];

	// file doesn't open
	if (!file.is_open())
		return;

	// read header
	for (i = 0; i < 54; i++)
		file.get();

	// read pixel
	for (i = 0; i < 10304; i++){
		f_image[i][0] = file.get();
		f_image[i][1] = file.get();
		f_image[i][2] = file.get();
		pixel[i] = f_image[i][0] * 0.114 + f_image[i][1] * 0.587 + f_image[i][2] * 0.299;
	}

	image[k] = pixel;	// input image in matrix

	file.close();				// file close
	memory_free(f_image);		// free memory
}
void PCA(double** t_image, double* a_vector, double** t_vector, double** t_matrix, double** c_matrix, double* e_value, double** e_vector, double** V){
	double** v;							// 10304x10304 eigenvector
	double m = 0;						// average vector
	double mul = 0;						// for multiply
	double covariance = 0;				// for covariance matrix
	int i, j, l;						// index

	// memory allocation
	v = memory_alloc(10304, 48);

	///////////////////////// average vector /////////////////////////
	for (i = 0; i < 10304; i++){
		m = 0;
		for (j = 0; j < 48; j++)
			m += t_image[j][i];
		a_vector[i] = m / 48;
	}

	///////////////////////// training vector(10304x48) /////////////////////////
	for (i = 0; i < 10304; i++){
		for (j = 0; j < 48; j++)
			t_vector[i][j] = t_image[j][i] - a_vector[i];
	}

	//////////////////////// matrix transpose (symmetric matrix) (48x10304) ////////////////////////////////////
	for (i = 0; i < 48; i++){
		for (j = 0; j < 10304; j++){
			t_matrix[i][j] = t_vector[j][i];
		}
	}

	//////////////////////// covariance matrix(48x48) ///////////////////////////////
	for (i = 0; i < 48; i++){
		for (j = 0; j < 48; j++){
			covariance = 0;
			for (l = 0; l < 10304; l++)
				covariance += t_vector[l][j] * t_matrix[j][l];
			c_matrix[i][j] = covariance;
		}
	}

	/////////////////////////// eigenvalues & eigenvector(numerical recipes) /////////////////////////////////
	svdcmp(c_matrix, 48, 48, e_value, e_vector);

	//////////////////////// 10304x48 eigenvector ////////////////////////////////
	for (i = 0; i < 10304; i++){
		for (j = 0; j < 48; j++){
			mul = 0;
			for (l = 0; l < 48; l++)
				mul += t_vector[i][l] * e_vector[l][j];
			v[i][j] = mul;
		}
	}

	////////////////////////// 10304X10 eigenvector ////////////////////////////////////
	for (i = 0; i < 10304; i++){
		for (j = 0; j < 10; j++)
			V[i][j] = v[i][j];
	}
}
void E_distance(double** y_t_image, double** y_t_vector, double** V, double** t_V, double* a_vector,double** t_vector, double** r_vector, double** y_r_vector, double** e_distance){
	double mul;			// mul
	int i, j, l;		// index

	///////////////////////// y's training vector(10304x12) /////////////////////////
	for (i = 0; i < 10304; i++){
		for (j = 0; j < 12; j++)
			y_t_vector[i][j] = y_t_image[j][i] - a_vector[i];
	}

	////////////////////////// V transpose (10X10304)	///////////////////////////////////
	for (i = 0; i < 10; i++){
		for (j = 0; j < 10304; j++)
			t_V[i][j] = V[j][i];
	}

	///////////////////////// x's representative vector(10x48) //////////////////////////
	for (i = 0; i < 48; i++){
		for (j = 0; j < 10; j++){
			mul = 0;
			for (l = 0; l < 10304; l++)
				mul += t_V[j][l] * t_vector[l][i];
			r_vector[j][i] = mul;
		}
	}

	//////////////////////// y's representative vector(10x12) //////////////////////////
	for (i = 0; i < 12; i++){
		for (j = 0; j < 10; j++){
			mul = 0;
			for (l = 0; l < 10304; l++)
				mul += t_V[j][l] * y_t_vector[l][i];
			y_r_vector[j][i] = mul;
		}
	}

	/////////////////////// calculate euclidean distance(12x48) ////////////////////////////
	for (i = 0; i < 12; i++){
		for (j = 0; j < 48; j++){
			mul = 0;
			for (l = 0; l < 10; l++)
				mul += pow(r_vector[l][j] - y_r_vector[l][i], 2.0);
			e_distance[i][j] = sqrt(mul);
		}
	}
}
void find_picture(double** e_distance){
	double compare;	// compare
	int p_num;		// picture number
	int i, j;		// index

	///////////////////// find same picture ///////////////////
	for (i = 0; i < 12; i++){
		p_num = 0;
		compare = e_distance[i][0];
		for (j = 0; j < 48; j++){
			if (compare>e_distance[i][j]){
				compare = e_distance[i][j];
				p_num = j;
			}
		}
		if (p_num < 8)
			cout << i + 1 << "th picture is simular to set_1" << endl;
		else if (p_num<16)
			cout << i + 1 << "th picture is simular to set_2" << endl;
		else if (p_num<24)
			cout << i + 1 << "th picture is simular to set_3" << endl;
		else if (p_num<32)
			cout << i + 1 << "th picture is simular to set_4" << endl;
		else if (p_num<40)
			cout << i + 1 << "th picture is simular to set_5" << endl;
		else
			cout << i + 1 << "th picture is simular to set_6" << endl;
	}
}
int main(void){
	double** t_image, **y_t_image;		// training image
	double** t_vector, **y_t_vector;	// training vector
	double* a_vector;					// average vector
	double** c_matrix;					// covariance matrix
	double** t_matrix;					// matrix transpose
	double** e_vector;					// 48x48 eigenvector
	double* e_value;					// eigenvalue
	double** V;							// 48x10 eigenvector
	double** t_V;						// V transepose
	double** r_vector, **y_r_vector;	// representative vector
	double** e_distance;				// euclidean distance
	char path[1024] = {};				// path
	int i, j, k = 0, l = 0;				// index

	// memory allocation
	t_image = memory_alloc(48, 10304);
	t_vector = memory_alloc(10304, 48);
	t_matrix = memory_alloc(48, 10304);
	c_matrix = memory_alloc(48, 48);
	e_vector = memory_alloc(48, 48);
	V = memory_alloc(10304, 10);
	t_V = memory_alloc(10, 10304);
	y_t_image = memory_alloc(12, 10304);
	y_t_vector = memory_alloc(10304, 12);
	r_vector = memory_alloc(10, 48);
	y_r_vector = memory_alloc(10, 12);
	e_distance = memory_alloc(12, 48);
	a_vector = new double[10304];
	e_value = new double[48];

	//////////////////////////////// training image ///////////////////////////////////////
	// convert bmp to raw
	for (i = 1; i < 7; i++){								// forder
		for (j = 1; j < 9; j++){							// file
			sprintf(path, "Training/set_%d/%d.bmp", i, j);
			bmp_to_raw(path, t_image, k++);
		}
	}
	
	// PCA
	PCA(t_image, a_vector, t_vector, t_matrix, c_matrix, e_value, e_vector, V);

	////////////////////////////// test image //////////////////////////////////
	// convert bmp to raw
	for (i = 1; i < 13; i++){
		sprintf(path, "Test/%d.bmp", i, j);
		bmp_to_raw(path, y_t_image, l++);
	}

	// Euclidean distance
	E_distance(y_t_image, y_t_vector, V, t_V, a_vector, t_vector, r_vector, y_r_vector, e_distance);

	// find same picture
	find_picture(e_distance);

	// free memory
	memory_free(t_image);
	memory_free(y_t_image);
	memory_free(t_vector);
	memory_free(y_t_vector);
	memory_free(c_matrix);
	memory_free(t_matrix);
	memory_free(e_vector);
	memory_free(V);
	memory_free(t_V);
	memory_free(r_vector);
	memory_free(y_r_vector);
	memory_free(e_distance);
	free(a_vector);
	free(e_value);

	return 0;
}