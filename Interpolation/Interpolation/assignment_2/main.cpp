#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <algorithm>
#include <fstream>
#include <cmath>

using namespace std;

uint8_t** memory_alloc2D(uint32_t height, uint32_t width){
	uint8_t** ppMem2D = 0;
	uint32_t j;
	
	// array of pointer
	ppMem2D = (uint8_t**)calloc(sizeof(uint8_t*), height);
	if (ppMem2D == 0)
		return 0;

	*ppMem2D = (uint8_t*)calloc(sizeof(uint8_t), width*height);
	if ((*ppMem2D) == 0){
		// free the memory of array of pointer
		free(ppMem2D);
		return 0;
	}

	// allocate memory
	for (j = 1; j < height; j++)
		ppMem2D[j] = ppMem2D[j - 1] + width;

	return ppMem2D;
}
int memory_free2D(uint8_t** ppMemAllocated){
	if (ppMemAllocated == 0)
		return -1;
	
	// free memory
	free(ppMemAllocated[0]);
	free(ppMemAllocated);

	return 0;
}
void PSNR_function(double rms){
	double psnr;	// PSNR

	// calculate PSNR
	psnr = 20 * log10(255 / rms);

	cout << "PSNR : " << psnr << endl;
}
void RMS_function(uint8_t** image, uint8_t** ppInputImageBuffer){
	double rms = 0;	// RMS
	int i, j;		// index;
	
	// calculate RMS
	for (i = 0; i < 512; i++){
		for (j = 0; j < 512; j++){
			rms += (ppInputImageBuffer[i][j] - image[i][j])*(ppInputImageBuffer[i][j] - image[i][j]);
		}
	}
	rms = sqrt(rms / (512 * 512));

	cout << "RMS : " << rms << "\t";
	PSNR_function(rms);	// calculate PSNR
}
void fill_pixel(uint8_t** image, int index, int direction, int how){
	int left = 2, right = 6;	// left, right
	int l;						// left
	float L0, L1, L2, L3;		// Lagrange values
	float L_f;					// Lagrange function's values
	int sum = 0;				// sum;
	int i, j;					// index

	// if direction is 0, the line is filled.
	// else, the row is filled.

	// if how is 0, the bilinear interpolation is executed.
	// els if how is 1, the Lagrand interpolation is executed.
	// else, the six_tab interpolation is execued.

	if (how == 0){	// bilinear interpolation
		while (left < 510){
			for (j = 0; j < 3; j++){	// fill the values which is between left and right
				if (direction == 0)	// fill the line
					image[index][left + j + 1] = (image[index][left] * (3 - j) + image[index][right] * (j + 1)) / 4;
				else				// fill the row
					image[left + j + 1][index] = (image[left][index] * (3 - j) + image[right][index] * (j + 1)) / 4;
			}
			// move left and right
			left = right;
			right += 4;
		}
	}
	else if (how == 1){	// lagrange interpolation
		while (left < 500){
			l = left;	// initailize l
			for (i = 0; i < 3; i++, l += 4){
				for (j = l + 1; j < l + 4; j++){
					// calculate Lagrange
					L0 = ((float)(j - left - 4)*(j - left - 8)*(j - left - 12))/(-384);
					L1 = ((float)(j - left)*(j - left - 8)*(j - left - 12))/128;
					L2 = ((float)(j - left)*(j - left - 4)*(j - left - 12))/(-128);
					L3 = ((float)(j - left)*(j - left - 4)*(j - left - 8)) / 384;
					if (direction == 0){	// fill the line
						L_f = L0*image[index][left] + L1*image[index][left + 4] + L2*image[index][left + 8] + L3*image[index][left + 12];
						if (L_f < 0) L_f = 0;			// Lagrange function's value is less than 0
						else if (L_f>255) L_f = 255;	// Lagrange function's value is more than 255
						image[index][j] = (uint8_t)L_f;
					}
					else{						// fille the row
						L_f = L0*image[left][index] + L1*image[left + 4][index] + L2*image[left + 8][index] + L3*image[left + 12][index];
						if (L_f < 0) L_f = 0;			// Lagrange function's value is less than 0
						else if (L_f>255) L_f = 255;	// Lagrange function's value is more than 255
						image[j][index] = (uint8_t)L_f;
					}
				}
			}
			left += 4; // move to left plus 4
		}
	}
	else{				// six_tab interpolation
		right = 510;	// initialize right
		for (i = left + 2; i < 510; i += 4){
			if (i == left + 2){	// the standard point is 4 
				if (direction == 0) // fill the line
					sum = image[index][left] * 16 + image[index][left + 4] * 20 - image[index][left + 8] * 5 + image[index][left + 12];
				else				// fill the row
					sum = image[left][index] * 16 + image[left + 4][index] * 20 - image[left + 8][index] * 5 + image[left + 12][index];
			}
			else if (i == left + 6){	// the standard point is 8
				if (direction == 0)		// fill the line
					sum = image[index][left] * (-4) + image[index][left + 4] * 20 + image[index][left + 8] * 20 - image[index][left + 12] * 5 + image[index][left + 16];
				else				// fill the row
					sum = image[left][index] * (-4) + image[left + 4][index] * 20 + image[left + 8][index] * 20 - image[left + 12][index] * 5 + image[left + 16][index];
			}
			else if (i == right - 2){	// the standard point is 508
				if (direction == 0)		// fill the line
					sum = image[index][right] * 16 + image[index][right - 4] * 20 - image[index][right - 8] * 5 + image[index][right - 12];
				else				// fill the row
					sum = image[right][index] * 16 + image[right - 4][index] * 20 - image[right - 8][index] * 5 + image[right - 12][index];
			}
			else if (i == right - 6){	// the standard point is 504
				if (direction == 0)		// fill the line
					sum = image[index][right] * (-4) + image[index][right - 4] * 20 + image[index][right - 8] * 20 - image[index][right - 12] * 5 + image[index][right - 16];
				else				// fill the row
					sum = image[right][index] * (-4) + image[right - 4][index] * 20 + image[right - 8][index] * 20 - image[right - 12][index] * 5 + image[right - 16][index];
			}
			else{						// 8 < the standard point <504
				if (direction == 0)		// fill the line
					sum = image[index][i - 10] - image[index][i - 6] * 5 + image[index][i - 2] * 20 + image[index][i + 2] * 20 - image[index][i + 6] * 5 + image[index][i + 10];
				else					// fill the row
					sum = image[i - 10][index] - image[i - 6][index] * 5 + image[i - 2][index] * 20 + image[i + 2][index] * 20 - image[i + 6][index] * 5 + image[i + 10][index];
			}

			// insert values
			if (direction == 0){   // fill the line
				image[index][i] = sum / 32;
				image[index][i - 1] = (image[index][i] + image[index][i - 2]) / 2;
				image[index][i + 1] = (image[index][i] + image[index][i + 2]) / 2;
			}
			else{					// fill the row
				image[i][index] = sum / 32;
				image[i - 1][index] = (image[i][index] + image[i - 2][index]) / 2;
				image[i + 1][index] = (image[i][index] + image[i + 2][index]) / 2;
			}
		}
	}
}
void padding(uint8_t** image){
	int i, j;		// index

	// fill the the external parts of 512x512 images
	for (i = 0; i < 2; i++){
		for (j = 2; j < 510; j++){
			image[i][j] = image[2][j];
			image[i + 510][j] = image[509][j];
			image[j][i] = image[j][2];
			image[j][i + 510] = image[j][509];
		}
	}
}
void Bilinear(uint8_t** downsamplingImage, uint8_t** ppInputImageBuffer){
	FILE* fpBilinearImage = 0;					// File
	uint8_t** bilinearImage = 0;				// bilinear image
	uint32_t IMG_HEIGHT = 512, IMG_WIDTH = 512;	// image height, image width
	int i, j, l, r;								// index
	
	// output fileopen
	fpBilinearImage = fopen("bilinear.raw", "wb");

	// memory allocation
	bilinearImage = memory_alloc2D(IMG_HEIGHT, IMG_WIDTH);

	// fill the corner of image
	for (i = 0; i < 2; i++){
		for (j = 0; j < 2; j++){
			for (l = 0; l < 2; l++){
				for (r = 0; r < 2; r++){
					bilinearImage[i * 510 + l][j * 510 + r] = downsamplingImage[i * 127][j * 127];
				}
			}
		}
	}
	
	// fill the inside image
	for (i = 0; i < 128; i++){
		for (j = 0; j < 128; j++){
			bilinearImage[(i * 4) + 2][(j * 4) + 2] = downsamplingImage[i][j];
		}
		// fill pixel (line)
		fill_pixel(bilinearImage, (i * 4) + 2, 0, 0);
	}
	// fill pixel (row)
	for (i = 2; i < 511; i++)
		fill_pixel(bilinearImage, i, 1, 0);

	// padding
	padding(bilinearImage);

	// write the file
	fwrite(&bilinearImage[0][0], sizeof(uint8_t), IMG_HEIGHT*IMG_WIDTH, fpBilinearImage);

	cout << "<Bilinear interpolation>" << endl;
	RMS_function(bilinearImage, ppInputImageBuffer);	// go to RMS_function

	// free memory and close file
	memory_free2D(bilinearImage);
	fclose(fpBilinearImage);
}
void Lagrange(uint8_t** downsamplingImage, uint8_t** ppInputImageBuffer){
	FILE* fpLagrangeImage = 0;					// File
	uint8_t** lagrangeImage = 0;				// lagrange image
	uint32_t IMG_HEIGHT = 512, IMG_WIDTH = 512;	// image height, image width
	int i, j, l, r;								// index
	
	// output fileopen
	fpLagrangeImage = fopen("lagrange.raw", "wb");

	// memory allocation
	lagrangeImage = memory_alloc2D(IMG_HEIGHT, IMG_WIDTH);

	// fill the corner of image
	for (i = 0; i < 2; i++){
		for (j = 0; j < 2; j++){
			for (l = 0; l < 2; l++){
				for (r = 0; r < 2; r++){
					lagrangeImage[i * 510 + l][j * 510 + r] = downsamplingImage[i * 127][j * 127];
				}
			}
		}
	}

	// fill the inside image
	for (i = 0; i < 128; i++){
		for (j = 0; j < 128; j++){
			lagrangeImage[(i * 4) + 2][(j * 4) + 2] = downsamplingImage[i][j];
		}
		// fill pixel (line)
		fill_pixel(lagrangeImage, (i * 4) + 2, 0, 1);
	}
	// fill pixel (row)
	for (i = 2; i < 511; i++)
		fill_pixel(lagrangeImage, i, 1, 1);

	// padding
	padding(lagrangeImage);

	// write the file
	fwrite(&lagrangeImage[0][0], sizeof(uint8_t), IMG_HEIGHT*IMG_WIDTH, fpLagrangeImage);
	
	cout << "<Lagrange interpolation>" << endl;
	RMS_function(lagrangeImage, ppInputImageBuffer);	// go to RMS_function

	// free memory and close file
	memory_free2D(lagrangeImage);
	fclose(fpLagrangeImage);
}
void Six_tab(uint8_t** downsamplingImage, uint8_t** ppInputImageBuffer){
	FILE* fpSix_tabImage = 0;					// File
	uint8_t** six_tabImage = 0;					// six_tab image
	uint32_t IMG_HEIGHT = 512, IMG_WIDTH = 512;	// image height, image width
	int i, j, l, r;								// index

	// output fileopen
	fpSix_tabImage = fopen("six_tab.raw", "wb");

	// memory allocation
	six_tabImage = memory_alloc2D(IMG_HEIGHT, IMG_WIDTH);

	// fill the corner of image
	for (i = 0; i < 2; i++){
		for (j = 0; j < 2; j++){
			for (l = 0; l < 2; l++){
				for (r = 0; r < 2; r++){
					six_tabImage[i * 510 + l][j * 510 + r] = downsamplingImage[i * 127][j * 127];
				}
			}
		}
	}

	// fill the inside image
	for (i = 0; i < 128; i++){
		for (j = 0; j < 128; j++){
			six_tabImage[(i * 4) + 2][(j * 4) + 2] = downsamplingImage[i][j];
		}
		// fill pixel (line)
		fill_pixel(six_tabImage, (i * 4) + 2, 0, 2);
	}
	// fill pixel (row)
	for (i = 2; i < 512; i += 2)
		fill_pixel(six_tabImage, i, 1, 2);
	
	// fill pixel (line)
	for (i = 4; i < 512;i+=4)
		fill_pixel(six_tabImage, i, 0, 2);

	// fill pixel (diagonal line)
	for (i = 2; i < 510; i += 4){
		for (j = 2; j < 510; j += 4){
			six_tabImage[i + 1][j + 1] = (six_tabImage[i][j + 2] + six_tabImage[i + 2][j]) / 2;
			six_tabImage[i + 1][j + 3] = (six_tabImage[i][j + 2] + six_tabImage[i + 2][j + 4]) / 2;
			six_tabImage[i + 3][j + 1] = (six_tabImage[i + 2][j] + six_tabImage[i + 4][j + 2]) / 2;
			six_tabImage[i + 3][j + 3] = (six_tabImage[i + 2][j + 4] + six_tabImage[i + 4][j + 2]) / 2;
		}
	}
	
	// padding
	padding(six_tabImage);

	// write the file
	fwrite(&six_tabImage[0][0], sizeof(uint8_t), IMG_HEIGHT*IMG_WIDTH, fpSix_tabImage);

	cout << "<Six_tab interpolation>" << endl;
	RMS_function(six_tabImage, ppInputImageBuffer);	// go to RMS_function

	// free memory and close file
	memory_free2D(six_tabImage);
	fclose(fpSix_tabImage);
}
void DownSampling(uint8_t** ppInputImageBuffer){
	FILE* fpDownsamplingImage = 0;					// File
	uint8_t** downsamplingImage = 0;				// downsampling image
	uint8_t arr[16] = {};							// array
	uint32_t IMG_HEIGHT = 128, IMG_WIDTH = 128;		// image height, image width
	int row = 0, line = 0;							// row, line
	int i, j, l, r, index;							// index
	
	// output fileopen
	fpDownsamplingImage = fopen("07_gLenna_128_128.raw", "wb");

	// memory allocation
	downsamplingImage = memory_alloc2D(IMG_HEIGHT, IMG_WIDTH);
	
	// execute downsampling
	for (i = 0; i < 512; i += 4, line++){
		row = 0;		// initialize row
		for (j = 0; j < 512; j += 4, row++){
			index = 0;				// initialize index
			for (l = i; l < 4 + i; l++){
				for (r = j; r < 4 + j; r++){
					arr[index++] = ppInputImageBuffer[l][r];
				}
			}
			sort(&arr[0], &arr[16]);					// sort array
			downsamplingImage[line][row] = arr[15 / 2];	// allocate the value
		}
	}

	// write the file
	fwrite(&downsamplingImage[0][0], sizeof(uint8_t), IMG_HEIGHT*IMG_WIDTH, fpDownsamplingImage);

	Bilinear(downsamplingImage, ppInputImageBuffer);	// excute bilinear interpolation
	Lagrange(downsamplingImage, ppInputImageBuffer);	// excute lagrange interpolation
	Six_tab(downsamplingImage, ppInputImageBuffer);		// excute six_tab interpolation

	memory_free2D(downsamplingImage);
	fclose(fpDownsamplingImage);
}
int main(void){
	FILE* fpInputImage = 0;						// File
	uint8_t** ppInputImageBuffer = 0;			// input image buffer
	uint32_t IMG_HEIGHT = 512, IMG_WIDTH = 512;	// image height, image width

	// input file open
	fpInputImage = fopen("07_gLenna_512_512.raw", "rb");

	// memory allocation
	ppInputImageBuffer = memory_alloc2D(IMG_HEIGHT, IMG_WIDTH);

	// input file read to memory from the file
	fread(&ppInputImageBuffer[0][0], sizeof(uint8_t), IMG_WIDTH*IMG_HEIGHT, fpInputImage);

	// get downsampling Image
	DownSampling(ppInputImageBuffer);

	memory_free2D(ppInputImageBuffer);
	fclose(fpInputImage);

	return 0;
}