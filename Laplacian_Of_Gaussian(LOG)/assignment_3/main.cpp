#include <iostream>
#include <cstdint>
#include <cstdlib>

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
void LoG_func(uint8_t** image){
	FILE* fpOutputImage = 0;					// Output File
	FILE* fpLoGImage = 0;						// LoG file
	uint8_t** ppOutputImageBuffer = 0;			// output image buffer
	uint8_t** ppLoGImageBuffer = 0;				// LoG image buffer
	uint8_t** tempImage = 0;					// temp image
	uint8_t** g_Image = 0;						// gaussian image
	uint32_t IMG_HEIGHT = 512, IMG_WIDTH = 512;	// image height, image width
	int gaussian[3][3] = { { 1, 2, 1 }, { 2, 4, 2 }, { 1, 2, 1 } };		// gaussian mask
	int laplacian[5][5] = { { 1, 1, 1, 1, 1 }, { 1, 1, 1, 1, 1 }, { 1, 1, -24, 1, 1 }, { 1, 1, 1, 1, 1 }, { 1, 1, 1, 1, 1 } };	// laplacian mask
	double sum;									// sum
	int count;									// count
	int i, j, k, l, s;							// index

	// output file open
	fpOutputImage = fopen("result_image.raw", "wb");
	fpLoGImage = fopen("result_LoG_image.raw", "wb");

	// memory allocation
	ppLoGImageBuffer = memory_alloc2D(IMG_HEIGHT + 4, IMG_WIDTH + 4);
	g_Image = memory_alloc2D(IMG_HEIGHT + 4, IMG_WIDTH + 4);
	ppOutputImageBuffer = memory_alloc2D(IMG_HEIGHT, IMG_WIDTH);
	tempImage = memory_alloc2D(IMG_HEIGHT, IMG_WIDTH);
	
	// make 516x516 image
	for (i = 2; i < 514; i++){			// fill the middle image
		for (j = 2; j < 514; j++){
			ppLoGImageBuffer[i][j] = image[i - 2][j - 2];
		}
	}
	for (i = 0; i < 2; i++){			// fill the corner
		for (j = 0; j < 2; j++){
			for (k = 0; k < 2; k++){
				for (l = 0; l < 2; l++){
					ppLoGImageBuffer[i * 514 + k][j * 514 + l] = image[i * 510 + k][j * 510 + l];
				}
			}
		}
	}
	for (i = 0; i < 2; i++){		// fill the edge of image
		for (j = 2; j < 514; j++){
			for (k = 0; k < 2; k++){
				ppLoGImageBuffer[i * 514 + k][j] = image[i * 510 + k][j - 2];
				ppLoGImageBuffer[j][i * 514 + k] = image[j - 2][i * 510 + k];
			}
		}
	}
	for (i = 0; i < 516; i++){		// fill the gausian image
		for (j = 0; j < 516; j++)
			g_Image[i][j] = ppLoGImageBuffer[i][j];
	}

	/////////////////////////// start Laplancian of Gaussian//////////////////////////
	for (i = 0; i < 514; i++){			// Gaussian filter
		for (j = 0; j < 514; j++){
			sum = 0;	// initialize sum 
			for (k = 0; k < 3; k++){
				for (l = 0; l < 3; l++){
					sum += gaussian[k][l] * ppLoGImageBuffer[i + k][j + l];
				}
			}
			sum /= 16;
			if (sum < 0)					// sum is less than 0
				g_Image[i + 2][j + 2] = 0;
			else if (sum>255)				// sum is more than 255
				g_Image[i + 2][j + 2] = 255;
			else							// 0 <= sum <= 255
				g_Image[i + 2][j + 2] = (int)sum;
		}
	}
	for (i = 0; i < 512; i++){		// Laplacian filter
		for (j = 0; j < 512; j++){
			sum = 0;	// initialize sum
			for (k = 0; k < 5; k++){
				for (l = 0; l < 5; l++){
					sum += laplacian[k][l] * g_Image[i + k][j + l];
				}
			}
			if (sum < 0)					// sum is less than 0
				tempImage[i][j] = 0;
			else if (sum>255)				// sum is more than 255
				tempImage[i][j] = 255;
			else							// 0 <= sum <= 255
				tempImage[i][j] = (int)sum;
		}
	}
	// write the file
	fwrite(&tempImage[0][0], sizeof(uint8_t), IMG_HEIGHT*IMG_WIDTH, fpLoGImage);

	/////////////////// execute LoG + alpa//////////////////////
	// execute zero cossing
	for (i = 0; i < 512; i++){
		for (j = 0; j < 512; j++){
			if (tempImage[i][j] > 50) // compare pixel's value and threshold
				tempImage[i][j] = 255;
		}
	}
	// calculate 255 pixel
	for (s = 0; s < 2; s++){
		for (i = 0; i < 510; i++){
			for (j = 0; j < 510; j++){
				count = 0;
				// count 255 pixel
				for (k = i; k < i + 3; k++){
					for (l = j; l < j + 3; l++){
						if (tempImage[k][l] == 255)
							count++;
					}
				}
				// save pixel's value
				if (s == 0){
					if (count < 5)
						ppOutputImageBuffer[i + 1][j + 1] = 0;
					else
						ppOutputImageBuffer[i + 1][j + 1] = 255;
				}
				else{
					if (count < 6)
						ppOutputImageBuffer[i + 1][j + 1] = 0;
				}
			}
		}
	}

	// write the file
	fwrite(&ppOutputImageBuffer[0][0], sizeof(uint8_t), IMG_HEIGHT*IMG_WIDTH, fpOutputImage);

	// free memory and close file
	memory_free2D(ppLoGImageBuffer);
	memory_free2D(ppOutputImageBuffer);
	memory_free2D(tempImage);
	memory_free2D(g_Image);
	fclose(fpOutputImage);
	fclose(fpLoGImage);
}
int main(void){
	FILE* fpInputImage = 0;						// Input File
	uint8_t** ppInputImageBuffer = 0;			// input image buffer
	uint32_t IMG_HEIGHT = 512, IMG_WIDTH = 512;	// image height, image width

	// input file open
	fpInputImage = fopen("07_gLenna_512_512.raw", "rb");

	// memory allocation
	ppInputImageBuffer = memory_alloc2D(IMG_HEIGHT, IMG_WIDTH);

	// input file read to memory from the file
	fread(&ppInputImageBuffer[0][0], sizeof(uint8_t), IMG_WIDTH*IMG_HEIGHT, fpInputImage);

	LoG_func(ppInputImageBuffer);

	// free memory and close file
	memory_free2D(ppInputImageBuffer);
	fclose(fpInputImage);

	return 0;
}