#include <stdio.h>
#include <stdlib.h>

void printMap(unsigned char* A, int w, int h) {
	int x, y;
	for (y = 0; y < h; y++) {
		for (x = 0; x < w; x++) {
			printf("%hhu ", A[y*w + x]);
		}
		printf("\n");
	}
}

void neighbours(int x, int y, int w, int h, int* coords) {
	int xl = (x-1) < 0 ? w-1 : x-1;
	int xh = (x+1) >= w ? 0 : x+1;
	int yl = (y-1) < 0 ? h-1 : y-1;
	int yh = (y+1) >= h ? 0 : y+1;
	coords[0]  = xl; coords[1]  = yl;
	coords[2]  = xl; coords[3]  = y;
	coords[4]  = xl; coords[5]  = yh;
	coords[6]  = xh; coords[7]  = yl;
	coords[8]  = xh; coords[9]  = y;
	coords[10] = xh; coords[11] = yh;
	coords[12] = x;  coords[13] = yl;
	coords[14] = x;  coords[15] = yh;
}

int main() {
	int x, y, i, j;
	int w, h;
	scanf("%d %d", &w, &h);
	unsigned char *A = (unsigned char*) malloc(w*h*sizeof(unsigned char));
	unsigned char *B = (unsigned char*) malloc(w*h*sizeof(unsigned char));
	int *coords = (int*) malloc(8*2*sizeof(int));

	// Read shit
	for (y = 0; y < h; y++) {
		for (x = 0; x < w; x++) {
			scanf("%hhu", A + y*w + x);
		}
	}

	for (j = 0; j < 30; j++) {
		for (y = 0; y < h; y++) {
			for (x = 0; x < w; x++) {
				neighbours(x, y, w, h, coords);
				int count = 0;
				for (i = 0; i < 16; i += 2) {
					count += A[coords[i+1]*w + coords[i]];
				}

				if (A[y*w+x] == 1) {		// Is alive
					if (count == 2 || count == 3) {
						B[y*w+x] = 1;
					} else {
						B[y*w+x] = 0;
					}
				} else {					// Is dead
					if (count == 3) {
						B[y*w+x] = 1;
					} else {
						B[y*w+x] = 0;
					}
				}
			}
		}
		printMap(B, w, h);
		unsigned char *tmp = A;
		A = B;
		B = tmp;
	}

	


	free(A);
	free(B);
	free(coords);
}