#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define NC 78		// Neighbor coordinates
#define MAXITER 30

void generateRandom(int w, int h, int d, double p, unsigned char *A) {
	int x, y, z;
	for (z = 0; z < d; z++) {
		for (y = 0; y < h; y++) {
			for (x = 0; x < w; x++) {
				double r = (double) rand() / RAND_MAX;
				if (r <= p) {
					A[z*w*h + y*w + x] = 1;
				} else {
					A[z*w*h + y*w + x] = 0;
				}
			}
		}
	}
}

void printMap(unsigned char* A, int w, int h, int d) {
	int x, y, z;
	//printf("%d %d %d\n", w, h, d);
	for (z = 0; z < d; z++) {
		for (y = 0; y < h; y++) {
			for (x = 0; x < w; x++) {
				printf("%hhu ", A[z*w*h + y*w + x]);
			}
			printf("\n");
		}
		//printf("\n");
	}
}

void neighbours(int x, int y, int z, int w, int h, int d, int* coords) {
	int xl = (x-1) < 0 ? w-1 : x-1;
	int xh = (x+1) >= w ? 0 : x+1;
	int yl = (y-1) < 0 ? h-1 : y-1;
	int yh = (y+1) >= h ? 0 : y+1;
	int zl = (z-1) < 0 ? d-1 : z-1;
	int zh = (z+1) >= d ? 0 : z+1;
	
	int xs[] = {xl, x, xh};
	int ys[] = {yl, y, yh};
	int zs[] = {zl, z, zh};
	int i, j, k, l = 0;
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			for (k = 0; k < 3; k++) {
				if (i == 1 && j == 1 && k == 1) {
					continue;
				}
				coords[l] = xs[i];
				coords[l+1] = ys[j];
				coords[l+2] = zs[k];
				l += 3;
			}
		}
	}
}

int countAlive(int x, int y, int z, int w, int h, int d, unsigned char *A, int *coords) {
	int count = 0, i;
	neighbours(x, y, z, w, h, d, coords);
	for (i = 0; i < NC; i += 3) {
		count += A[coords[i+2]*w*h + coords[i+1]*w + coords[i]];
	}
	return count;
}

int main() {
	int x, y, z, i, j, it;
	int w, h, d;
	
	scanf("%d %d %d", &w, &h, &d);
	//w = 100; h = 100; d = 100;
	unsigned char *A = (unsigned char*) malloc(w*h*d*sizeof(unsigned char));
	unsigned char *B = (unsigned char*) malloc(w*h*d*sizeof(unsigned char));
	int *coords = (int*) malloc(NC*sizeof(int));

	srand(time(NULL));
	//generateRandom(w, h, d, 0.3, A);

	// Read shit
	
	for (z = 0; z < d; z++) {
		for (y = 0; y < h; y++) {
			for (x = 0; x < w; x++) {
				scanf("%hhu", A + z*w*h + y*w + x);
			}
		}
	}
	

	printMap(A, w, h, d);

	for (it = 0; it < MAXITER; it++) {
		//fprintf(stderr, "%d\n", it);
		
		for (z = 0; z < d; z++) {
			for (y = 0; y < h; y++) {
				for (x = 0; x < w; x++) {
		
					int count = countAlive(x, y, z, w, h, d, A, coords);
					int idx = z*w*h + y*w + x;

					if (A[idx] == 1) {		// Is alive
						if (count == 4 || count == 5) {
							B[idx] = 1;
						} else {
							B[idx] = 0;
						}
					} else {					// Is dead
						if (count == 4) {
							B[idx] = 1;
						} else {
							B[idx] = 0;
						}
					}
				}
			}
		}

		printMap(B, w, h, d);
		unsigned char *tmp = A;
		A = B;
		B = tmp;
	}

	


	free(A);
	free(B);
	free(coords);
}