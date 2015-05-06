/* Use MPI */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "updateFunctions.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))

#define iterations 5 /* number of iterations */

/* Number of processors in each dimension */
#define Px 2
#define Py 2
#define Pz 2

/* Periodicity of boundaries, should always be 1 */
#define xBC 1
#define yBC 1
#define zBC 1

/* Size of each dimension */
#define xN 10
#define yN 10
#define zN 10

/**/
#define X 0
#define Y 1
#define Z 2

/* Check function. */
void check(int rc) {
	if (rc != MPI_SUCCESS) {
		exit(-1);
	}
}

void copyData(int Ix, int Iy, int Iz, unsigned char data[Ix+2][Iy+2][Iz+2],
				   unsigned char newData[Ix][Iy][Iz], 
				  int t, unsigned char totalProcessorResults[Ix][Iy][Iz][iterations+1]) {
	int x, y, z;
	for (x=0;x<Ix;x++) {
		for (y=0;y<Iy;y++) {
			for (z=0;z<Iz;z++) {
				data[x+1][y+1][z+1] = newData[x][y][z];
				totalProcessorResults[x][y][z][t] = newData[x][y][z];
			}
		}
	}
}





int main(int argc, char *argv[]){
	//printf("start\n");
	/* Constants we need */
	int rank, size, i;
	int Nx, Lx, Rx, Ix;
	int Ny, Ly, Ry, Iy;
	int Nz, Lz, Rz, Iz;
	int Xforward, Xbackward, Yforward, Ybackward, Zforward, Zbackward; 
	int  FORWARD = 1, BACKWARDS = -1;
    int dim[3], period[3], reorder;
    int coord[3], id;
    int doMPI = 1;
    MPI_Comm TORUS_COMM;

    /* initialize MPI and create virtual 2d torus topology*/

    check ( MPI_Init(&argc, &argv) );
    check ( MPI_Comm_rank(MPI_COMM_WORLD, &rank) );
    check ( MPI_Comm_size(MPI_COMM_WORLD, &size) );

    dim[0] = Px; dim[1] = Px; dim[2] = Px;
	period[0] = xBC; period[1] = yBC; period[2] = zBC;
	reorder = 0;
	//printf("create cart\n");
	check ( MPI_Cart_create(MPI_COMM_WORLD, 3, dim, period, reorder, &TORUS_COMM) );
	check ( MPI_Cart_coords(TORUS_COMM, rank, 3, coord) );
	
	/* 3d linear data distribution */
	Nx = xN;
	Lx = Nx / Px;
	Rx = Nx % Px;
	Ix = ( Nx + Px - coord[0]-1 ) / Px;

	Ny = yN;
	Ly = Ny / Py;
	Ry = Ny % Py;
	Iy = ( Ny + Py - coord[1]-1 ) / Py;

	Nz = zN;
	Lz = Nz / Pz;
	Rz = Nz % Pz;
	Iz = ( Nz + Pz - coord[2]-1 ) / Pz;

	//printf("declerations\n");
	/* Arrays */

	unsigned char totalProcessorResults[Ix][Iy][Iz][iterations+1];

	unsigned char data[Ix+2][Iy+2][Iz+2];
	memset(data, 0, sizeof data);
	unsigned char newData[Ix][Iy][Iz];
	memset(newData, 0, sizeof newData);

	/* Receive buffers */
	/* Sides */
	unsigned char XforwardGhosts[Iy][Iz];
	memset(XforwardGhosts, 0, sizeof XforwardGhosts);
	unsigned char XbackwardGhosts[Iy][Iz];
	memset(XbackwardGhosts, 0, sizeof XbackwardGhosts);

	unsigned char YforwardGhosts[Iz][Ix];
	memset(YforwardGhosts, 0, sizeof YforwardGhosts);
	unsigned char YbackwardGhosts[Iz][Ix];
	memset(YbackwardGhosts, 0, sizeof YbackwardGhosts);

	unsigned char ZforwardGhosts[Ix][Iy];
	memset(ZforwardGhosts, 0, sizeof ZforwardGhosts);
	unsigned char ZbackwardGhosts[Ix][Iy];
	memset(ZbackwardGhosts, 0, sizeof ZbackwardGhosts);

	/* Edges */
	unsigned char XYGhosts[4][Iz];
	memset(XYGhosts, 0, sizeof XYGhosts);
	unsigned char YZGhosts[4][Ix];
	memset(YZGhosts, 0, sizeof YZGhosts);
	unsigned char ZXGhosts[4][Iy];
	memset(ZXGhosts, 0, sizeof ZXGhosts);

	/* corners are enumerated in binary as their coordinates [z, y, x]  */
	unsigned char cornerGhosts[8] = {0};
	memset(cornerGhosts, 0, sizeof cornerGhosts);

	/* Send buffers */
	/* Sides */
	unsigned char XforwardBuffers[Iy][Iz];
	memset(XforwardBuffers, 0, sizeof XforwardBuffers);
	unsigned char XbackwardBuffers[Iy][Iz];
	memset(XbackwardBuffers, 0, sizeof XbackwardBuffers);
	
	unsigned char YforwardBuffers[Iz][Ix];
	memset(YforwardBuffers, 0, sizeof YforwardBuffers);
	unsigned char YbackwardBuffers[Iz][Ix];
	memset(YbackwardBuffers, 0, sizeof YbackwardBuffers);

	unsigned char ZforwardBuffers[Ix][Iy];
	memset(ZforwardBuffers, 0, sizeof ZforwardBuffers);
	unsigned char ZbackwardBuffers[Ix][Iy];
	memset(ZbackwardBuffers, 0, sizeof ZbackwardBuffers);


	/* Edges */
	unsigned char XYBuffers[4][Iz];
	memset(XYBuffers, 0, sizeof XYBuffers);
	unsigned char YZBuffers[4][Ix];
	memset(YZBuffers, 0, sizeof YZBuffers);
	unsigned char ZXBuffers[4][Iy];
	memset(ZXBuffers, 0, sizeof ZXBuffers);

	/* corners are enumerated in binary as their coordinates [x, y, y]  */
	unsigned char cornerBuffers[8];
	memset(cornerBuffers, 0, sizeof cornerBuffers);

	//printf("generate data\n");
	// inital image
	srandom(rank+1);
	int x,y,z;
	// for (y=0;y<Iy;y++) {
	// 	for (x=0;x<Ix;x++) {
	// 		for (z=0;z<Iz;z++) {
	// 			data[x][y][z] = random() % 2;
	// 		}
	// 	}
	// }

	if (rank == 0) {

		// data[2][2][2] = 1;
		// data[3][3][2] = 1;
		// data[2][3][3] = 1;
		// data[3][2][3] = 1;
	}
	if (rank == 1) {
		// XforwardBuffers[1][1] = 1;
		// XforwardBuffers[1][2] = 1;
		// XforwardBuffers[1][3] = 1;
		// XforwardBuffers[2][2] = 1;
		// data[4][1][1] = 1;
		// data[4][1][2] = 2;
		// data[4][1][3] = 3;
		// data[4][2][2] = 4;
		 unsigned char uop = 0;
		// for (y=0;y<Iy;y++) {
			
		// 		for (z=0;z<Iz;z++) {
		// 			newData[4][y][z] = ++uop;
		// 			newData[0][y][z] = ++uop;
		// 		}
			
		// }
		// for (x=0;x<Ix;x++) {
			
		// 		for (z=0;z<Iz;z++) {
		// 			newData[x][4][z] = ++uop;
		// 			newData[x][0][z] = ++uop;
		// 		}
			
		// }
		// for (x=0;x<Ix;x++) {
			
		// 		for (y=0;y<Iy;y++) {
		// 			newData[x][y][4] = ++uop;
		// 			newData[x][y][0] = ++uop;
		// 		}
			
		// }
	}


	newData[2][2][0] = 1;
	newData[3][3][0] = 1;
	newData[2][3][4] = 1;
	newData[3][2][4] = 1;

	newData[2][0][2] = 1;
	newData[3][4][2] = 1;
	newData[2][4][3] = 1;
	newData[3][0][3] = 1;

	newData[0][2][2] = 1;
	newData[4][3][2] = 1;
	newData[0][3][3] = 1;
	newData[4][2][3] = 1;

	copyData(Ix,Iy,Iz, data, newData,
					 0, totalProcessorResults );
	// for (x=0;x<Ix;x++) {
	// 	for (y=0;y<Iy;y++) {
	// 		for (z=0;z<Iz;z++) {
	// 			newData[x][y][z] = data[x+1][y+1][z+1]
	// 			totalProcessorResults[x][y][z][0] = newData[x][y][z];
	// 		}
	// 	}
	// }

	//printf("side ranks\n");
	/* Side ranks */
	check ( MPI_Cart_shift(TORUS_COMM, X, FORWARD, &Xforward, &Xbackward) );
	check ( MPI_Cart_shift(TORUS_COMM, Y, FORWARD, &Yforward, &Ybackward) );
	check ( MPI_Cart_shift(TORUS_COMM, Z, FORWARD, &Zforward, &Zbackward) );

	
	//printf("createn data type\n");
	/* Define new datatypes for communications */
	MPI_Datatype  Xside, XsideRecv, Yside, YsideRecv, Zside, ZforwardRecv, ZbackwardRecv; 
	// int MPI_Type_vector(int count, int blocklength, int stride, MPI_Datatype oldtype, MPI_Datatype *newtype)


	/* Custom datatypes for sending and receiving X and Y sides */
	MPI_Type_vector(Iy, Iz, Ix, MPI_UNSIGNED_CHAR, &Xside);
	MPI_Type_commit(&Xside);

	MPI_Type_vector(Iy, Iz, Ix+2, MPI_UNSIGNED_CHAR, &XsideRecv);
	MPI_Type_commit(&XsideRecv);

	MPI_Type_vector(Ix, Iz, Iz*Iy, MPI_UNSIGNED_CHAR, &Yside);
	MPI_Type_commit(&Yside);

	MPI_Type_vector(Ix, Iz, (Iz+2)*(Iy+2), MPI_UNSIGNED_CHAR, &YsideRecv);
	MPI_Type_commit(&YsideRecv);

	/* Custom datatypes for sending and receiving Z sides.  The indexing
	 * of the 3D array maces it impossible to use a vector type for receiving.
	 * We need indexed type for to account for the variable stride.  We even
	 * need different indexed types for the forward and backwards pass.
	 */
	int SIZE = Ix*Iy;
	int blockLengths[Ix*Iy];
	int displacementsBackward[Ix*Iy];
	int displacementsForward[Ix*Iy];
	int disp = (Iy+2)*(Iz+2);
	int index = 0;

	for (x=1;x<Ix+1;x++) {
		for (y=1;y<Iy+1;y++,index++) {
			blockLengths[index] = 1;
			displacementsBackward[index] = x*(Iz+2)*(Iy+2) + y*(Iz+2);
			displacementsForward[index] = x*(Iz+2)*(Iy+2) + y*(Iz+2) + Iz+1;
		}
	}

	MPI_Type_indexed(Ix*Iz, blockLengths, displacementsForward, MPI_UNSIGNED_CHAR, &ZforwardRecv);
	MPI_Type_indexed(Ix*Iz, blockLengths, displacementsBackward, MPI_UNSIGNED_CHAR, &ZbackwardRecv);
	MPI_Type_vector(Ix*Iy, 1, Iy, MPI_UNSIGNED_CHAR, &Zside);

	MPI_Type_commit(&Zside);
	MPI_Type_commit(&ZbackwardRecv);
	MPI_Type_commit(&ZforwardRecv);


	/* Request arrays for the wait command 
	 * [0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25]
	 * [xf xb yf yb zf zb] */
	MPI_Request sendRequests[26];
	MPI_Request recvRequests[26];
	if (rank == 1) {
		// for (x=0;x<Ix;x++) {
		// 	for (y=0;y<Iy;y++) {
		// 		for (z=0;z<Iz;z++) {
		// 			//printf("%d, ",XbackwardGhosts[y][z]);
		// 			printf("%d, ",newData[x][y][z]);
		// 		}
		// 		printf("\n");
		// 	}
		// 	printf("\n");
		// }
		// for (y=0;y<Iy;y++) {
		// 			for (x=0;x<Ix;x++) {
		// 				for (z=0;z<Iz;z++) {
		// 					// printf("%d, ",XbackwardGhosts[y][z]);
		// 					printf("%d, ",newData[x][y][z]);
		// 				}
		// 				printf("\n");
		// 			}
		// 			printf("\n");
		// 		}

		for (z=0;z<Iz;z++) {
			for (y=0;y<Iy;y++) {
				for (x=0;x<Ix;x++) {
				
					// printf("%d, ",XbackwardGhosts[y][z]);
					printf("%d, ",newData[x][y][z]);
				}
				printf("\n");
			}
			printf("\n");
		}
		for (z=0;z<Iz+2;z++) {
			for (y=0;y<Iy+2;y++) {
				for (x=0;x<Ix+2;x++) {
				
					// printf("%d, ",XbackwardGhosts[y][z]);
					printf("%d, ",data[x][y][z]);
				}
				printf("\n");
			}
			printf("\n");
		}
	}
	
	for (i=1;i<=iterations;i++) {
		if (doMPI==1) {
			
			// if (rank==1 || rank==5) {
			// printf("X: Rank %d: (%d,%d,%d) sends to %d, receives from %d\n",rank,coord[0],coord[1],coord[2],Xforward, Xbackward);
			// printf("Y: Rank %d: (%d,%d,%d) sends to %d, receives from %d\n",rank,coord[0],coord[1],coord[2],Yforward, Ybackward);
			// printf("Z: Rank %d: (%d,%d,%d) sends to %d, receives from %d\n",rank,coord[0],coord[1],coord[2],Zforward, Zbackward);

			/* X side communications */
			MPI_Isend(&newData[Ix-1][0][0], 1, Xside, Xforward, 1, TORUS_COMM, &sendRequests[0]);
			MPI_Irecv(&data[0][1][1], 1, XsideRecv, Xbackward, 1, TORUS_COMM, &recvRequests[1]);

			MPI_Isend(&newData[0][0][0], 1, Xside, Xbackward, 1, TORUS_COMM, &sendRequests[1]);
			MPI_Irecv(&data[Ix+1][1][1], 1, XsideRecv, Xforward, 1, TORUS_COMM, &recvRequests[0]);

			/* Y side communications */
			MPI_Isend(&newData[0][Iy-1][0], 1, Yside, Yforward, 1, TORUS_COMM, &sendRequests[2]);
			MPI_Irecv(&data[1][0][1], 1, YsideRecv, Ybackward, 1, TORUS_COMM, &recvRequests[3]);

			MPI_Isend(&newData[0][0][0], 1, Yside, Ybackward, 1, TORUS_COMM, &sendRequests[3]);
			MPI_Irecv(&data[1][Iy+1][1], 1, YsideRecv, Yforward, 1, TORUS_COMM, &recvRequests[2]);

			/* Z side communications */
			MPI_Isend(&newData[0][0][Iz-1], 1, Zside, Zforward, 1, TORUS_COMM, &sendRequests[4]);
			MPI_Irecv(&data[0][0][0], 1, ZbackwardRecv, Zbackward, 1, TORUS_COMM, &recvRequests[5]);

			MPI_Isend(&newData[0][0][0], 1, Zside, Zbackward, 1, TORUS_COMM, &sendRequests[5]);
			MPI_Irecv(&data[0][0][0], 1, ZforwardRecv, Zforward, 1, TORUS_COMM, &recvRequests[4]);
			

			updateCenter(Ix,Iy,Iz, data, newData, rank);

			
			
			

			MPI_Wait(&sendRequests[0], MPI_STATUS_IGNORE);
			MPI_Wait(&sendRequests[1], MPI_STATUS_IGNORE);
			MPI_Wait(&recvRequests[0], MPI_STATUS_IGNORE);
			MPI_Wait(&recvRequests[1], MPI_STATUS_IGNORE);

			updateXside(Ix, Iy, Iz, data, newData, 0, rank);
			updateXside(Ix, Iy, Iz, data, newData, Ix-1, rank);

			MPI_Wait(&sendRequests[2], MPI_STATUS_IGNORE);
			MPI_Wait(&sendRequests[3], MPI_STATUS_IGNORE);
			MPI_Wait(&recvRequests[2], MPI_STATUS_IGNORE);
			MPI_Wait(&recvRequests[3], MPI_STATUS_IGNORE);

			updateYside(Ix, Iy, Iz, data, newData, 0, rank);
			updateYside(Ix, Iy, Iz, data, newData, Iy-1, rank);

			MPI_Wait(&sendRequests[4], MPI_STATUS_IGNORE);
			MPI_Wait(&sendRequests[5], MPI_STATUS_IGNORE);
			MPI_Wait(&recvRequests[4], MPI_STATUS_IGNORE);
			MPI_Wait(&recvRequests[5], MPI_STATUS_IGNORE);

			updateZside(Ix, Iy, Iz, data, newData, 0, rank);
			updateZside(Ix, Iy, Iz, data, newData, Iz-1, rank);

			
			if (rank == 0) {
				printf("receiver\n");
				// for (x=0;x<Ix+2;x++) {
				// 	for (y=0;y<Iy+2;y++) {
				// 		for (z=0;z<Iz+2;z++) {
				// 			// printf("%d, ",XbackwardGhosts[y][z]);
				// 			printf("%d, ",data[x][y][z]);
				// 		}
				// 		printf("\n");
				// 	}
				// 	printf("\n");
				// }
				// for (y=0;y<Iy+2;y++) {
				// 	for (x=0;x<Ix+2;x++) {
				// 		for (z=0;z<Iz+2;z++) {
				// 			// printf("%d, ",XbackwardGhosts[y][z]);
				// 			printf("%d, ",data[x][y][z]);
				// 		}
				// 		printf("\n");
				// 	}
				// 	printf("\n");
				// }
				for (z=0;z<Iz+2;z++) {
					for (y=0;y<Iy+2;y++) {
						for (x=0;x<Ix+2;x++) {
						
							// printf("%d, ",XbackwardGhosts[y][z]);
							printf("%d, ",data[x][y][z]);
						}
						printf("\n");
					}
					printf("\n");
				}
			}

			//printf("Rank %d: (%d,%d,%d)  got through 1X\n",rank,coord[0],coord[1],coord[2]);



			copyData(Ix,Iy,Iz, data, newData,
					 i, totalProcessorResults );
			if (rank == 5) {
				// printf("newData\n");
				// for (y=0;y<Iy;y++) {
				// 	for (z=0;z<Iz;z++) {
				// 		printf("%d, ",newData[0][y][z]);
				// 	}
				// 	printf("\n");
				// }
				// printf("data\n");
				// for (y=0;y<Iy;y++) {
				// 	for (z=0;z<Iz;z++) {
				// 		printf("%d, ",data[0][y][z]);
				// 	}
				// 	printf("\n");
				// }
			}
		}
	}


	FILE *fp;
	if (rank==0) {
		char str[25];
		sprintf(str, "data/metainfodata.txt");
		fp = fopen(str, "w");
		fprintf(fp, "%d %d %d %d %d \n",size,iterations+1, xN,yN,zN);
		fclose(fp);

	}	
	
	printf("\n");
	//	 Save the file. 
	char str[25];
	sprintf(str, "data/%ddata.txt",rank);
	fp = fopen(str, "w");
	int t;
	
	int dx = coord[0]*Lx + MIN(coord[0], Rx);
	int dy = coord[1]*Ly + MIN(coord[1], Ry);
	int dz = coord[2]*Lz + MIN(coord[2], Rz);
	for (t=0;t<=iterations;t++) {
		for (x=0;x<Ix;x++) {
			for (y=0;y<Iy;y++) {
				for (z=0;z<Iz;z++) {
					//fprintf(fp, "%d ", totalProcessorResults[x][y][z][t]);
					if ( totalProcessorResults[x][y][z][t]) {
						fprintf(fp, "%d %d %d ", x+dx, y+dy, z+dz);
					}
				}
			}
		}
		fprintf(fp, "\n");
	}

	fclose(fp);
	
	
	MPI_Finalize();
	exit(0);
}














































