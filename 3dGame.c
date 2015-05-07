/* Use MPI */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "updateFunctions.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))

#define iterations 1 /* number of iterations */

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
				   unsigned char newData[Ix+2][Iy+2][Iz+2], 
				  int t, unsigned char totalProcessorResults[Ix][Iy][Iz][iterations+1]) {
	int x, y, z;
	for (x=0;x<Ix+2;x++) {
		for (y=0;y<Iy+2;y++) {
			for (z=0;z<Iz+2;z++) {
				data[x][y][z] = newData[x][y][z];
			}
		}
	}
}

void saveData(int Ix, int Iy, int Iz, unsigned char data[Ix+2][Iy+2][Iz+2],
				   unsigned char newData[Ix+2][Iy+2][Iz+2], 
				  int t, unsigned char totalProcessorResults[Ix][Iy][Iz][iterations+1]) {
	int x, y, z;
	for (x=0;x<Ix;x++) {
		for (y=0;y<Iy;y++) {
			for (z=0;z<Iz;z++) {
				totalProcessorResults[x][y][z][t] = newData[x][y][z];
			}
		}
	}
}





int main(int argc, char *argv[]){

	/* Constants we need for data distribution */
	int Nx, Lx, Rx, Ix;
	int Ny, Ly, Ry, Iy;
	int Nz, Lz, Rz, Iz;

	/* Ranks of the processrs to communicate with */
	int Xforward, Xbackward, Yforward, Ybackward, Zforward, Zbackward; 
	int XYprocs[4];
	int YZprocs[4];
	int ZXprocs[4];
	int corners[8];

    /* helper variables */
	int dx, dy, dz;
	int x, y, z, i;
    
    /* Define new datatypes for communications and virtual topolgy */
    int rank, size;
    int FORWARD = 1, BACKWARDS = -1;
    int dim[3], period[3], reorder;
    int coord[3];
    int doMPI = 1;
	MPI_Comm TORUS_COMM;
	MPI_Datatype Xside, Yside, Zside;
	MPI_Datatype XYedge, YZedge, ZXedge;

	/* Request arrays for the wait command 
	 * [0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25]
	 * [xf xb yf yb zf zb and so on ... ] 
	 */
	MPI_Request sendRequests[26];
	MPI_Request recvRequests[26];

    /* initialize MPI and create virtual 2d torus topology*/

    check ( MPI_Init(&argc, &argv) );
    check ( MPI_Comm_rank(MPI_COMM_WORLD, &rank) );
    check ( MPI_Comm_size(MPI_COMM_WORLD, &size) );

    dim[0] = Px; dim[1] = Px; dim[2] = Px;
	period[0] = xBC; period[1] = yBC; period[2] = zBC;
	reorder = 0;

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

	/* Data arrays */
	unsigned char totalProcessorResults[Ix][Iy][Iz][iterations+1];
	unsigned char data[Ix+2][Iy+2][Iz+2];
	memset(data, 0, sizeof data);
	unsigned char newData[Ix+2][Iy+2][Iz+2];
	memset(newData, 0, sizeof newData);

	/* inital data */
	srandom(rank+1);
	
	// for (y=0;y<Iy;y++) {
	// 	for (x=0;x<Ix;x++) {
	// 		for (z=0;z<Iz;z++) {
	// 			data[x][y][z] = random() % 2;
	// 		}
	// 	}
	// }


	// if (rank == 0) {
	// 	unsigned char uop = 0;
	// 	// for (y=1;y<Iy+1;y++) {
			
	// 	// 		for (z=1;z<Iz+1;z++) {
	// 	// 			newData[5][y][z] = ++uop;
	// 	// 			newData[1][y][z] = ++uop;
	// 	// 		}
			
	// 	// }
	// 	// for (x=1;x<Ix+1;x++) {
			
	// 	// 		for (z=1;z<Iz+1;z++) {
	// 	// 			newData[x][5][z] = ++uop;
	// 	// 			newData[x][1][z] = ++uop;
	// 	// 		}
			
	// 	// }
	// 	for (x=1;x<Ix+1;x++) {
			
	// 			for (y=1;y<Iy+1;y++) {
	// 				newData[x][y][5] = ++uop;
	// 				newData[x][y][1] = ++uop;
	// 			}
			
	// 	}
	// }

		//side oscillating
	// newData[2][2][0] = 1;
	// newData[3][3][0] = 1;
	// newData[2][3][4] = 1;
	// newData[3][2][4] = 1;

	// newData[2][0][2] = 1;
	// newData[3][4][2] = 1;
	// newData[2][4][3] = 1;
	// newData[3][0][3] = 1;

	// newData[0][2][2] = 1;
	// newData[4][3][2] = 1;
	// newData[0][3][3] = 1;
	// newData[4][2][3] = 1;

	// newData[0][0][0] = 1;
	// newData[0][1][0] = 1;
	// newData[0][2][0] = 1;
	// newData[0][3][0] = 1;
	// newData[0][4][0] = 1;

	// newData[0][0][4] = 2;
	// newData[0][1][4] = 2;
	// newData[0][2][4] = 2;
	// newData[0][3][4] = 2;
	// newData[0][4][4] = 2;

	// newData[4][0][0] = 3;
	// newData[4][1][0] = 3;
	// newData[4][2][0] = 3;
	// newData[4][3][0] = 3;
	// newData[4][4][0] = 3;

	// newData[4][0][4] = 4;
	// newData[4][1][4] = 4;
	// newData[4][2][4] = 4;
	// newData[4][3][4] = 4;
	// newData[4][4][4] = 4;

	// newData[1][1][1] = 1;
	// newData[1][1][5] = 2;
	// newData[1][5][1] = 3;
	// newData[1][5][5] = 6;
	// newData[5][1][1] = 5;
	// newData[5][1][5] = 6;
	// newData[5][5][1] = 7;
	// newData[5][5][5] = 8;
	
	


	copyData(Ix,Iy,Iz, data, newData,
					 0, totalProcessorResults );
	saveData(Ix,Iy,Iz, data, newData,
					 0, totalProcessorResults );


	/* Side ranks */
	check ( MPI_Cart_shift(TORUS_COMM, X, FORWARD, &Xforward, &Xbackward) );
	check ( MPI_Cart_shift(TORUS_COMM, Y, FORWARD, &Yforward, &Ybackward) );
	check ( MPI_Cart_shift(TORUS_COMM, Z, FORWARD, &Zforward, &Zbackward) );

	/* get ranks of diagonal processors in the XY plane */
	int procCoords[3] = {coord[0]-1, coord[1]-1, coord[2]};
	check ( MPI_Cart_rank(TORUS_COMM, procCoords, &XYprocs[0]) );
	procCoords[0] = coord[0]-1; procCoords[1] = coord[1]+1;
	check ( MPI_Cart_rank(TORUS_COMM, procCoords, &XYprocs[1]) );
	procCoords[0] = coord[0]+1; procCoords[1] = coord[1]-1;
	check ( MPI_Cart_rank(TORUS_COMM, procCoords, &XYprocs[2]) );
	procCoords[0] = coord[0]+1; procCoords[1] = coord[1]+1;
	check ( MPI_Cart_rank(TORUS_COMM, procCoords, &XYprocs[3]) );

	/* get ranks of diagonal processors in the YZ plane */
	procCoords[1] = coord[1]-1; procCoords[2] = coord[2]-1;
	check ( MPI_Cart_rank(TORUS_COMM, procCoords, &YZprocs[0]) );
	procCoords[1] = coord[1]-1; procCoords[2] = coord[2]+1;
	check ( MPI_Cart_rank(TORUS_COMM, procCoords, &YZprocs[1]) );
	procCoords[1] = coord[1]+1; procCoords[2] = coord[2]-1;
	check ( MPI_Cart_rank(TORUS_COMM, procCoords, &YZprocs[2]) );
	procCoords[1] = coord[1]+1; procCoords[2] = coord[2]+1;
	check ( MPI_Cart_rank(TORUS_COMM, procCoords, &YZprocs[3]) );

	/* get ranks of diagonal processors in the ZX plane */
	procCoords[0] = coord[0]-1; procCoords[2] = coord[2]-1;
	check ( MPI_Cart_rank(TORUS_COMM, procCoords, &ZXprocs[0]) );
	procCoords[0] = coord[0]-1; procCoords[2] = coord[2]+1;
	check ( MPI_Cart_rank(TORUS_COMM, procCoords, &ZXprocs[1]) );
	procCoords[0] = coord[0]+1; procCoords[2] = coord[2]-1;
	check ( MPI_Cart_rank(TORUS_COMM, procCoords, &ZXprocs[2]) );
	procCoords[0] = coord[0]+1; procCoords[2] = coord[2]+1;
	check ( MPI_Cart_rank(TORUS_COMM, procCoords, &ZXprocs[3]) );

	/* Get the ranks of 8 corner processors */
	i=0;
	for (dx=-1;dx<=1;dx=dx+2) {
		for (dy=-1;dy<=1;dy=dy+2) {
			for (dz=-1;dz<=1;dz=dz+2) {
				procCoords[0] = coord[0]+dx; procCoords[1] = coord[1]+dy; procCoords[2] = coord[2]+dz; 
				check( MPI_Cart_rank(TORUS_COMM, procCoords, &corners[i]) );
				i++;
			}
		}	
	}

	/* Custom datatypes for sending and receiving X and Y sides */
	// MPI_Type_vector(int count, int blocklength, int stride, MPI_Datatype oldtype, MPI_Datatype *newtype)
	check ( MPI_Type_vector(Iy, Iz, Ix+2, MPI_UNSIGNED_CHAR, &Xside) );
	check ( MPI_Type_commit(&Xside) );

	check ( MPI_Type_vector(Ix, Iz, (Iz+2)*(Iy+2), MPI_UNSIGNED_CHAR, &Yside) );
	check ( MPI_Type_commit(&Yside) );

	/* Custom datatypes for sending and receiving Z sides.  The indexing
	 * of the 3D array maces it impossible to use a vector type for receiving.
	 * We need indexed type for to account for the variable stride.
	 */
	
	int blockLengths[Ix*Iy];
	int displacements[Ix*Iy];
	int index = 0;

	for (x=0;x<Ix;x++) {
		for (y=0;y<Iy;y++,index++) {
			blockLengths[index] = 1;
			displacements[index] = x*(Iz+2)*(Iy+2) + y*(Iz+2);
		}
	}

	check ( MPI_Type_indexed(Ix*Iz, blockLengths, displacements, MPI_UNSIGNED_CHAR, &Zside) );
	check ( MPI_Type_commit(&Zside) );

	/* Custom data types for edges */
	/* XY edge */
	check ( MPI_Type_vector(1, Iz, Iz+2, MPI_UNSIGNED_CHAR, &XYedge) );
	check ( MPI_Type_commit(&XYedge) );

	/* YZ edge */
	check ( MPI_Type_vector(Ix, 1, (Iz+2)*(Iy+2), MPI_UNSIGNED_CHAR, &YZedge) );
	check ( MPI_Type_commit(&YZedge) );

	/* ZX edge */
	check ( MPI_Type_vector(Iy, 1, (Iz+2), MPI_UNSIGNED_CHAR, &ZXedge) );
	check ( MPI_Type_commit(&ZXedge) );

	if (rank == 0) {
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

		// 		for (x=0;x<Ix+2;x++) {
		// for (z=0;z<Iz+2;z++) {
		// 	for (y=0;y<Iy+2;y++) {
				
		// 			// printf("%d, ",XbackwardGhosts[y][z]);
		// 			printf("%d, ",newData[x][y][z]);
		// 		}
		// 		printf("\n");
		// 	}
		// 	printf("\n");
		// }

		// for (z=0;z<Iz+2;z++) {
		// 	for (y=0;y<Iy+2;y++) {
		// 		for (x=0;x<Ix+2;x++) {
				
		// 			// printf("%d, ",XbackwardGhosts[y][z]);
		// 			printf("%d, ",data[x][y][z]);
		// 		}
		// 		printf("\n");
		// 	}
		// 	printf("\n");
		// }
	}
	
	for (i=1;i<=iterations;i++) {
		if (doMPI==1) {
			
			
			/* X side communications */
			check ( MPI_Isend(&data[Ix][1][1], 1, Xside, Xforward, 0, TORUS_COMM, &sendRequests[0]) );
			check ( MPI_Irecv(&data[0][1][1], 1, Xside, Xbackward, 0, TORUS_COMM, &recvRequests[1]) );

			check ( MPI_Isend(&data[1][1][1], 1, Xside, Xbackward, 1, TORUS_COMM, &sendRequests[1]) );
			check ( MPI_Irecv(&data[Ix+1][1][1], 1, Xside, Xforward, 1, TORUS_COMM, &recvRequests[0]) );

			/* Y side communications */
			check ( MPI_Isend(&data[1][Iy][1], 1, Yside, Yforward, 2, TORUS_COMM, &sendRequests[2]) );
			check ( MPI_Irecv(&data[1][0][1], 1, Yside, Ybackward, 2, TORUS_COMM, &recvRequests[3]) );

			check ( MPI_Isend(&data[1][1][1], 1, Yside, Ybackward, 3, TORUS_COMM, &sendRequests[3]) );
			check ( MPI_Irecv(&data[1][Iy+1][1], 1, Yside, Yforward, 3, TORUS_COMM, &recvRequests[2]) );


			/* Z side communications */
			check ( MPI_Isend(&data[1][1][Iz], 1, Zside, Zforward, 4, TORUS_COMM, &sendRequests[4]) );
			check ( MPI_Irecv(&data[1][1][0], 1, Zside, Zbackward, 4, TORUS_COMM, &recvRequests[5]) );

			check ( MPI_Isend(&data[1][1][1], 1, Zside, Zbackward, 5, TORUS_COMM, &sendRequests[5]) );
			check ( MPI_Irecv(&data[1][1][Iz+1], 1, Zside, Zforward, 5, TORUS_COMM, &recvRequests[4]) );

			/* XY edges communication */
			check ( MPI_Isend(&data[1][1][1], 1, XYedge, XYprocs[0], 1, TORUS_COMM, &sendRequests[6]) );
			check ( MPI_Isend(&data[1][Iy][1], 1, XYedge, XYprocs[1], 1, TORUS_COMM, &sendRequests[7]) );
			check ( MPI_Isend(&data[Ix][1][1], 1, XYedge, XYprocs[2], 1, TORUS_COMM, &sendRequests[8]) );
			check ( MPI_Isend(&data[Ix][Iy][1], 1, XYedge, XYprocs[3], 1, TORUS_COMM, &sendRequests[9]) );

			check ( MPI_Irecv(&data[Ix+1][Iy+1][1], 1, XYedge, XYprocs[3], 1, TORUS_COMM, &recvRequests[6]) );
			check ( MPI_Irecv(&data[Ix+1][0][1], 1, XYedge, XYprocs[2], 1, TORUS_COMM, &recvRequests[7]) );
			check ( MPI_Irecv(&data[0][Iy+1][1], 1, XYedge, XYprocs[1], 1, TORUS_COMM, &recvRequests[8]) );
			check ( MPI_Irecv(&data[0][0][1], 1, XYedge, XYprocs[0], 1, TORUS_COMM, &recvRequests[9]) );

			/* YZ edge commuinication */
			check ( MPI_Isend(&data[1][1][1], 1, YZedge, YZprocs[0], 1, TORUS_COMM, &sendRequests[10]) );
			check ( MPI_Isend(&data[1][1][Iz], 1, YZedge, YZprocs[1], 1, TORUS_COMM, &sendRequests[11]) );
			check ( MPI_Isend(&data[1][Iy][1], 1, YZedge, YZprocs[2], 1, TORUS_COMM, &sendRequests[12]) );
			check ( MPI_Isend(&data[1][Iy][Iz], 1, YZedge, YZprocs[3], 1, TORUS_COMM, &sendRequests[13]) );

			check ( MPI_Irecv(&data[1][Iy+1][Iz+1], 1, YZedge, YZprocs[3], 1, TORUS_COMM, &recvRequests[10]) );
			check ( MPI_Irecv(&data[1][Iy+1][0], 1, YZedge, YZprocs[2], 1, TORUS_COMM, &recvRequests[11]) );
			check ( MPI_Irecv(&data[1][0][Iz+1], 1, YZedge, YZprocs[1], 1, TORUS_COMM, &recvRequests[12]) );
			check ( MPI_Irecv(&data[1][0][0], 1, YZedge, YZprocs[0], 1, TORUS_COMM, &recvRequests[13]) );

			/* ZX edge communication */
			check ( MPI_Isend(&data[1][1][1], 1, ZXedge, ZXprocs[0], 1, TORUS_COMM, &sendRequests[14]) );
			check ( MPI_Isend(&data[1][1][Iz], 1, ZXedge, ZXprocs[1], 1, TORUS_COMM, &sendRequests[15]) );
			check ( MPI_Isend(&data[Ix][1][1], 1, ZXedge, ZXprocs[2], 1, TORUS_COMM, &sendRequests[16]) );
			check ( MPI_Isend(&data[Ix][1][Iz], 1, ZXedge, ZXprocs[3], 1, TORUS_COMM, &sendRequests[17]) );

			chech ( MPI_Irecv(&data[Ix+1][1][Iz+1], 1, ZXedge, ZXprocs[3], 1, TORUS_COMM, &recvRequests[14]) );
			check ( MPI_Irecv(&data[Ix+1][1][0], 1, ZXedge, ZXprocs[2], 1, TORUS_COMM, &recvRequests[15]) );
			check ( MPI_Irecv(&data[0][1][Iz+1], 1, ZXedge, ZXprocs[1], 1, TORUS_COMM, &recvRequests[16]) );
			check ( MPI_Irecv(&data[0][1][0], 1, ZXedge, ZXprocs[0], 1, TORUS_COMM, &recvRequests[17]) );

			/* Corner communications */
			check ( MPI_Isend(&data[1][1][1], 1, MPI_UNSIGNED_CHAR, corners[0], 18, TORUS_COMM, &sendRequests[18]) );
			check ( MPI_Isend(&data[1][1][Iz], 1, MPI_UNSIGNED_CHAR, corners[1], 19, TORUS_COMM, &sendRequests[19]) );
			check ( MPI_Isend(&data[1][Iy][1], 1, MPI_UNSIGNED_CHAR, corners[2], 20, TORUS_COMM, &sendRequests[20]) );
			check ( MPI_Isend(&data[1][Iy][Iz], 1, MPI_UNSIGNED_CHAR, corners[3], 21, TORUS_COMM, &sendRequests[21]) );
			check ( MPI_Isend(&data[Ix][1][1], 1, MPI_UNSIGNED_CHAR, corners[4], 22, TORUS_COMM, &sendRequests[22]) );
			check ( MPI_Isend(&data[Ix][1][Iz], 1, MPI_UNSIGNED_CHAR, corners[5], 23, TORUS_COMM, &sendRequests[23]) );
			check ( MPI_Isend(&data[Ix][Iy][1], 1, MPI_UNSIGNED_CHAR, corners[6], 24, TORUS_COMM, &sendRequests[24]) );
			check ( MPI_Isend(&data[Ix][Iy][Iz], 1, MPI_UNSIGNED_CHAR, corners[7], 25, TORUS_COMM, &sendRequests[25]) );

			check ( MPI_Irecv(&data[Ix+1][Iy+1][Iz+1], 1, MPI_UNSIGNED_CHAR, corners[7], 18, TORUS_COMM, &recvRequests[18]) );
			check ( MPI_Irecv(&data[Ix+1][Iy+1][0], 1, MPI_UNSIGNED_CHAR, corners[6], 19, TORUS_COMM, &recvRequests[19]) );
			check ( MPI_Irecv(&data[Ix+1][0][Iz+1], 1, MPI_UNSIGNED_CHAR, corners[5], 20, TORUS_COMM, &recvRequests[20]) );
			check ( MPI_Irecv(&data[Ix+1][0][0], 1, MPI_UNSIGNED_CHAR, corners[4], 21, TORUS_COMM, &recvRequests[21]) );
			check ( MPI_Irecv(&data[0][Iy+1][Iz+1], 1, MPI_UNSIGNED_CHAR, corners[3], 22, TORUS_COMM, &recvRequests[22]) );
			check ( MPI_Irecv(&data[0][Iy+1][0], 1, MPI_UNSIGNED_CHAR, corners[2], 23, TORUS_COMM, &recvRequests[23]) );
			check ( MPI_Irecv(&data[0][0][Iz+1], 1, MPI_UNSIGNED_CHAR, corners[1], 24, TORUS_COMM, &recvRequests[24]) );
			check ( MPI_Irecv(&data[0][0][0], 1, MPI_UNSIGNED_CHAR, corners[0], 25, TORUS_COMM, &recvRequests[25]) );

			/* Start updating the center */
			updateCenter(Ix,Iy,Iz, data, newData, rank);

			/* Wait until a forward Xside has been received and update it */
			check ( MPI_Wait(&recvRequests[0], MPI_STATUS_IGNORE);
			updateXside(Ix, Iy, Iz, data, newData, Ix, rank);

			/* Wait until a backward Xside has been received and update it */
			check ( MPI_Wait(&recvRequests[1], MPI_STATUS_IGNORE);
			updateXside(Ix, Iy, Iz, data, newData, 1, rank);

			/* Wait until a forward Yside has been received and update it */
			check ( MPI_Wait(&recvRequests[2], MPI_STATUS_IGNORE);
			updateYside(Ix, Iy, Iz, data, newData, Iy, rank);

			/* Wait until a backward Yside has been received and update it */
			check ( MPI_Wait(&recvRequests[3], MPI_STATUS_IGNORE);
			updateYside(Ix, Iy, Iz, data, newData, 1, rank);

			/* Wait until a forward Zside has been received and update it */
			check ( MPI_Wait(&recvRequests[4], MPI_STATUS_IGNORE);
			updateZside(Ix, Iy, Iz, data, newData, Iy, rank);

			/* Wait until a backward Zside has been received and update it */
			check ( MPI_Wait(&recvRequests[5], MPI_STATUS_IGNORE);
			updateZside(Ix, Iy, Iz, data, newData, 1, rank);

			/* Wait and Update the XY edges */

			check ( MPI_Wait(&recvRequests[6], MPI_STATUS_IGNORE) );
			updateXYedge(Ix, Iy, Iz, data, newData, Ix, Iy, rank);

			check ( MPI_Wait(&recvRequests[7], MPI_STATUS_IGNORE) );
			updateXYedge(Ix, Iy, Iz, data, newData, Ix, 1, rank);

			check ( MPI_Wait(&recvRequests[8], MPI_STATUS_IGNORE) );
			updateXYedge(Ix, Iy, Iz, data, newData, 1, Iy, rank);

			check ( MPI_Wait(&recvRequests[9], MPI_STATUS_IGNORE) );
			updateXYedge(Ix, Iy, Iz, data, newData, 1, 1, rank);

			/* Wait and Update the YZ edges */

			check ( MPI_Wait(&recvRequests[10], MPI_STATUS_IGNORE) );
			updateYZedge(Ix, Iy, Iz, data, newData, Iy, Iz, rank);

			check ( MPI_Wait(&recvRequests[11], MPI_STATUS_IGNORE) );
			updateYZedge(Ix, Iy, Iz, data, newData, Iy, 1, rank);

			check ( MPI_Wait(&recvRequests[12], MPI_STATUS_IGNORE) );
			updateYZedge(Ix, Iy, Iz, data, newData, 1, Iz, rank);

			check ( MPI_Wait(&recvRequests[13], MPI_STATUS_IGNORE) );
			updateYZedge(Ix, Iy, Iz, data, newData, 1, 1, rank);

			/* Wait and Update the ZX edges */

			check ( MPI_Wait(&recvRequests[14], MPI_STATUS_IGNORE) );
			updateZXedge(Ix, Iy, Iz, data, newData, Iz, Ix, rank);

			check ( MPI_Wait(&recvRequests[15], MPI_STATUS_IGNORE) );
			updateZXedge(Ix, Iy, Iz, data, newData, 1, Ix, rank);

			check ( MPI_Wait(&recvRequests[16], MPI_STATUS_IGNORE) );
			updateZXedge(Ix, Iy, Iz, data, newData, Iz, 1, rank);

			check ( MPI_Wait(&recvRequests[17], MPI_STATUS_IGNORE) );
			updateZXedge(Ix, Iy, Iz, data, newData, 1, 1, rank);
			
			/* Wait for the corners to come in before updating them */
			
			check ( MPI_Wait(&recvRequests[18], MPI_STATUS_IGNORE) );
			check ( MPI_Wait(&recvRequests[19], MPI_STATUS_IGNORE) );
			check ( MPI_Wait(&recvRequests[20], MPI_STATUS_IGNORE) );
			check ( MPI_Wait(&recvRequests[21], MPI_STATUS_IGNORE) );
			check ( MPI_Wait(&recvRequests[22], MPI_STATUS_IGNORE) );
			check ( MPI_Wait(&recvRequests[23], MPI_STATUS_IGNORE) );
			check ( MPI_Wait(&recvRequests[24], MPI_STATUS_IGNORE) );
			check ( MPI_Wait(&recvRequests[25], MPI_STATUS_IGNORE) );

			updateCorners(Ix, Iy, Iz, data, newData, rank)

			check ( MPI_Waitall(26, sendRequests,  MPI_STATUSES_IGNORE) );

			
			// if (rank == 1) {
			// 	printf("receiver\n");
			// 	// for (x=0;x<Ix+2;x++) {
			// 	// 	for (y=0;y<Iy+2;y++) {
			// 	// 		for (z=0;z<Iz+2;z++) {
			// 	// 			// printf("%d, ",XbackwardGhosts[y][z]);
			// 	// 			printf("%d, ",data[x][y][z]);
			// 	// 		}
			// 	// 		printf("\n");
			// 	// 	}
			// 	// 	printf("\n");
			// 	// }
			// 	// for (y=0;y<Iy+2;y++) {
			// 	// 	for (x=0;x<Ix+2;x++) {
			// 	// 		for (z=0;z<Iz+2;z++) {
			// 	// 			// printf("%d, ",XbackwardGhosts[y][z]);
			// 	// 			printf("%d, ",data[x][y][z]);
			// 	// 		}
			// 	// 		printf("\n");
			// 	// 	}
			// 	// 	printf("\n");
			// 	// }
			// 	for (z=0;z<Iz+2;z++) {
			// 		for (y=0;y<Iy+2;y++) {
			// 			for (x=0;x<Ix+2;x++) {
						
			// 				// printf("%d, ",XbackwardGhosts[y][z]);
			// 				printf("%d, ",data[x][y][z]);
			// 			}
			// 			printf("\n");
			// 		}
			// 		printf("\n");
			// 	}
			// }



			// copyData(Ix,Iy,Iz, data, newData,
			// 		 i, totalProcessorResults );
			unsigned char * tmp = data;
			data = newData;
			newData = tmp;
			saveData(Ix,Iy,Iz, data, newData,
					 i, totalProcessorResults );

			// if (rank == 0) {
			// 	printf("receiver\n");
			// 	// for (x=0;x<Ix+2;x++) {
			// 	// 	for (y=0;y<Iy+2;y++) {
			// 	// 		for (z=0;z<Iz+2;z++) {
			// 	// 			// printf("%d, ",XbackwardGhosts[y][z]);
			// 	// 			printf("%d, ",data[x][y][z]);
			// 	// 		}
			// 	// 		printf("\n");
			// 	// 	}
			// 	// 	printf("\n");
			// 	// }
			// 	// for (y=0;y<Iy+2;y++) {
			// 	// 	for (x=0;x<Ix+2;x++) {
			// 	// 		for (z=0;z<Iz+2;z++) {
			// 	// 			// printf("%d, ",XbackwardGhosts[y][z]);
			// 	// 			printf("%d, ",data[x][y][z]);
			// 	// 		}
			// 	// 		printf("\n");
			// 	// 	}
			// 	// 	printf("\n");
			// 	// }
			// 	for (z=0;z<Iz+2;z++) {
			// 		for (y=0;y<Iy+2;y++) {
			// 			for (x=0;x<Ix+2;x++) {
						
			// 				// printf("%d, ",XbackwardGhosts[y][z]);
			// 				printf("%d, ",data[x][y][z]);
			// 			}
			// 			printf("\n");
			// 		}
			// 		printf("\n");
			// 	}
			// }

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
	
	dx = coord[0]*Lx + MIN(coord[0], Rx);
	dy = coord[1]*Ly + MIN(coord[1], Ry);
	dz = coord[2]*Lz + MIN(coord[2], Rz);
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














































