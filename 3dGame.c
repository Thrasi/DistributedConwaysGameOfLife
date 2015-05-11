/* Use MPI */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "updateFunctions.h"
#include "optimization.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))

#define iterations 60 /* number of iterations */

/* Number of processors in each dimension */
// #define Px 2
// #define Py 2
// #define Pz 2

/* Periodicity of boundaries, should always be 1 */
#define xBC 1
#define yBC 1
#define zBC 1

/* Size of each dimension */
#define xN 100
#define yN 100
#define zN 100

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
				   unsigned char newData[Ix+2][Iy+2][Iz+2]) {
	int x, y, z;
	for (x=0;x<Ix+2;x++) {
		for (y=0;y<Iy+2;y++) {
			for (z=0;z<Iz+2;z++) {
				data[x][y][z] = newData[x][y][z];
			}
		}
	}
}

void saveData(int Ix, int Iy, int Iz, unsigned char newData[Ix+2][Iy+2][Iz+2], 
				  int t, unsigned char totalProcessorResults[Ix][Iy][Iz][iterations+1]) {
	int x, y, z;
	for (x=0;x<Ix;x++) {
		for (y=0;y<Iy;y++) {
			for (z=0;z<Iz;z++) {
				// printf("aaa %d %d %d %d %d\n",x,y,z,t,newData[x+1][y+1][z+1]);
				totalProcessorResults[x][y][z][t] = newData[x+1][y+1][z+1];
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

    /* initialize MPI and create virtual 3d torus topology*/
    double start = MPI_Wtime();
    check ( MPI_Init(&argc, &argv) );
    check ( MPI_Comm_rank(MPI_COMM_WORLD, &rank) );
    check ( MPI_Comm_size(MPI_COMM_WORLD, &size) );
    struct triple opt = optimize(xN, yN, zN, size);
   
   	int Px = opt.p, Py = opt.q, Pz = opt.r;
   	printf("%d %d %d\n",Px, Py, Pz);
    // if (rank==0)
    // {
    // printf("rank: %d, size: %d, == %d\n",rank, size, Px*Py*Pz);
    // printf("rank: %d, size: %d, == %d\n",rank, size, Px*Py*Pz);
    // }

    dim[0] = Px; dim[1] = Py; dim[2] = Pz;
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
	// printf("rank: %d, %d,%d,%d",rank, Ix,Iy,Iz);
	// MPI_Finalize();
	// exit(0);
	/* Data arrays */
	unsigned char totalProcessorResults[Ix][Iy][Iz][iterations+1];
	unsigned char data[Ix+2][Iy+2][Iz+2];
	memset(data, 0, sizeof data);
	unsigned char newData[Ix+2][Iy+2][Iz+2];
	memset(newData, 0, sizeof newData);

	/* inital data */
	srandom(rank+1);
	if (rank==0) {
		for (y=1;y<Iy+1;y++) {
			for (x=1;x<Ix+1;x++) {
				for (z=1;z<Iz+1;z++) {
					newData[x][y][z] = random() % 2;
					//newData[x][y][1] = 1;

				}
			}
		}

		
		

		
		newData[1][5][2] = 1;
		newData[2][5][2] = 1;

		newData[1][4][1] = 1;
		newData[2][4][1] = 1;

		newData[3][5][1] = 1;	

		// newData[2][3][1] = 1;
		// newData[2][3][2] = 1;

		// newData[6][6][6] = 1;
		// newData[6][6][7] = 1;

		// newData[6][9][6] = 1;
		// newData[6][9][7] = 1;

		// newData[6][7][8] = 1;
		// newData[6][8][8] = 1;

		// newData[7][7][6] = 1;
		// newData[7][7][7] = 1;
		// newData[7][8][6] = 1;
		// newData[7][8][7] = 1;
	}

	if (rank==1) {
		newData[1][1][Iz] = 1;

	}
	if (rank==2) {
		newData[1][Iy][1] = 1;

		newData[1][1][2] = 1;
		newData[2][1][2] = 1;

		newData[1][2][1] = 1;
		newData[2][2][1] = 1;

		newData[3][1][1] = 1;
	}
	if (rank==4) {
		newData[Ix][1][1] = 1;
	}
	if (rank==7) {
		newData[Ix][Iy][Iz] = 1;
	}
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
	//if (rank == 1)
		//side oscillating
	// newData[1][1][1] = 1;
	// newData[5][5][1] = 1;
	// newData[1][5][5] = 1;
	// newData[5][1][5] = 1;

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
	
	


	copyData(Ix,Iy,Iz, data, newData );
	// printf("rank: %d, [1][1][1]: %d\n",rank, newData[1][1][1]);
	saveData(Ix,Iy,Iz, newData, 0, totalProcessorResults );
	// printf("rank: %d, [1][1][1]: %d\n",rank, newData[1][1][1]);
	// printf("rank: %d, THIS: %d\n",rank, totalProcessorResults[0][0][0][0]);


	/* Side ranks */
	check ( MPI_Cart_shift(TORUS_COMM, X, BACKWARDS, &Xforward, &Xbackward) );
	check ( MPI_Cart_shift(TORUS_COMM, Y, BACKWARDS, &Yforward, &Ybackward) );
	check ( MPI_Cart_shift(TORUS_COMM, Z, BACKWARDS, &Zforward, &Zbackward) );

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
	procCoords[0] = coord[0]; procCoords[1] = coord[1]; procCoords[2] = coord[2];  // re initialize the offset coordinates
	procCoords[1] = coord[1]-1; procCoords[2] = coord[2]-1;
	check ( MPI_Cart_rank(TORUS_COMM, procCoords, &YZprocs[0]) );
	procCoords[1] = coord[1]-1; procCoords[2] = coord[2]+1;
	check ( MPI_Cart_rank(TORUS_COMM, procCoords, &YZprocs[1]) );
	procCoords[1] = coord[1]+1; procCoords[2] = coord[2]-1;
	check ( MPI_Cart_rank(TORUS_COMM, procCoords, &YZprocs[2]) );
	procCoords[1] = coord[1]+1; procCoords[2] = coord[2]+1;
	check ( MPI_Cart_rank(TORUS_COMM, procCoords, &YZprocs[3]) );

	/* get ranks of diagonal processors in the ZX plane */
	procCoords[0] = coord[0]; procCoords[1] = coord[1]; procCoords[2] = coord[2];
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
	check ( MPI_Type_vector(Iy, Iz, Iz+2, MPI_UNSIGNED_CHAR, &Xside) );
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

	check ( MPI_Type_indexed(Ix*Iy, blockLengths, displacements, MPI_UNSIGNED_CHAR, &Zside) );
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

		// for (z=0;z<Iz+2;z++) {
		// 	for (y=0;y<Iy+2;y++) {
		// 		for (x=0;x<Ix+2;x++) {
				
		// 			// printf("%d, ",XbackwardGhosts[y][z]);
		// 			printf("%d, ",newData[x][y][z]);
		// 		}
		// 		printf("\n");
		// 	}
		// 	printf("\n");
		// }

		// 		for (x=0;x<Ix+2;x++) {
		// for (z=0;z<Iz+2;z++) {
		// 	for (y=0;y<Iy+2;y++) {
				
		// 			// printf("%d, ",XbackwardGhosts[y][z]);
		// 			printf("%d, ",data[x][y][z]);
		// 		}
		// 		printf("\n");
		// 	}
		// 	printf("\n");
		// }
	}
	// if (rank == 0) {
	// 	printf("rank: %d, Before iter: %d\n",rank, totalProcessorResults[0][0][0][0]);	
	// }
	
	for (i=1;i<=iterations;i++) {
		if (doMPI==1) {
			
			
			/* X side communications */
			check ( MPI_Isend(&data[Ix][1][1], 1, Xside, Xforward, 0, TORUS_COMM, &sendRequests[0]) );
			check ( MPI_Irecv(&data[0][1][1], 1, Xside, Xbackward, 0, TORUS_COMM, &recvRequests[1]) );
			// if (i==1) {
			// 	printf("rank %d, Xforward: %d, Xbackward: %d\n",rank,Xforward, Xbackward);
			// }
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

			check ( MPI_Irecv(&data[Ix+1][1][Iz+1], 1, ZXedge, ZXprocs[3], 1, TORUS_COMM, &recvRequests[14]) );
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
			
			

			/* Wait until a forward Xside has been received and update it */
			check ( MPI_Wait(&recvRequests[0], MPI_STATUS_IGNORE) );
			

			/* Wait until a backward Xside has been received and update it */
			check ( MPI_Wait(&recvRequests[1], MPI_STATUS_IGNORE) );
			

			/* Wait until a forward Yside has been received and update it */
			check ( MPI_Wait(&recvRequests[2], MPI_STATUS_IGNORE) );
			

			/* Wait until a backward Yside has been received and update it */
			check ( MPI_Wait(&recvRequests[3], MPI_STATUS_IGNORE) );
			

			/* Wait until a forward Zside has been received and update it */
			check ( MPI_Wait(&recvRequests[4], MPI_STATUS_IGNORE) );
			

			/* Wait until a backward Zside has been received and update it */
			check ( MPI_Wait(&recvRequests[5], MPI_STATUS_IGNORE) );
			



			/* Wait and Update the XY edges */

			check ( MPI_Wait(&recvRequests[6], MPI_STATUS_IGNORE) );
			

			check ( MPI_Wait(&recvRequests[7], MPI_STATUS_IGNORE) );
			

			check ( MPI_Wait(&recvRequests[8], MPI_STATUS_IGNORE) );
			

			check ( MPI_Wait(&recvRequests[9], MPI_STATUS_IGNORE) );
			

			/* Wait and Update the YZ edges */

			check ( MPI_Wait(&recvRequests[10], MPI_STATUS_IGNORE) );
			

			check ( MPI_Wait(&recvRequests[11], MPI_STATUS_IGNORE) );
			

			check ( MPI_Wait(&recvRequests[12], MPI_STATUS_IGNORE) );
			

			check ( MPI_Wait(&recvRequests[13], MPI_STATUS_IGNORE) );
			

			/* Wait and Update the ZX edges */

			check ( MPI_Wait(&recvRequests[14], MPI_STATUS_IGNORE) );
			

			check ( MPI_Wait(&recvRequests[15], MPI_STATUS_IGNORE) );
			

			check ( MPI_Wait(&recvRequests[16], MPI_STATUS_IGNORE) );
			

			check ( MPI_Wait(&recvRequests[17], MPI_STATUS_IGNORE) );
			
			
			/* Wait for the corners to come in before updating them */
			
			check ( MPI_Wait(&recvRequests[18], MPI_STATUS_IGNORE) );
			check ( MPI_Wait(&recvRequests[19], MPI_STATUS_IGNORE) );
			check ( MPI_Wait(&recvRequests[20], MPI_STATUS_IGNORE) );
			check ( MPI_Wait(&recvRequests[21], MPI_STATUS_IGNORE) );
			check ( MPI_Wait(&recvRequests[22], MPI_STATUS_IGNORE) );
			check ( MPI_Wait(&recvRequests[23], MPI_STATUS_IGNORE) );
			check ( MPI_Wait(&recvRequests[24], MPI_STATUS_IGNORE) );
			check ( MPI_Wait(&recvRequests[25], MPI_STATUS_IGNORE) );
			// if (rank==0){
			// 	printf("rank: %d, before update: %d\n",rank, totalProcessorResults[0][0][0][0]);	
			// }
			
			
			updateCenter(Ix,Iy,Iz, data, newData, rank);
			
			updateXside(Ix, Iy, Iz, data, newData, Ix, rank);
			updateXside(Ix, Iy, Iz, data, newData, 1, rank);
			updateYside(Ix, Iy, Iz, data, newData, Iy, rank);
			updateYside(Ix, Iy, Iz, data, newData, 1, rank);
			updateZside(Ix, Iy, Iz, data, newData, Iz, rank);
			updateZside(Ix, Iy, Iz, data, newData, 1, rank);

			updateXYedge(Ix, Iy, Iz, data, newData, Ix, Iy, rank);
			updateXYedge(Ix, Iy, Iz, data, newData, Ix, 1, rank);
			updateXYedge(Ix, Iy, Iz, data, newData, 1, Iy, rank);
			updateXYedge(Ix, Iy, Iz, data, newData, 1, 1, rank);
			updateYZedge(Ix, Iy, Iz, data, newData, Iy, Iz, rank);
			updateYZedge(Ix, Iy, Iz, data, newData, Iy, 1, rank);
			updateYZedge(Ix, Iy, Iz, data, newData, 1, Iz, rank);
			updateYZedge(Ix, Iy, Iz, data, newData, 1, 1, rank);
			updateZXedge(Ix, Iy, Iz, data, newData, Iz, Ix, rank);
			updateZXedge(Ix, Iy, Iz, data, newData, 1, Ix, rank);
			updateZXedge(Ix, Iy, Iz, data, newData, Iz, 1, rank);
			updateZXedge(Ix, Iy, Iz, data, newData, 1, 1, rank);
			updateCorners(Ix, Iy, Iz, data, newData, rank);

			// check ( MPI_Wait(&sendRequests[0], MPI_STATUS_IGNORE) );
			// check ( MPI_Wait(&sendRequests[1], MPI_STATUS_IGNORE) );
			// check ( MPI_Wait(&sendRequests[2], MPI_STATUS_IGNORE) );
			// check ( MPI_Wait(&sendRequests[3], MPI_STATUS_IGNORE) );
			// check ( MPI_Wait(&sendRequests[4], MPI_STATUS_IGNORE) );
			// check ( MPI_Wait(&sendRequests[5], MPI_STATUS_IGNORE) );
			// check ( MPI_Wait(&sendRequests[6], MPI_STATUS_IGNORE) );
			// check ( MPI_Wait(&sendRequests[7], MPI_STATUS_IGNORE) );
			// check ( MPI_Wait(&sendRequests[8], MPI_STATUS_IGNORE) );
			// check ( MPI_Wait(&sendRequests[9], MPI_STATUS_IGNORE) );
			// check ( MPI_Wait(&sendRequests[10], MPI_STATUS_IGNORE) );
			// check ( MPI_Wait(&sendRequests[11], MPI_STATUS_IGNORE) );
			// check ( MPI_Wait(&sendRequests[12], MPI_STATUS_IGNORE) );
			// check ( MPI_Wait(&sendRequests[13], MPI_STATUS_IGNORE) );
			// check ( MPI_Wait(&sendRequests[14], MPI_STATUS_IGNORE) );
			// check ( MPI_Wait(&sendRequests[15], MPI_STATUS_IGNORE) );
			// check ( MPI_Wait(&sendRequests[16], MPI_STATUS_IGNORE) );
			// check ( MPI_Wait(&sendRequests[17], MPI_STATUS_IGNORE) );
			// check ( MPI_Wait(&sendRequests[18], MPI_STATUS_IGNORE) );
			// check ( MPI_Wait(&sendRequests[19], MPI_STATUS_IGNORE) );
			// check ( MPI_Wait(&sendRequests[20], MPI_STATUS_IGNORE) );
			// check ( MPI_Wait(&sendRequests[21], MPI_STATUS_IGNORE) );
			// check ( MPI_Wait(&sendRequests[22], MPI_STATUS_IGNORE) );
			// check ( MPI_Wait(&sendRequests[23], MPI_STATUS_IGNORE) );
			// check ( MPI_Wait(&sendRequests[24], MPI_STATUS_IGNORE) );
			// check ( MPI_Wait(&sendRequests[25], MPI_STATUS_IGNORE) );
			check ( MPI_Waitall(26, sendRequests,  MPI_STATUSES_IGNORE) );

			
			if (rank == 0) {
				// printf("receiver\n");
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
				// printf("DATA\n");
				// 		for (x=0;x<Ix+2;x++) {
				// for (z=0;z<Iz+2;z++) {
				// 	for (y=0;y<Iy+2;y++) {
						
				// 			// printf("%d, ",XbackwardGhosts[y][z]);
				// 			printf("%d, ",data[x][y][z]);
				// 		}
				// 		printf("\n");
				// 	}
				// 	printf("\n");
				// }
				// printf("NEWDATA\n");
				// for (z=0;z<Iz+2;z++) {
				// 	for (y=0;y<Iy+2;y++) {
				// 		for (x=0;x<Ix+2;x++) {
						
				// 			// printf("%d, ",XbackwardGhosts[y][z]);
				// 			printf("%d, ",newData[x][y][z]);
				// 		}
				// 		printf("\n");
				// 	}
				// 	printf("\n");
				// }
			}



			copyData(Ix,Iy,Iz, data, newData );
			// unsigned char *tmp[Ix][Iy][Iz] = &data[0][0][0];
			// if (rank==5) {
			// 	printf("%d",data[0][0][0]);
			// 	printf("%d",tmp[0][0][0]);
			// 	printf("%d",data[5][0][0]);
			// 	printf("%d",tmp[5][0][0]);
			// 	printf("%d",data[3][0][3]);
			// 	printf("%d",tmp[3][0][3]);
			// }
			//data = newData;
			//newData = tmp;
			saveData(Ix,Iy,Iz, newData, i, totalProcessorResults );

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

		}
	}

	double end = MPI_Wtime();

	printf("P: %d\n(Px, Py, Pz): (%d, %d, %d)\nDatasize (x,y,z): (%d,%d,%d) => %d\nRunTime: %e\n", size, Px, Py, Pz, xN, yN, zN, xN*yN*zN, end-start);


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
	// if (rank==0){
	// 	printf("%d %d %d\n", dx, dy,dz);
	// 	printf("the element: %d", totalProcessorResults[0][0][0][0]);
	// }
	printf("THE END\n");
	for (t=0;t<=iterations;t++) {
		for (x=0;x<Ix;x++) {
			for (y=0;y<Iy;y++) {
				for (z=0;z<Iz;z++) {
					//fprintf(fp, "%d ", totalProcessorResults[x][y][z][t]);
					if ( totalProcessorResults[x][y][z][t] == 1 ) {
						//printf("%d %d %d ", x+dx, y+dy, z+dz);
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














































