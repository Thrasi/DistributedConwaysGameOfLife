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

void copyData(int Ix, int Iy, int Iz, unsigned char data[Ix][Iy][Iz],
				   unsigned char newData[Ix][Iy][Iz], 
				  int t, unsigned char totalProcessorResults[Ix][Iy][Iz][iterations+1]) {
	int x, y, z;
	for (x=0;x<Ix;x++) {
		for (y=0;y<Iy;y++) {
			for (z=0;z<Iz;z++) {
				data[x][y][z] = newData[x][y][z];
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

	unsigned char data[Ix][Iy][Iz];
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

		data[2][2][2] = 1;
		data[3][3][2] = 1;
		data[2][3][3] = 1;
		data[3][2][3] = 1;
	}
	if (rank == 1) {
		XforwardBuffers[1][1] = 1;
		XforwardBuffers[1][2] = 1;
		XforwardBuffers[1][3] = 1;
		XforwardBuffers[2][2] = 1;
	}


	for (x=0;x<Ix;x++) {
		for (y=0;y<Iy;y++) {
			for (z=0;z<Iz;z++) {
				totalProcessorResults[x][y][z][0] = data[x][y][z];
			}
		}
	}

	//printf("side ranks\n");
	/* Side ranks */
	check ( MPI_Cart_shift(TORUS_COMM, X, FORWARD, &Xforward, &Xbackward) );
	check ( MPI_Cart_shift(TORUS_COMM, Y, FORWARD, &Yforward, &Ybackward) );
	check ( MPI_Cart_shift(TORUS_COMM, Z, FORWARD, &Zforward, &Zbackward) );

	
	//printf("createn data type\n");
	/* Define new datatypes for communications */
	MPI_Datatype  Xside, Yside, Zside; 
	MPI_Type_contiguous(Iy*Iz, MPI_UNSIGNED_CHAR, &Xside);
	MPI_Type_contiguous(Iz*Ix, MPI_UNSIGNED_CHAR, &Yside);
	MPI_Type_contiguous(Ix*Iy, MPI_UNSIGNED_CHAR, &Zside);
	MPI_Type_commit(&Xside);
	MPI_Type_commit(&Yside);
	MPI_Type_commit(&Zside);

	MPI_Request rXforward, rXbackward, rYforward, rYbackward, rZforward, rZbackward;
	if (rank == 1) {
		for (y=0;y<Iy;y++) {
			for (z=0;z<Iz;z++) {
				printf("%d, ",XforwardBuffers[y][z]);
			}
			printf("\n");
		}
	}
	printf("\n");
	for (i=1;i<=iterations;i++) {
		//updateBuffers();
		if (doMPI==1) {
			// if (rank==1 || rank==5) {
			// printf("X: Rank %d: (%d,%d,%d) sends to %d, receives from %d\n",rank,coord[0],coord[1],coord[2],Xforward, Xbackward);
			// printf("Y: Rank %d: (%d,%d,%d) sends to %d, receives from %d\n",rank,coord[0],coord[1],coord[2],Yforward, Ybackward);
			// printf("Z: Rank %d: (%d,%d,%d) sends to %d, receives from %d\n",rank,coord[0],coord[1],coord[2],Zforward, Zbackward);

			MPI_Isend(&XforwardBuffers, 1, Xside, Xforward, 1, TORUS_COMM, &rXforward);
			MPI_Irecv(&XbackwardGhosts, 1, Xside, Xbackward, 1, TORUS_COMM, &rXbackward);

			// MPI_Isend(&YforwardBuffers, 1, Yside, Yforward, 2, TORUS_COMM, &rYforward);
			// MPI_Irecv(&YbackwardGhosts, 1, Yside, Ybackward, 2, TORUS_COMM, &rYbackward);

			// MPI_Isend(&ZforwardBuffers, 1, Zside, Zforward, 3, TORUS_COMM, &rZforward);
			// MPI_Irecv(&ZbackwardGhosts, 1, Zside, Zbackward, 3, TORUS_COMM, &rZbackward);
			// if (rank==0) {
			// 	printf("step %d",i);
			// }
			updateCenter(Ix,Iy,Iz, data, newData, rank);

			MPI_Wait(&rXforward, MPI_STATUS_IGNORE);
			MPI_Wait(&rXbackward, MPI_STATUS_IGNORE);

			//updateXside(Ix, Iy, Iz, XbackwardGhosts, data, newData, BACKWARDS, rank);
			
			if (rank == 5) {
				updateXside(Ix, Iy, Iz, XbackwardGhosts, data, newData, BACKWARDS, rank);

				for (y=0;y<Iy;y++) {
					for (z=0;z<Iz;z++) {
						printf("%d, ",XbackwardGhosts[y][z]);
					}
					printf("\n");
				}
			}
			// MPI_Wait(&rYforward, MPI_STATUS_IGNORE);
			// MPI_Wait(&rYbackward, MPI_STATUS_IGNORE);
			// MPI_Wait(&rZforward, MPI_STATUS_IGNORE);
			// MPI_Wait(&rZbackward, MPI_STATUS_IGNORE);

			//printf("Rank %d: (%d,%d,%d)  got through 1X\n",rank,coord[0],coord[1],coord[2]);

			/*sendRecvEdges() // XY
			sendRecv() // YZ
			sendRecv() // ZX

			sendRecv() // corners*/

			copyData(Ix,Iy,Iz, data, newData,
					 i, totalProcessorResults );
			if (rank == 5) {
				printf("newData\n");
				for (y=0;y<Iy;y++) {
					for (z=0;z<Iz;z++) {
						printf("%d, ",newData[0][y][z]);
					}
					printf("\n");
				}
				printf("data\n");
				for (y=0;y<Iy;y++) {
					for (z=0;z<Iz;z++) {
						printf("%d, ",data[0][y][z]);
					}
					printf("\n");
				}
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














































