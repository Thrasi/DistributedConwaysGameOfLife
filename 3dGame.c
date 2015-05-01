/* Use MPI */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MIN(a,b) ((a) < (b) ? (a) : (b))

#define iterations 2 /* number of iterations */

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

/* Check function. */
void check(int rc) {
	if (rc != MPI_SUCCESS) {
		exit(-1);
	}
}

/* Count the number of active cells around position (x,y,z) */
unsigned char count(int ox, int oy, int oz, unsigned char data[ox][oy][oz], int x, int y, int z) {
	unsigned char c = 0;
    int dx,dy,dz;
    for (dx=-1;dx<=1;dx++) {
    	for (dy=-1;dy<=1;dy++) {
    		for (dy=-1;dy<=1;dy++) {
    			c += data[x+dx][y+dy][z+dz];

    		}
    	}	
    }
    //c -= data[x][y][z];
	return c;
}

void copyData(int ox, int oy, int oz, unsigned char data[ox][oy][oz],
				  int Ix, int Iy, int Iz, unsigned char newData[Ix][Iy][Iz], 
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

/* Given an old cell value and a count of neighbours calculate the new one */
unsigned char transition(unsigned char c, unsigned char old) {
	if ( old == 1) {
		if ( c == 4 || c == 5 ) {
			return 1;
		} else {
			return 0;
		}
	} else if ( c == 4 ) {
	        return 1;
	} else {
	        return 0;
	}
}

void updateCenter(int ox, int oy, int oz, unsigned char data[ox][oy][oz],
				  int Ix, int Iy, int Iz, unsigned char newData[Ix][Iy][Iz], int rank) {

	int x, y, z;
	unsigned char c, old;
	int xLim = Ix-1;
	int yLim = Iy-1;
	int zLim = Iz-1;

	for (x=1;x<xLim;x++) {
		for (y=1;y<yLim;y++) {
			for (z=1;z<zLim;z++) {
				c = count(ox, oy, oz, data, x+1, y+1, z+1);
				if (rank==0){
					printf("%d,%d,%d: %hhu\n",x+1,y+1,z+1,c);
				}
				old = data[x+1][y+1][z+1];
				newData[x][y][z] = transition(c, old);
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
	int X = 0, Y = 1, Z = 2, FORWARD = 1, BACKWARDS = -1;
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
	unsigned char *picture;
	unsigned char *picData;

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

		data[2][2][2] = 1;
		data[3][3][2] = 1;
		data[2][3][3] = 1;
		data[3][2][3] = 1;
	}

	for (x=0;x<Ix;x++) {
		for (y=0;y<Iy;y++) {
			for (z=0;z<Iz;z++) {
				totalProcessorResults[x][y][z][0] = data[x+1][y+1][z+1];
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

	for (i=1;i<=iterations;i++) {
		//updateBuffers();
		if (doMPI==1) {
			// if (rank==1 || rank==5) {
			// printf("X: Rank %d: (%d,%d,%d) sends to %d, receives from %d\n",rank,coord[0],coord[1],coord[2],Xforward, Xbackward);
			// printf("Y: Rank %d: (%d,%d,%d) sends to %d, receives from %d\n",rank,coord[0],coord[1],coord[2],Yforward, Ybackward);
			// printf("Z: Rank %d: (%d,%d,%d) sends to %d, receives from %d\n",rank,coord[0],coord[1],coord[2],Zforward, Zbackward);

   //          MPI_Isend(&XforwardBuffers, 1, Xside, Xforward, 1, TORUS_COMM, &rXforward);
			// MPI_Irecv(&XbackwardGhosts, 1, Xside, Xbackward, 1, TORUS_COMM, &rXbackward);

			// MPI_Isend(&YforwardBuffers, 1, Yside, Yforward, 2, TORUS_COMM, &rYforward);
			// MPI_Irecv(&YbackwardGhosts, 1, Yside, Ybackward, 2, TORUS_COMM, &rYbackward);

			// MPI_Isend(&ZforwardBuffers, 1, Zside, Zforward, 3, TORUS_COMM, &rZforward);
			// MPI_Irecv(&ZbackwardGhosts, 1, Zside, Zbackward, 3, TORUS_COMM, &rZbackward);

			updateCenter(Ix+2,Iy+2,Iz+2, data, Ix,Iy,Iz, newData, rank);

			// MPI_Wait(&rXforward, MPI_STATUS_IGNORE);
			// MPI_Wait(&rXbackward, MPI_STATUS_IGNORE);
			// MPI_Wait(&rYforward, MPI_STATUS_IGNORE);
			// MPI_Wait(&rYbackward, MPI_STATUS_IGNORE);
			// MPI_Wait(&rZforward, MPI_STATUS_IGNORE);
			// MPI_Wait(&rZbackward, MPI_STATUS_IGNORE);

			printf("Rank %d: (%d,%d,%d)  got through 1X\n",rank,coord[0],coord[1],coord[2]);

			/*sendRecvEdges() // XY
			sendRecv() // YZ
			sendRecv() // ZX

			sendRecv() // corners*/

			copyData(Ix+2,Iy+2,Iz+2, data, 
					 Ix,Iy,Iz, newData,
					 i, totalProcessorResults );
		}
	}

	if (rank==0) {
	printf("\n");
	FILE *fp;
	//	 Save the file. 
	char str[15];
	sprintf(str, "%ddata.txt",rank);
	fp = fopen(str, "w");
	int t;
	fprintf(fp, "%d %d %d \n",Ix,Iy,Iz);
	for (t=0;t<=iterations;t++) {
		for (x=0;x<Ix;x++) {
			for (y=0;y<Iy;y++) {
				for (z=0;z<Iz;z++) {
					fprintf(fp, "%d ", totalProcessorResults[x][y][z][t]);
				}
			}
		}
		fprintf(fp, "\n");
	}

	fclose(fp);
	
	}
	MPI_Finalize();
	exit(0);


}














































