/* Use MPI */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

/* Check function. */
void check(int rc) {
	if (rc != MPI_SUCCESS) {
		exit(-1);
	}
}

/* This calls sendRecv to communicate the edge ghost values */
void sendRecv(void *sendbuf, void *recvbuf, int sendcount, int dest, int source, MPI_Comm TORUS_COMM) {
	check ( MPI_Sendrecv(sendbuf, sendcount, MPI_UNSIGNED_CHAR, dest, 0,
                			 recvbuf, sendcount, MPI_UNSIGNED_CHAR, source, 0,
                			 TORUS_COMM, MPI_STATUS_IGNORE) );
}

void sendRecvSide(void *sendbuf, void *recvbuf, MPI_Datatype type, int dest, int source, MPI_Comm TORUS_COMM)  {
	int sendCount = 1;
	int tag1 = 0;
	int tag2 = 1;
	check (
		 MPI_Sendrecv(sendbuf, sendCount, type, dest, tag1,
                			 recvbuf, sendCount, type, source, tag2,
                			 TORUS_COMM, MPI_STATUS_IGNORE)

		);
}

/* Count the number of active cells around position (x,y,z) */
int count(int ox, int oy, int oz, unsigned char data[ox][oy][oz], int x, int y, int z) {
	unsigned char c = 0;
    int dx,dy,dz;
    for (dx=-1;dx<=1;dx++) {
    	for (dy=-1;dy<=1;dy++) {
    		for (dy=-1;dy<=1;dy++) {
    			c += data[x+dx][y+dy][z+dz];
    		}
    	}	
    }
    c -= data[x][y][z];
	return c;
}

/* Given an old cell value and a count of neighbours calculate the new one */
unsigned char transition(unsigned char c, unsigned char old) {
	if ( old == 1) {
		if ( c == 4 || c == 5 ) {
			return 1;
		} else {
			return 0;
		}
	} else if ( c == 5 ) {
	        return 1;
	} else {
	        return 0;
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
	unsigned char YZnorthBuffers[4][Ix];
	memset(YZnorthBuffers, 0, sizeof YZnorthBuffers);
	unsigned char ZXnorthBuffers[4][Iy];
	memset(ZXnorthBuffers, 0, sizeof ZXnorthBuffers);

	/* corners are enumerated in binary as their coordinates [x, y, y]  */
	unsigned char cornerBuffers[8];
	memset(cornerBuffers, 0, sizeof cornerBuffers);

	//printf("generate data\n");
	// inital image
	srandom(rank+1);
	int x,y,z;
	for (y=0;y<Iy;y++) {
		for (x=0;x<Ix;x++) {
			for (z=0;z<Iz;z++) {
				data[x][y][z] = random() % 2;
			}
		}
	}
	//printf("side ranks\n");
	/* Side ranks */
	check ( MPI_Cart_shift(TORUS_COMM, X, FORWARD, &Xforward, &Xbackward) );
	check ( MPI_Cart_shift(TORUS_COMM, Y, FORWARD, &Yforward, &Ybackward) );
	check ( MPI_Cart_shift(TORUS_COMM, Z, FORWARD, &Zforward, &Zbackward) );

	unsigned char *totalProcessorResults = (unsigned char*) malloc(Ix*Iy*Iz*(iterations+1)*sizeof(unsigned char));
	//printf("createn data type\n");
	/* Define new datatypes for communications */
	MPI_Datatype  Xside, Yside, Zside; 
	MPI_Type_contiguous(Iy*Iz, MPI_UNSIGNED_CHAR, &Xside);
	MPI_Type_contiguous(Iz*Ix, MPI_UNSIGNED_CHAR, &Yside);
	MPI_Type_contiguous(Ix*Iy, MPI_UNSIGNED_CHAR, &Zside);
	MPI_Type_commit(&Xside);
	MPI_Type_commit(&Yside);
	MPI_Type_commit(&Zside);

	MPI_Datatype tmp;
	MPI_Type_contiguous(4, MPI_UNSIGNED_CHAR, &tmp);
	MPI_Type_commit(&tmp);

	int a=1;
	int b=2;
	unsigned char c[] = {1,1};
	unsigned char d[] = {2,2};
	unsigned char e[2][2];
	unsigned char f[2][2];
	memset(e, 0, sizeof e);
	memset(f, 0, sizeof f);
	MPI_Request r1, r2;
	for (i=0;i<iterations+1;i++) {
		//updateBuffers();
		if (doMPI==1) {
			if (rank==1 || rank==5) {
			printf("Rank %d: (%d,%d,%d) sends to %d, receives from %d\n",rank,coord[0],coord[1],coord[2],Xforward, Xbackward);
			//MPI_SendRecv(&a, &b, MPI_INT, Xforward, Xbackward, TORUS_COMM); // X
			//sendRecvSide(XforwardBuffers, XbackwardGhosts, Xside, Xforward, Xbackward, TORUS_COMM); // X
			/*check ( MPI_Sendrecv(&a, 1, MPI_INT, Xforward, 0,
                			 &b, 1, MPI_INT, Xbackward, 0,
                			 TORUS_COMM, MPI_STATUS_IGNORE) );*/
			// MPI_Isend(&a, 1, MPI_INT, Xforward, 11, MPI_COMM_WORLD, &r1);
			// MPI_Irecv(&b, 1, MPI_INT, Xbackward, 11, MPI_COMM_WORLD, &r2);

			// MPI_Isend(&c, 1, tmp, Xforward, 11, MPI_COMM_WORLD, &r1);
			// MPI_Irecv(&d, 1, tmp, Xbackward, 11, MPI_COMM_WORLD, &r2);

			// MPI_Isend(&e, 1, tmp, Xforward, 11, MPI_COMM_WORLD, &r1);
			// MPI_Irecv(&f, 1, tmp, Xbackward, 11, MPI_COMM_WORLD, &r2);

            MPI_Isend(&XforwardBuffers, 1, Xside, Xforward, 11, TORUS_COMM, &r1);
			MPI_Irecv(&XbackwardGhosts, 1, Xside, Xbackward, 11, TORUS_COMM, &r2);

			MPI_Wait(&r1,MPI_STATUS_IGNORE);
			MPI_Wait(&r2,MPI_STATUS_IGNORE);

			printf("Rank %d: (%d,%d,%d)  got through 1X\n",rank,coord[0],coord[1],coord[2]);
			}
			/*sendRecvSide(XbackwardBuffers, XforwardGhosts, Xside, Xbackward, Xforward, TORUS_COMM); // X
			printf("Rank %d: (%d,%d)  got through 2X\n",rank);
			sendRecvSide(YforwardBuffers, YbackwardGhosts, Yside, Yforward, Ybackward, TORUS_COMM); // Y
			printf("Rank %d: (%d,%d)  got through 1Y\n",rank);
			sendRecvSide(YbackwardBuffers, YforwardGhosts, Yside, Ybackward, Yforward, TORUS_COMM); // Y
			printf("Rank %d: (%d,%d)  got through 2Y\n",rank);
			sendRecvSide(ZforwardBuffers, ZbackwardGhosts, Zside, Zforward, Zbackward, TORUS_COMM); // Z
			printf("Rank %d: (%d,%d)  got through 1Z\n",rank);
			sendRecvSide(ZbackwardBuffers, ZforwardGhosts, Zside, Zbackward, Zforward, TORUS_COMM); // Z
			printf("Rank %d: (%d,%d)  got through 2Z\n",rank);*/

			/*sendRecvEdges() // XY
			sendRecv() // YZ
			sendRecv() // ZX

			sendRecv() // corners*/

			/* copy buffers into data */

		}
	}

	MPI_Finalize();
	exit(0);


}

//void copySideToData(unsigned char )
void updateCenter(int ox, int oy, int oz, unsigned char data[ox][oy][oz],
				  int Ix, int Iy, int Iz, unsigned char newData[Ix][Iy][Iz]) {

	int x, y, z;
	unsigned char c, old;
	int xLim = Ix-1;
	int yLim = Iy-1;
	int zLim = Iz-1;

	for (x=1;x<xLim;x++) {
		for (y=1;y<yLim;y++) {
			for (z=1;z<zLim;z++) {
				c = count(ox, oy, oz, data, x+1, y+1, z+1);
				old = data[x+1][y+1][z+1];
				newData[x][y][z] = transition(c, old);
			}
		}
	}
}













































