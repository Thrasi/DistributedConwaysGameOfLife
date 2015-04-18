/* Use MPI */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MIN(a,b) ((a) < (b) ? (a) : (b))

#define iterations 400 /* number of iterations */

#define xDim 2
#define yDim 2
#define xBC 1
#define yBC 1
#define N 4096
#define sqrtN 64


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


			
/* Calls sendrecv to communicate the corner cells */
void sendRecvCorners(void *sendbuf, void *recvbuf, int sendcount,int *coord, int dx, int dy,  MPI_Comm TORUS_COMM) {
	
	int dest, source;
	int destCoord[2] = {coord[0]+dx, coord[1]+dy};
	int sourceCoord[2] = {coord[0]-dx, coord[1]-dy};

	check ( MPI_Cart_rank(TORUS_COMM, destCoord, &dest) );
	check ( MPI_Cart_rank(TORUS_COMM, sourceCoord, &source) );

	//printf("(%d,%d), NW: %d, SE: %d",coord[0],coord[1], dest, source);
	
	check ( MPI_Sendrecv(sendbuf, sendcount, MPI_UNSIGNED_CHAR, dest, 0,
                			 recvbuf, sendcount, MPI_UNSIGNED_CHAR, source, 0,
                			 TORUS_COMM, MPI_STATUS_IGNORE) );
}

/*Count the number of active cells around position (i,j)*/
int count(int x, int y, int X, int Y,  unsigned char *data) {
	unsigned char c = 0;
    int y1 = y-1, y2 = y+1;
    int x1 = x-1, x2 = x+1;
    
	c += data[y1*X + x1];
	c += data[y1*X + x];
	c += data[y1*X + x2];

	c += data[y*X + x1];
	c += data[y*X + x2];

	c += data[y2*X + x1];
	c += data[y2*X + x];
	c += data[y2*X + x2];

	return c;
}

/* Given an old cell value and a count of neighbours calculate the new one*/
unsigned char transition(unsigned char c, unsigned char old) {
	unsigned char ret;
	if ( old == 1) {
        	if ( c < 2 ) {
                	ret = 0;
	        } else if ( c <= 3 ) {
        	        ret = 1;
        	} else {
        	        ret = 0;
        	}
	} else if ( c == 3 ) {
	        ret = 1;
	} else {
	        ret = 0;
	}
	return ret;
}

/*update the send buffers for side communications*/
void updateBuffers(unsigned char *sendBufEast, unsigned char *sendBufWest, unsigned char *sendBufNorth, 
			unsigned char *sendBufSouth, unsigned char *data, int I) {
	int i;
	int II = I*I;
	for (i=0;i<I;i++) {
		sendBufEast[i] = data[i*I + (I-1)];
		sendBufWest[i] = data[i*I];
		sendBufNorth[i] = data[i];
		sendBufSouth[i] = data[II-I+i];
	}
}

/* calculate the new values of the center cells*/
void updateCenter(int X, int Y, unsigned char *newData, unsigned char *data) {
	int x, y;
	unsigned char c, old;
	for (y=1;y<Y-1;y++) {
	    for (x=1;x<X-1;x++) {
			c = count(x, y, X, Y, data);
			old = data[y*X + x];
			newData[y*X + x] = transition(c, old);
		}
	}
}

/* calculate the new values on the vertical edges */
void updateSide (int Y, unsigned char *ghosts, int col, int inner, unsigned char *newData, unsigned char *data) {
	int y, y1, y2;
	int X = Y;
    
	unsigned char c, old;
	for (y=1;y<Y-1;y++) {
	    y1 = y-1, y2 = y+1;
		c = 0;
		c += ghosts[y1];
		c += ghosts[y];
		c += ghosts[y2];

		c += data[y1*X + col];
		c += data[y2*X + col];

		c += data[y1*X + inner];
		c += data[y*X + inner];
		c += data[y2*X + inner];

		old = data[y*X + col];
		newData[y*X + col] = transition(c, old);
	}
}

/* calculate the new values on the horizontal edges */
void updateRow (int X, unsigned char *ghosts, int row, int inner, unsigned char *newData, unsigned char *data) {
	int x, x1, x2;
	unsigned char c, old;
	for (x=1;x<X-1;x++) {
	    x1 = x-1, x2 = x+1;
        c = 0;
        c += ghosts[x1];
        c += ghosts[x];
        c += ghosts[x2];

        c += data[row*X + x1];
        c += data[row*X + x2];

        c += data[inner*X + x1];
        c += data[inner*X + x];
        c += data[inner*X + x2];

		old = data[row*X + x];
        newData[row*X + x] = transition(c, old);
    }
}
/*This updates the corner cells.*/
void updateCorners ( int Ix, int Iy, unsigned char *newData, unsigned char *data, 
					unsigned char *eastGhosts, unsigned char *northGhosts,
					unsigned char *westGhosts, unsigned char *southGhosts,  
					unsigned char northWest, unsigned char southEast, 
					unsigned char southWest, unsigned char northEast ) {
	
	unsigned char c = 0, old;
	// NW:
	c += northWest;
	c += northGhosts[0] + northGhosts[1];
	c += westGhosts[0] + westGhosts[1];
	c += data[1] + data[Ix] + data[Ix+1];
	old = data[0];
	newData[0] = transition(c, old);

	// NE:
	c = 0;
	c += northEast;
	c += northGhosts[Ix-1] + northGhosts[Ix-2];
	c += eastGhosts[0] + eastGhosts[1];
	c += data[Ix-2] + data[2*Ix-1] + data[2*Ix-2];
	old = data[Ix-1];
	newData[Ix-1] = transition(c, old);

	// SW:
	c = 0;
	c += southWest;
	c += westGhosts[Iy-1] + westGhosts[Iy-2];
	c += southGhosts[0] + southGhosts[1];
	c += data[(Ix-2)*Iy] + data[(Ix-2)*Iy+1] + data[(Ix-1)*Iy+1];
	old = data[(Ix-1)*Iy];
	newData[(Ix-1)*Iy] = transition(c, old);

	// SE:
	c = 0;
	c += southEast;
	c += eastGhosts[Iy-1] + eastGhosts[Iy-2];
	c += southGhosts[Ix-1] + southGhosts[Ix-2];
	c += data[(Ix-1)*Iy-1] + data[(Ix-1)*Iy-2] + data[Ix*Iy-2];
	old = data[Ix*Iy-1];
	newData[Ix*Iy-1] = transition(c, old);
}


int main(int argc, char *argv[]){
	int rank, size, i, EAST, WEST, NORTH, SOUTH, Nx, Lx, Rx, Ix, Ny, Ly, Ry, Iy;
	int VERTICAL = 0, HORIZONTAL = 1, FORWARD = 1, BACKWARDS = -1;
    int dim[2], period[2], reorder;
    int coord[2], id;
    int doMPI = 0;
    MPI_Comm TORUS_COMM;

    /* initialize MPI and create virtual 2d torus topology*/

    check ( MPI_Init(&argc, &argv) );
    check ( MPI_Comm_rank(MPI_COMM_WORLD, &rank) );
    check ( MPI_Comm_size(MPI_COMM_WORLD, &size) );
//    printf("Rank %d:  Initialized mpi\n",rank);

	dim[0] = xDim; dim[1] = yDim;
	period[0] = xBC; period[1] = yBC;
	reorder = 0;

	check ( MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &TORUS_COMM) );
	check ( MPI_Cart_coords(TORUS_COMM, rank, 2, coord) );
//	printf("Rank %d:  cartesian created\n",rank);

    /* 2d linear data distribution */
	Nx = sqrtN;
	Lx = Nx / xDim;
	Rx = Nx % xDim;
	Ix = ( Nx + xDim - coord[0]-1 ) / xDim;

	Ny = sqrtN;
	Ly = Ny / yDim;
	Ry = Ny % yDim;
	Iy = ( Ny + yDim - coord[0]-1 ) / yDim;
//    printf("Rank %d:  Ix: %d, Iy: %d\n", rank, Ix, Iy);
    
	// declare independence
	unsigned char *picture, *data, *newData, *northGhosts, *southGhosts, *eastGhosts, *westGhosts;
	unsigned char northWest=0, southWest=0, northEast=0, southEast=0;
	unsigned char *sendBufSouth, *sendBufNorth, *sendBufEast, *sendBufWest;
	unsigned char *picData;
	
	unsigned char *totalProcessorResults;

	// initialize data
	data = (unsigned char*) calloc(Ix*Iy,sizeof(unsigned char));
	newData = (unsigned char*) calloc(Ix*Iy,sizeof(unsigned char));
	totalProcessorResults = (unsigned char*) malloc(Ix*Iy*(iterations+1)*sizeof(unsigned char));
	if ( rank == 0) {
		printf("size of picture: %d",N*iterations);
		picture = (unsigned char*) malloc(N*iterations*sizeof(unsigned char));
	}
	if ( rank == 0) {
		data[(Iy-10)*Ix + 5] = 1;
		data[(Iy-10)*Ix + 6] = 1;
		data[(Iy-10)*Ix + 7] = 1;
		data[(Iy-11)*Ix + 7] = 1;
		data[(Iy-12)*Ix + 6] = 1;
		//data[11] = 3;
	}
/*	if ( rank == 2) {

		data[3] = 1;
		//data[11] = 3;
	}
	if ( rank == 3) {
		data[0] = 1; data[1] = 1;
		//data[11] = 3;

	}*/
	// ghost points
	northGhosts = (unsigned char*) calloc(Ix, sizeof(unsigned char));
	southGhosts = (unsigned char*) calloc(Ix, sizeof(unsigned char));
	eastGhosts = (unsigned char*) calloc(Iy, sizeof(unsigned char));
	westGhosts = (unsigned char*) calloc(Iy, sizeof(unsigned char));
	
	sendBufNorth = (unsigned char*) calloc(Ix, sizeof(unsigned char));
	sendBufSouth = (unsigned char*) calloc(Ix, sizeof(unsigned char));
	sendBufEast = (unsigned char*) calloc(Iy, sizeof(unsigned char));
	sendBufWest = (unsigned char*) calloc(Iy, sizeof(unsigned char));

	// inital image
	srandom(rank+1);
	int x,y;
	for (y=0;y<Iy;y++) {
		for (x=0;x<Ix;x++) {
			//data[y*Ix + x] = random() % 2;
		}
	}

	check ( MPI_Cart_shift(TORUS_COMM, HORIZONTAL, FORWARD, &EAST, &WEST) );
	check ( MPI_Cart_shift(TORUS_COMM, VERTICAL, FORWARD, &SOUTH, &NORTH) );
	//printf("Rank %d:  Cart shifted\n",rank);
	memcpy(newData, data, Ix*Iy);
	//printf("Rank %d: before loop\n",rank);
	for (i=0;i<iterations;i++) {
		// update the buffers with the data to send
		updateBuffers(sendBufEast, sendBufWest, sendBufNorth, sendBufSouth, data, Ix);
		

		if ( doMPI == 0 ) {
			// EAST-WEST
//			printf("did mpi\n");
			sendRecv ( sendBufEast, westGhosts, Iy, EAST, WEST, TORUS_COMM );
    			// MAYBE HAVE SEPARATE BUFFERS FOR OPTIMIZATION
    			sendRecv ( sendBufWest, eastGhosts, Iy, WEST, EAST, TORUS_COMM );

			// NORTH-SOUTH
			sendRecv ( sendBufNorth, southGhosts, Ix, NORTH, SOUTH, TORUS_COMM );
			sendRecv ( sendBufSouth, northGhosts, Ix, SOUTH, NORTH, TORUS_COMM );

			// send corners
			// NW
			int dx=-1, dy=-1;
			sendRecvCorners ( &data[0], &southEast, 1, coord, dx, dy, TORUS_COMM);
			// SE
			dx=1, dy=1;
			sendRecvCorners ( &data[Ix*Iy-1], &northWest, 1, coord, dx, dy, TORUS_COMM);
			// SW
			dx=-1, dy=1;
			sendRecvCorners ( &data[Ix*Iy-Ix], &northEast, 1, coord, dx, dy, TORUS_COMM);
			// NE
			dx=1, dy=-1;
			sendRecvCorners ( &data[Ix-1], &southWest, 1, coord, dx, dy, TORUS_COMM);
		}

		// Update the grid
		updateCenter ( Ix, Iy, newData, data );
		
		// here we should wait for communications to end.
		updateSide ( Iy, eastGhosts, Iy-1, Iy-2, newData, data );
		updateRow ( Ix, northGhosts, 0, 1, newData, data );
		updateSide ( Iy, westGhosts, 0, 1, newData, data );
		updateRow ( Ix, southGhosts, Ix-1, Ix-2, newData, data );

		updateCorners ( Ix, Iy, newData, data, eastGhosts, northGhosts, westGhosts, southGhosts,  
					northWest, southEast, southWest, northEast );

		//unsigned char sum = 0;
		for (y=0;y<Iy;y++) {
			for (x=0;x<Ix;x++) {
				totalProcessorResults[Ix*Iy*i + y*Ix + x] = data[y*Ix + x];
			//	data[y*Ix + x] = newData[y*Ix + x];
			}
			if (rank==1) {
			//	printf("\n");
			}
		}
		memcpy(data, newData, Ix*Iy);
    	// copy the data into a 1d array for sending
		
		
		//printf("Rank %d, iteration: %d, sum: %d\n",rank, i, sum);
	}
	if (rank==0) {
		//printf("Rank %d has N:%d, S:%d, W:%d, E:%d\n",rank, NORTH,SOUTH,WEST,EAST);
		//printf("Ix*Iy*iterations %d, iteration: %d, sum: %d\n",rank, i, sum);
		for (x=0;x<Ix*Iy*iterations;x++) {
		//	printf("%hhu ", totalProcessorResults[x]);
		}
		//printf("\n");

	}
	check(	MPI_Gather(totalProcessorResults, Ix*Iy*iterations, MPI_UNSIGNED_CHAR, picture,
				 Ix*Iy*iterations, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD) );

	if (rank == 0) {
		printf("\n");
		FILE *fp;
		//	 Save the file. 
		char str[15];
		sprintf(str, "fullData.txt");
		fp = fopen(str, "w");
		int r;
		printf("sqrtN: %d, Ix*Iy*iterations: %d\n",sqrtN,Ix*Iy*iterations);
		for (r=0;r<size;r++) {
			for (x=0;x<Ix*Iy*iterations;x++) {
				fprintf(fp, "%hhu ", picture[Ix*Iy*iterations*r + x]);
				//printf("%hhu ", picture[Ix*Iy*iterations*r + x]);
			}
			fprintf(fp, "\n");
			//printf("\n");

		}
		fclose(fp);
	}

	check( MPI_Finalize() );
	//printf("Rank %d: finished\n",rank);
	exit(0);
}
