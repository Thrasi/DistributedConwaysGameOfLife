
#include "updateFunctions.h"


/* Count the number of active cells around position (x,y,z) */
unsigned char count(int Ix, int Iy, int Iz, unsigned char data[Ix+2][Iy+2][Iz+2], int x, int y, int z) {
	unsigned char c = 0;
    int dx,dy,dz;
    for (dx=-1;dx<=1;dx++) {
    	for (dy=-1;dy<=1;dy++) {
    		for (dz=-1;dz<=1;dz++) {
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
	} else if ( c == 4 ) {
	        return 1;
	} else {
	        return 0;
	}
}

/* calculate the new values of the center cells*/
void updateCenter(int Ix, int Iy, int Iz, unsigned char data[Ix+2][Iy+2][Iz+2],
				   unsigned char newData[Ix+2][Iy+2][Iz+2], int rank) {

	int x, y, z;
	unsigned char c, old;
	int xLim = Ix+1;
	int yLim = Iy+1;
	int zLim = Iz+1;

	for (x=1;x<xLim;x++) {
		for (y=1;y<yLim;y++) {
			for (z=1;z<zLim;z++) {
				c = count(Ix, Iy, Iz, data, x, y, z);
				old = data[x][y][z];
				newData[x][y][z] = transition(c, old);
			}
		}
	}
}

/* calculate the new values of a side with index X */
void updateXside(int Ix, int Iy, int Iz, unsigned char data[Ix+2][Iy+2][Iz+2],
				 unsigned char newData[Ix+2][Iy+2][Iz+2], int X, int rank) {
	int x, y, z;
	x = X+1;
	unsigned char c, old;
	int yLim = Iy+1;
	int zLim = Iz+1;

	for (y=1;y<yLim;y++) {
		for (z=1;z<zLim;z++) {
			c = count(Ix, Iy, Iz, data, x, y, z);
			old = data[x][y][z];
			newData[X][y][z] = transition(c, old);
		}
	}
}

/* calculate the new values of a side with index Y */
void updateYside(int Ix, int Iy, int Iz, unsigned char data[Ix+2][Iy+2][Iz+2],
				 unsigned char newData[Ix+2][Iy+2][Iz+2], int Y, int rank) {
	int x, y, z;
	y = Y+1;
	unsigned char c, old;
	int xLim = Ix+1;
	int zLim = Iz+1;

	for (x=1;x<xLim;x++) {
		for (z=1;z<zLim;z++) {
			c = count(Ix, Iy, Iz, data, x, y, z);
			old = data[x][y][z];
			newData[x][Y][z] = transition(c, old);
		}
	}
}

/* calculate the new values of a side with index Z */
void updateZside(int Ix, int Iy, int Iz, unsigned char data[Ix+2][Iy+2][Iz+2],
				 unsigned char newData[Ix+2][Iy+2][Iz+2], int Z, int rank) {
	int x, y, z;
	z = Z+1;
	unsigned char c, old;
	int yLim = Iy+1;
	int xLim = Ix+1;

	for (y=1;y<yLim;y++) {
		for (x=1;x<xLim;x++) {
			c = count(Ix, Iy, Iz, data, x, y, z);
			old = data[x][y][z];
			newData[x][y][Z] = transition(c, old);
		}
	}
}

/* Calculate the new values of an edge normal to the XY plane */
void updateXYedge(int Ix, int Iy, int Iz, unsigned char data[Ix+2][Iy+2][Iz+2],
				 unsigned char newData[Ix+2][Iy+2][Iz+2], int X, int Y, int rank) {
	int x, y, z;
	x = X+1;
	y = Y+1;
	unsigned char c, old;
	int zLim = Iz+1;

	for (z=1;z<zLim;z++) {
		c = count(Ix, Iy, Iz, data, x, y, z);
		old = data[x][y][z];
		newData[X][Y][z] = transition(c, old);
	}
}

/* Calculate the new values of an edge normal to the YZ plane */
void updateYZedge(int Ix, int Iy, int Iz, unsigned char data[Ix+2][Iy+2][Iz+2],
				 unsigned char newData[Ix+2][Iy+2][Iz+2], int Y, int Z, int rank) {
	int x, y, z;
	z = Z+1;
	y = Y+1;
	unsigned char c, old;
	int xLim = Ix+1;

	for (x=1;x<xLim;x++) {
		c = count(Ix, Iy, Iz, data, x, y, z);
		old = data[x][y][z];
		newData[x][Y][Z] = transition(c, old);
	}
}

/* Calculate the new values of an edge normal to the ZX plane */
void updateZXedge(int Ix, int Iy, int Iz, unsigned char data[Ix+2][Iy+2][Iz+2],
				 unsigned char newData[Ix+2][Iy+2][Iz+2], int Z, int X, int rank) {
	int x, y, z;
	x = X+1;
	z = Z+1;
	unsigned char c, old;
	int yLim = Iy+1;

	for (y=1;y<yLim;y++) {
		c = count(Ix, Iy, Iz, data, x, y, z);
		old = data[x][y][z];
		newData[X][y][Z] = transition(c, old);
	}
}


/* A helper function to make update corners less cluttered */
void updateIndex(int Ix, int Iy, int Iz, unsigned char data[Ix+2][Iy+2][Iz+2],
				 unsigned char newData[Ix+2][Iy+2][Iz+2], int x, int y, int z) {

	unsigned char c, old;
	c = count(Ix, Iy, Iz, data, x, y, z);
	old = data[x][y][z];
	newData[x][y][z] = transition(c, old);
}

/* Calculate the new values of a corner */
void updateCorners(int Ix, int Iy, int Iz, unsigned char data[Ix+2][Iy+2][Iz+2],
				 unsigned char newData[Ix+2][Iy+2][Iz+2], int rank) {

	
	int x=1; int y=1;int z=1;
	updateIndex(Ix, Iy, Iz, data, newData, x, y, z);

	x=1;y=1;z=Iz;
	updateIndex(Ix, Iy, Iz, data, newData, x, y, z);

	x=1;y=Iy;z=1;
	updateIndex(Ix, Iy, Iz, data, newData, x, y, z);

	x=1;y=Iy;z=Iz;
	updateIndex(Ix, Iy, Iz, data, newData, x, y, z);

	x=Ix;y=1;z=1;
	updateIndex(Ix, Iy, Iz, data, newData, x, y, z);

	x=Ix;y=1;z=Iz;
	updateIndex(Ix, Iy, Iz, data, newData, x, y, z);

	x=Ix;y=Iy;z=1;
	updateIndex(Ix, Iy, Iz, data, newData, x, y, z);

	x=Ix;y=Iy;z=Iz;
	updateIndex(Ix, Iy, Iz, data, newData, x, y, z);
}

