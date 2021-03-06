
#include "updateFunctions.h"

#define xoff (Iy+2)*(Iz+2)
#define yoff (Iz+2)


/* Count the number of active cells around position (x,y,z) */
inline unsigned char count(int Ix, int Iy, int Iz, unsigned char *data, int x, int y, int z) {
	unsigned char c = 0;
    int dx,dy,dz;
    for (dx=-1;dx<=1;dx++) {
    	for (dy=-1;dy<=1;dy++) {
    		for (dz=-1;dz<=1;dz++) {
    			c += data[xoff*(x+dx) + yoff*(y+dy) + z+dz];
    		}
    	}	
    }
    c -= data[xoff*x + yoff*y + z];
	return c;
}

/* Given an old cell value and a count of neighbours calculate the new one */
// 55 54
inline unsigned char transition(unsigned char c, unsigned char old) {
	// if ( old == 1) {
	// 	if ( c == 4 || c == 5) {
	// 		return 1;
	// 	} else {
	// 		return 0;
	// 	}
	// } else if ( c == 5 ) {
	//         return 1;
	// } else {
	//         return 0;
	// }
	if ( old == 1) {
		if ( c == 4 || c == 5) {
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

/* calculate the new values of the center cells*/
inline void updateCenter(int Ix, int Iy, int Iz, unsigned char *data,
				   unsigned char *newData, int rank) {

	int x, y, z;
	unsigned char c, old;
	int xLim = Ix;
	int yLim = Iy;
	int zLim = Iz;

	for (x=2;x<xLim;x++) {
		for (y=2;y<yLim;y++) {
			for (z=2;z<zLim;z++) {
				c = count(Ix, Iy, Iz, data, x, y, z);
				old = data[xoff*x + yoff*y + z];
				newData[xoff*x + yoff*y + z] = transition(c, old);
			}
		}
	}
}

/* calculate the new values of a side with index X */
inline void updateXside(int Ix, int Iy, int Iz, unsigned char *data,
				 unsigned char *newData, int x, int rank) {
	int y, z;
	unsigned char c, old;

	for (y=2;y<Iy;y++) {
		for (z=2;z<Iz;z++) {
			c = count(Ix, Iy, Iz, data, x, y, z);
			old = data[xoff*x + yoff*y + z];
			newData[xoff*x + yoff*y + z] = transition(c, old);
		}
	}
}

/* calculate the new values of a side with index Y */
inline void updateYside(int Ix, int Iy, int Iz, unsigned char *data,
				 unsigned char *newData, int y, int rank) {
	int x, z;
	unsigned char c, old;

	for (x=1;x<Ix;x++) {
		for (z=1;z<Iz;z++) {
			c = count(Ix, Iy, Iz, data, x, y, z);
			old = data[xoff*x + yoff*y + z];
			newData[xoff*x + yoff*y + z] = transition(c, old);
		}
	}
}

/* calculate the new values of a side with index Z */
inline void updateZside(int Ix, int Iy, int Iz, unsigned char *data,
				 unsigned char *newData, int z, int rank) {
	int x, y;
	unsigned char c, old;

	for (y=2;y<Iy;y++) {
		for (x=2;x<Ix;x++) {
			c = count(Ix, Iy, Iz, data, x, y, z);
			old = data[xoff*x + yoff*y + z];
			newData[xoff*x + yoff*y + z] = transition(c, old);
		}
	}
}

/* Calculate the new values of an edge normal to the XY plane */
inline void updateXYedge(int Ix, int Iy, int Iz, unsigned char *data,
				 unsigned char *newData, int x, int y, int rank) {
	int z;
	unsigned char c, old;

	for (z=2;z<Iz;z++) {
		c = count(Ix, Iy, Iz, data, x, y, z);
		old = data[xoff*x + yoff*y + z];
		newData[xoff*x + yoff*y + z] = transition(c, old);
	}
}

/* Calculate the new values of an edge normal to the YZ plane */
inline void updateYZedge(int Ix, int Iy, int Iz, unsigned char *data,
				 unsigned char *newData, int y, int z, int rank) {
	int x;
	unsigned char c, old;

	for (x=2;x<Ix;x++) {
		c = count(Ix, Iy, Iz, data, x, y, z);
		old = data[xoff*x + yoff*y + z];
		newData[xoff*x + yoff*y + z] = transition(c, old);
	}
}

/* Calculate the new values of an edge normal to the ZX plane */
inline void updateZXedge(int Ix, int Iy, int Iz, unsigned char *data,
				 unsigned char *newData, int z, int x, int rank) {
	int y;
	unsigned char c, old;

	for (y=2;y<Iy;y++) {
		c = count(Ix, Iy, Iz, data, x, y, z);
		old = data[xoff*x + yoff*y + z];
		newData[xoff*x + yoff*y + z] = transition(c, old);
	}
}


/* A helper function to make update corners less cluttered */
inline void updateIndex(int Ix, int Iy, int Iz, unsigned char *data,
				 unsigned char *newData, int x, int y, int z) {

	unsigned char c, old;
	c = count(Ix, Iy, Iz, data, x, y, z);
	old = data[xoff*x + yoff*y + z];
	newData[xoff*x + yoff*y + z] = transition(c, old);
}

/* Calculate the new values of a corner */
inline void updateCorners(int Ix, int Iy, int Iz, unsigned char *data,
				 unsigned char *newData, int rank) {

	
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

