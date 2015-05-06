
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
				   unsigned char newData[Ix][Iy][Iz], int rank) {

	int x, y, z;
	unsigned char c, old;
	int xLim = Ix-1;
	int yLim = Iy-1;
	int zLim = Iz-1;

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
				 unsigned char newData[Ix][Iy][Iz], int X, int rank) {
	int x, y, z;
	x = X+1;
	unsigned char c, old;
	int yLim = Iy+1;
	int zLim = Iz+1;

	for (y=1;y<yLim;y++) {
		for (z=1;z<zLim;z++) {
			c = count(Ix, Iy, Iz, data, x, y, z);
			old = data[x][y][z];
			newData[X][y-1][z-1] = transition(c, old);
		}
	}
}

/* calculate the new values of a side with index Y */
void updateYside(int Ix, int Iy, int Iz, unsigned char data[Ix+2][Iy+2][Iz+2],
				 unsigned char newData[Ix][Iy][Iz], int Y, int rank) {
	int x, y, z;
	y = Y+1;
	unsigned char c, old;
	int xLim = Ix+1;
	int zLim = Iz+1;

	for (x=1;x<xLim;x++) {
		for (z=1;z<zLim;z++) {
			c = count(Ix, Iy, Iz, data, x, y, z);
			old = data[x][y][z];
			newData[x-1][Y][z-1] = transition(c, old);
		}
	}
}

/* calculate the new values of a side with index Z */
void updateZside(int Ix, int Iy, int Iz, unsigned char data[Ix+2][Iy+2][Iz+2],
				 unsigned char newData[Ix][Iy][Iz], int Z, int rank) {
	int x, y, z;
	z = Z+1;
	unsigned char c, old;
	int yLim = Iy+1;
	int xLim = Ix+1;

	for (y=1;y<yLim;y++) {
		for (x=1;x<xLim;x++) {
			c = count(Ix, Iy, Iz, data, x, y, z);
			old = data[x][y][z];
			// if (rank == 1) {
			// 	if (c > 0) {
			// 		printf("",)
			// 	}
			// }
			newData[x-1][y-1][Z] = transition(c, old);
		}
	}
}

// void updateSide(int Ix, int Iy, int Iz, int d1, int d2, unsigned char side[d1][d2], unsigned char data[Ix][Iy][Iz],
// 					unsigned char newData[Ix][Iy][Iz], int direction, int dimension, int rank) {
// 	if (dimension == 0) {

// 	} else if (dimension == 1) {

// 	} else {

// 	}
// }

/* calculate the new values of the sides in the X direction*/
// void updateXside(int Ix, int Iy, int Iz, unsigned char Xside[Iy][Iz], unsigned char data[Ix][Iy][Iz],
// 					unsigned char newData[Ix][Iy][Iz], int direction, int rank) {
// 	int x,y,z;
// 	if (direction > 1) {
// 		x = Ix-2;
// 	} else {
// 		x = 0;
// 	}

// 	int yLim = Iy-1;
// 	int zLim = Iz-1;
// 	unsigned char c, old;
//     int dx,dy,dz;
// 	for (y=1;y<yLim;y++) {
// 		for (z=1;z<zLim;z++) {
// 			c = 0;
// 			// 18 from data
// 		    for (dx=0;dx<=1;dx++) {
// 		    	for (dy=-1;dy<=1;dy++) {
// 		    		for (dz=-1;dz<=1;dz++) {
// 		    			c += data[x+dx][y+dy][z+dz];
// 		    		}
// 		    	}	
// 		    }
// 		    // 9 from 
// 		    for (dy=-1;dy<=1;dy++) {
// 	    		for (dz=-1;dz<=1;dz++) {
// 	    			c += Xside[y+dy][z+dz];
// 	    		}
// 	    	}

// 		    c -= data[x][y][z];

// 			old = data[x][y][z];
// 			// if (rank==5){
// 			// 	printf("%d,%d,%d: %hhu\n",x,y,z,c);
// 			// }
// 			newData[x][y][z] = transition(c, old);
// 		}
// 	}
// }

// /* calculate the new values of the sides in the Y direction*/
// void updateYside(int Ix, int Iy, int Iz, unsigned char Yside[Iz][Ix], unsigned char data[Ix][Iy][Iz],
// 					unsigned char newData[Ix][Iy][Iz], int direction, int rank) {
// 	int x,y,z;
// 	if (direction > 1) {
// 		y = Iy-2;
// 	} else {
// 		y = 0;
// 	}

// 	int xLim = Ix-1;
// 	int zLim = Iz-1;
// 	unsigned char c, old;
//     int dx,dy,dz;
// 	for (x=1;x<xLim;x++) {
// 		for (z=1;z<zLim;z++) {
// 			c = 0;
// 			// 18 from data
// 		    for (dx=-1;dx<=1;dx++) {
// 		    	for (dy=0;dy<=1;dy++) {
// 		    		for (dz=-1;dz<=1;dz++) {
// 		    			c += data[x+dx][y+dy][z+dz];
// 		    		}
// 		    	}	
// 		    }
// 		    // 9 from 
// 		    for (dx=-1;dx<=1;dx++) {
// 	    		for (dz=-1;dz<=1;dz++) {
// 	    			c += Yside[z+dz][x+dx];
// 	    		}
// 	    	}

// 		    c -= data[x][y][z];

// 			old = data[x][y][z];
// 			newData[x][y][z] = transition(c, old);
// 		}
// 	}
// }

// /* calculate the new values of the sides in the Z direction*/
// void updateZside(int Ix, int Iy, int Iz, unsigned char Zside[Ix][Iy], unsigned char data[Ix][Iy][Iz],
// 					unsigned char newData[Ix][Iy][Iz], int direction, int rank) {
// 	int x,y,z;
// 	if (direction > 1) {
// 		z = Iz-2;
// 	} else {
// 		z = 0;
// 	}

// 	int xLim = Ix-1;
// 	int yLim = Iy-1;
// 	unsigned char c, old;
//     int dx,dy,dz;
// 	for (x=1;x<xLim;x++) {
// 		for (y=1;y<yLim;y++) {
// 			c = 0;
// 			// 18 from data
// 		    for (dx=-1;dx<=1;dx++) {
// 		    	for (dy=-1;dy<=1;dy++) {
// 		    		for (dz=0;dz<=1;dz++) {
// 		    			c += data[x+dx][y+dy][z+dz];
// 		    		}
// 		    	}	
// 		    }
// 		    // 9 from 
// 		    for (dx=-1;dx<=1;dx++) {
// 	    		for (dy=-1;dy<=1;dy++) {
// 	    			c += Zside[x+dx][y+dy];
// 	    		}
// 	    	}

// 		    c -= data[x][y][z];

// 			old = data[x][y][z];
// 			newData[x][y][z] = transition(c, old);
// 		}
// 	}
// }