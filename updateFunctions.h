#ifndef UPDATE_FUNCTIONS_H
#define UPDATE_FUNCTIONS_H

/* calculate the new values of the center cells*/
void updateCenter(int Ix, int Iy, int Iz, unsigned char data[Ix][Iy][Iz],
				   unsigned char newData[Ix][Iy][Iz], int rank);

/* calculate the new values of the sides in the X direction*/
void updateXside(int Ix, int Iy, int Iz, unsigned char data[Ix+2][Iy+2][Iz+2],
				 unsigned char newData[Ix][Iy][Iz], int X, int rank);

/* calculate the new values of the sides in the X direction*/
void updateYside(int Ix, int Iy, int Iz, unsigned char data[Ix+2][Iy+2][Iz+2],
				 unsigned char newData[Ix][Iy][Iz], int Y, int rank);

/* calculate the new values of the sides in the X direction*/
void updateZside(int Ix, int Iy, int Iz, unsigned char data[Ix+2][Iy+2][Iz+2],
				 unsigned char newData[Ix][Iy][Iz], int Z, int rank);

/* Calculate the new values of an edge normal to the XY plane */
void updateXYedge(int Ix, int Iy, int Iz, unsigned char data[Ix+2][Iy+2][Iz+2],
				 unsigned char newData[Ix][Iy][Iz], int X, int Y, int rank);

/* Calculate the new values of an edge normal to the YZ plane */
void updateYZedge(int Ix, int Iy, int Iz, unsigned char data[Ix+2][Iy+2][Iz+2],
				 unsigned char newData[Ix][Iy][Iz], int Y, int Z, int rank);

/* Calculate the new values of an edge normal to the ZX plane */
void updateZXedge(int Ix, int Iy, int Iz, unsigned char data[Ix+2][Iy+2][Iz+2],
				 unsigned char newData[Ix][Iy][Iz], int Z, int X, int rank);

#endif