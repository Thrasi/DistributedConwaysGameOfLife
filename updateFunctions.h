#ifndef UPDATE_FUNCTIONS_H
#define UPDATE_FUNCTIONS_H

/* calculate the new values of the center cells*/
void updateCenter(int Ix, int Iy, int Iz, unsigned char data[Ix][Iy][Iz],
				   unsigned char newData[Ix][Iy][Iz], int rank);

/* calculate the new values of the sides in the X direction*/
void updateXside(int Ix, int Iy, int Iz, unsigned char data[Ix+2][Iy+2][Iz+2],
				 unsigned char newData[Ix][Iy][Iz], int x, int rank);

/* calculate the new values of the sides in the X direction*/
void updateYside(int Ix, int Iy, int Iz, unsigned char data[Ix+2][Iy+2][Iz+2],
				 unsigned char newData[Ix][Iy][Iz], int y, int rank);

/* calculate the new values of the sides in the X direction*/
void updateZside(int Ix, int Iy, int Iz, unsigned char data[Ix+2][Iy+2][Iz+2],
				 unsigned char newData[Ix][Iy][Iz], int z, int rank);

/* Calculate the new values of an edge normal to the XY plane */
void updateXYedge(int Ix, int Iy, int Iz, unsigned char data[Ix+2][Iy+2][Iz+2],
				 unsigned char newData[Ix][Iy][Iz], int x, int y, int rank);

/* Calculate the new values of an edge normal to the YZ plane */
void updateYZedge(int Ix, int Iy, int Iz, unsigned char data[Ix+2][Iy+2][Iz+2],
				 unsigned char newData[Ix][Iy][Iz], int y, int z, int rank);

/* Calculate the new values of an edge normal to the ZX plane */
void updateZXedge(int Ix, int Iy, int Iz, unsigned char data[Ix+2][Iy+2][Iz+2],
				 unsigned char newData[Ix][Iy][Iz], int z, int x, int rank);

#endif