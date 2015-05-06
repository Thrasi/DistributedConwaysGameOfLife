#ifndef UPDATE_FUNCTIONS_H
#define UPDATE_FUNCTIONS_H

/* Count the number of active cells around position (x,y,z) */
unsigned char count(int Ix, int Iy, int Iz, unsigned char data[Ix][Iy][Iz], int x, int y, int z);

/* Given an old cell value and a count of neighbours calculate the new one */
unsigned char transition(unsigned char c, unsigned char old);

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

#endif