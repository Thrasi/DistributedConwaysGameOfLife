#ifndef UPDATE_FUNCTIONS_H
#define UPDATE_FUNCTIONS_H

/* calculate the new values of the center cells*/
inline void updateCenter(int Ix, int Iy, int Iz, unsigned char *data,
				   unsigned char *newData, int rank);

/* calculate the new values of the sides in the X direction*/
inline void updateXside(int Ix, int Iy, int Iz, unsigned char *data,
				 unsigned char *newData, int x, int rank);

/* calculate the new values of the sides in the X direction*/
inline void updateYside(int Ix, int Iy, int Iz, unsigned char *data,
				 unsigned char *newData, int y, int rank);

/* calculate the new values of the sides in the X direction*/
inline void updateZside(int Ix, int Iy, int Iz, unsigned char *data,
				 unsigned char *newData, int z, int rank);

/* Calculate the new values of an edge normal to the XY plane */
inline void updateXYedge(int Ix, int Iy, int Iz, unsigned char *data,
				 unsigned char *newData, int x, int y, int rank);

/* Calculate the new values of an edge normal to the YZ plane */
inline void updateYZedge(int Ix, int Iy, int Iz, unsigned char *data,
				 unsigned char *newData, int y, int z, int rank);

/* Calculate the new values of an edge normal to the ZX plane */
inline void updateZXedge(int Ix, int Iy, int Iz, unsigned char *data,
				 unsigned char *newData, int z, int x, int rank);

#endif