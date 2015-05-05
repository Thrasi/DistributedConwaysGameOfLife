#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

struct triple {
	int p, q, r;
};

struct triple optimize(int w, int h, int d, int R);

#endif