#include "optimization.h"
#include <float.h>

/* Takes the size of the data w, h and d and number of processors R
and returns the parameters (dimensions) of the optimal topology. */
struct triple optimize(int w, int h, int d, int R) {
	int p, q, r, bp, bq, br;
	double phi = DBL_MAX;
	struct triple optParams;

	for (p = 1; p <= R; p++) {
		if (R % p != 0) {
			continue;
		}
		int rest = R / p;
		
		for (q = 1; q <= rest; q++) {
			if (rest % q != 0) {
				continue;
			}
			r = rest / q;

			double dp = (double) p, dq = (double) q, dr = (double) r;
			double nphi = dp / w + dq / h + dr / d + 2 * (dp*dq/(w*h) + dp*dr/(w*d) + dq*dr/(h*d));
			if (nphi < phi) {
				phi = nphi;
				bp = p;
				bq = q;
				br = r;
			}
		}
	}

	optParams.p = bp;
	optParams.q = bq;
	optParams.r = br;
	return optParams;
}