/*
 * Includes.cpp
 *
 *  Created on: 17/12/2009
 *      Author: user
 */

#include "Includes.h"

_inline(int cmp)
(double x, double y) {
	return (x <= y + eps) ? (x + eps < y) ? -1 : 0 : 1;
}

_inline(double sqr)
(double x) {
	return (x * x);
}
