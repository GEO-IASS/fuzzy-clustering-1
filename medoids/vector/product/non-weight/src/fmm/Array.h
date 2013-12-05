/*
 * Array.h
 *
 *  Created on: 15/12/2009
 *      Author: user
 */

#ifndef ARRAY_H_
#define ARRAY_H_

#include "Includes.h"

template<class __Tp>
struct Array : public vector<__Tp> {
	__Tp& operator[](size_t __n) {
		assert(__n < this->size());
		return vector<__Tp>::operator[](__n);
	}
	const __Tp& operator[](size_t __n) const {
		assert(__n < this->size());
		return vector<__Tp>::operator[](__n);
	}
};

#endif /* ARRAY_H_ */
