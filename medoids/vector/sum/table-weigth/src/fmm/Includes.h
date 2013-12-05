/*
 * Include.h
 *
 *  Created on: 15/12/2009
 *      Author: user
 */

#ifndef INCLUDES_H_
#define INCLUDES_H_

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cfloat>
#include <climits>
#include <cctype>
#include <cmath>
#include <ctime>
#include <cassert>

#ifdef WINDOWS
#include <process.h>
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <vector>
#include <deque>
#include <queue>
#include <stack>
#include <map>
#include <set>
#include <bitset>
#include <stdexcept>
#include <exception>

#ifdef WINDOWS
#include <exception_defines.h>
#else
#include <bits/exception_defines.h>
#endif

using namespace std;

typedef long long ll;
typedef unsigned long long ull;

typedef istringstream is;
typedef ostringstream os;

#define eps 1e-10
#define inf 0x3f3f3f3f
#define pi acos(-1.0)

#define _inline(f...) f() __attribute__((always_inline)); f

#define mp(a,b) make_pair(a,b)

#define fr(x,y,z) for(int(x) =(y);(x)<(z);(x)++)
#define fi(n) fr(i,0,n)
#define fj(n) fr(j,0,n)
#define fk(n) fr(k,0,n)

#define cast(x,t) *({stringstream ss;static t __ret;ss<<x,ss>>__ret;&__ret;})

#define tr(container,it) for(typeof(container.begin()) it = container.begin(); it != container.end(); it++)

#define _dbg(x,y) x << #y << " == " << y << endl
#define _print(x,y) x << y << endl
#define dbg(x) // _dbg(cerr,x)
#define print(x) // _print(cerr,x)

#define dbg7(x)  _dbg(cerr,x)
#define print7(x)  _print(cerr,x)

int cmp(double x, double y = 0);
double sqr(double x);

enum Funcao {
	ABS, POW
};

enum Relacao {
	INDIVIDUAL, GRUPO
};

enum Coeficiente {
	INDIVIDUO, TABELA, TABELA_CLUSTER
};

#endif /* INCLUDES_H_ */
