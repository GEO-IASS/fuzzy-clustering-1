//#pragma comment(linker, "/STACK:1073741824")

/*
 * The_CARD_R_Algorithm.cpp
 *
 *  Created on: 08/03/2009
 *      Author: Filipe
 */
#include <cstdio>
#include <cstdlib>
#include <cmath>

//#define NDEBUG
#include <cassert>
#include <ctime>
#include <cfloat>
#include <climits>
#include <cstring>

#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>
#include <sstream>
#include <set>
#include <map>
#include <queue>
#include <stack>

#include "validation.h"

using namespace std;

#define _inline(f...) f() __attribute__((always_inline)); f
#define EPS 1e-7
_inline(int cmp)
(double x, double y = 0.0) {
	return (x <= y + EPS) ? (x + EPS < y) ? -1 : 0 : 1;
}
#undef EPS
#undef _inline

///////////////////////////////////////////////////////////
#define INF 0x3f3f3f3f
#define MAXV 11111
#define MAXE 11111
#define LIM (1<<20)

struct MinCostFlow {

	typedef int T;

	// Heap
	struct Heap {
		T no, dist;
		bool operator<(const Heap& h) const {
			return dist > h.dist;
		}
	} heap[LIM], atual;
	int Q;

	// Edge
	struct Edge {
		T u, v, cap, cost, ant;
	} edge[MAXE];
	int E;

	// Graph
	T adj[MAXV], pai[MAXV], dist[MAXV], pot[MAXV];
	int V;

	void init(int n) {
		memset(adj,-1,sizeof(adj));
		V = n, E = 0;
	}

	void add_edge(int u, int v, T cap, T cost, bool rev = false) {
		edge[E] = (Edge){ u, v, cap, cost, adj[u]};
		adj[u] = E++;

		edge[E] = (Edge){ v, u, 0, -cost, adj[v]};
		adj[v] = E++;

		if(rev) {
			add_edge(v,u,cap,cost);
		}
	}

	Heap& top() {
		pop_heap(heap,heap+Q);
		return heap[--Q];
	}

	void update(Heap& h, int e = -1) {
		if(h.dist < dist[h.no]) {
			dist[h.no] = h.dist, pai[h.no] = e;
			heap[Q++] = h;
			push_heap(heap,heap+Q);
		}
	}

	bool bellman_ford(int s) {
		memset(dist,INF,sizeof(dist));
		dist[s] = 0;
		for(int k = 0; k < V; k++) {
			for(int i = 0; i < V; i++) {
				for(int j = adj[i]; j != -1; j = edge[j].ant) {
					if(edge[j].cap <= 0) continue;
					T u = edge[j].u, v = edge[j].v, cost = edge[j].cost;
					if(dist[u] + cost < dist[v]) {
						dist[v] = dist[u] + cost;
					}
				}
			}
		}
		for(int i = 0; i < V; i++) {
			for(int j = adj[i]; j != -1; j = edge[j].ant) {
				if(edge[j].cap <= 0) continue;
				T u = edge[j].u, v = edge[j].v, cost = edge[j].cost;
				if(dist[u] + cost < dist[v]) { // negative cycle
					return true;
				}
			}
		}
		return false;
	}

	bool dijkstra(int s, int t) {
		memset(dist,INF,sizeof(dist)), Q = 0;
		memset(pai,-1,sizeof(pai));
		Heap h = {s,0};
		update(h);

		while(Q) {
			atual = top();
			if(atual.dist > dist[atual.no]) continue;
			for(int i = adj[atual.no]; i != -1; i = edge[i].ant) {
				if(edge[i].cap <= 0) continue;
				T u = edge[i].u, v = edge[i].v, cost = edge[i].cost;
				Heap h = { v, atual.dist + cost + pot[u] - pot[v] };
				update(h, i);
			}
		}

		return dist[t] < INF;
	}

	pair<T,T> solve(int s, int t) {

		memset(pot,0,sizeof(pot));

		/* dealing with negative costs */
		if(bellman_ford(s)) {
			assert(0); // negative cycle
		}

		for(int i = 0; i < V; i++) {
			pot[i] += dist[i];
		}
		/* dealing with negative costs */

		T flow = 0, cost = 0;
		while( dijkstra(s,t) ) {
			T vmin = INF;
			for(int i = pai[t]; i != -1; i = pai[edge[i].u]) {
				vmin = min(vmin, edge[i].cap);
			}
			flow += vmin;
			for(int i = pai[t]; i != -1; i = pai[edge[i].u]) {
				cost += (edge[i].cost * vmin);
				edge[i].cap -= vmin;
				edge[i^1].cap += vmin;
			}
			for(int i = 0; i < V; i++) {
				pot[i] += dist[i];
			}
		}

		return make_pair(flow,cost);
	}
} mcf;

#undef INF
#undef MAXV
#undef MAXE
#undef LIM
///////////////////////////////////////////////////////////


int var_classe;

struct Matrix {
	vector<vector<double> > v;

	Matrix(int row = 1, int column = 1) {
		resize(row, column);
	}

	Matrix(const char* filename) {
		readDissimilarityMatrix(filename);
	}

	virtual ~Matrix() {
	}

	int getNumRow() const {
		return v.size();
	}

	int getNumColumn() const {
		return v[0].size();
	}

	double getMaxElement() const {
		int r = getNumRow();
		int c = getNumColumn();

		double maior = DBL_MIN;
		for (int i = 0; i < r; i++) {
			for (int j = 0; j < c; j++) {
				maior = max(maior, (*this)(i, j));
			}
		}

		return maior;
	}

	void resize(int row, int column) {
		assert(row> 0 && column> 0);

		v.resize(row);
		for (int i = 0; i < row; i++) {
			v[i].resize(column);
		}
	}

	static double getEuclidianNorm(const Matrix& ma, const Matrix& mb) {
		int ra = ma.getNumRow(), ca = ma.getNumColumn();
		int rb = mb.getNumRow(), cb = mb.getNumColumn();

		assert(ra == rb && ca == cb);

		double norm = 0.0;
		for (int i = 0; i < ra; i++) {
			for (int j = 0; j < cb; j++) {
				norm += std::pow(ma(i, j) - ma(i, j), 2.0);
			}
		}

		return norm;
	}

	void readDissimilarityMatrix(const char* filename) {
		ifstream soda(filename, ios::in);
		string linha;
		int n;

		assert(soda != NULL);

		// read "indiv_nb"
		do {
			soda >> linha;
		} while (linha.compare("indiv_nb"));

		// read "="
		soda >> linha;
		assert(linha.compare("=") == 0);

		// read the individual number
		soda >> n;
		assert(n> 0);

		resize(n, n);

		do {
			getline(soda, linha);
		} while (linha.find("DIST_MATRIX= (") == string::npos);

		char c;
		double d;

		//size_t pos;
		for (int i = 0; i < n; i++) {
			linha.clear();
			while (soda >> c && c != '(')
				;
			linha.push_back(c);
			while (soda >> c && c != ')') {
				linha.push_back(c);
			}
			linha.push_back(c);
			while (soda >> c && c != ',') {
				linha.push_back(c);
			}

			linha[linha.find("(")] = ' ';
			linha[linha.find(")")] = ' ';

			/*while ((pos = linha.find(",")) != string::npos)
				linha[pos] = ' ';*/

			int tlinha = linha.size();
			for(int tt = 0; tt < tlinha; tt++) if(linha[tt] == ',') {
				linha[tt] = ' ';
			}

			istringstream str(linha);
			for (int j = 0; j <= i; j++) {
				if (!(str >> d)) {
					__throw_ios_failure("Arquivo invalido!!");
				}
				(*this)(i, j) = (*this)(j, i) = d;
			}
		}

		soda.close();
	}

	Matrix getRow(int row) const {
		int r = getNumRow();
		int c = getNumColumn();

		assert(0 <= row && row < r);

		Matrix m(1, c);

		for (int i = 0; i < c; i++) {
			m(0, i) = (*this)(row, i);
		}

		return m;
	}

	Matrix getColumn(int column) const {
		int r = getNumRow();
		int c = getNumColumn();

		assert(0 <= column && column < c);

		Matrix m(r, 1);

		for (int i = 0; i < r; i++) {
			m(i, 0) = (*this)(i, column);
		}

		return m;
	}

	Matrix& abs() {
		int r = getNumRow();
		int c = getNumColumn();

		for (int i = 0; i < r; i++) {
			for (int j = 0; j < c; j++) {
				(*this)(i, j) = std::fabs((*this)(i, j));
			}
		}

		return (*this);
	}

	static Matrix getIndentityMatrix(int n) {
		Matrix identity(n, n);

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i == j) {
					identity(i, j) = 1.0;
				} else {
					identity(i, j) = 0.0;
				}
			}
		}
		return identity;
	}

	static Matrix getAntiIndentityMatrix(int n) {
		Matrix antiIdentity(n, n);

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i == j) {
					antiIdentity(i, j) = 0.0;
				} else {
					antiIdentity(i, j) = 1.0;
				}
			}
		}

		return antiIdentity;
	}

	static Matrix getRandInitialU(int row, int column) {
		Matrix m(row, column);

		for (int i = 0; i < row; i++) {
			for (int j = 0; j < column; j++) {
				m(i, j) = std::rand();
			}
		}

		double sum;
		for (int i = 0; i < column; i++) {
			sum = 0.0;
			for (int j = 0; j < row; j++) {
				sum += m(j, i);
			}

			assert(cmp(sum)> 0);

			for (int j = 0; j < row; j++) {
				m(j, i) /= sum;
			}
		}

		return m;
	}

	double operator()(int row, int column) const {
		assert( 0 <= row && row < getNumRow() );
		assert( 0 <= column && column < getNumColumn() );

		return v[row][column];
	}

	double& operator()(int row, int column) {
		assert( 0 <= row && row < getNumRow() );
		assert( 0 <= column && column < getNumColumn() );

		return v[row][column];
	}

	Matrix operator+(const Matrix& m) const {
		int r = getNumRow();
		int c = getNumColumn();

		assert( r == m.getNumRow() );
		assert( c == m.getNumColumn() );

		Matrix sum(r, c);

		for (int i = 0; i < r; i++) {
			for (int j = 0; j < c; j++) {
				sum(i, j) = (*this)(i, j) + m(i, j);
			}
		}

		return sum;
	}

	Matrix operator-(const Matrix& m) const {
		int r = getNumRow();
		int c = getNumColumn();

		assert( r == m.getNumRow() );
		assert( c == m.getNumColumn() );

		Matrix subtraction(r, c);

		for (int i = 0; i < r; i++) {
			for (int j = 0; j < c; j++) {
				subtraction(i, j) = (*this)(i, j) - m(i, j);
			}
		}

		return subtraction;
	}

	Matrix operator*(double val) const {
		int r = getNumRow();
		int c = getNumColumn();

		Matrix product = *this;

		for (int i = 0; i < r; i++) {
			for (int j = 0; j < c; j++) {
				product(i, j) *= val;
			}
		}

		return product;
	}

	Matrix operator/(double val) const {
		assert(cmp(val) != 0);

		int r = getNumRow();
		int c = getNumColumn();

		Matrix division = *this;

		for (int i = 0; i < r; i++) {
			for (int j = 0; j < c; j++) {
				division(i, j) /= val;
			}
		}

		return division;
	}

	Matrix operator~() const {
		int r = getNumRow();
		int c = getNumColumn();

		Matrix trans(c, r);

		for (int i = 0; i < r; i++) {
			for (int j = 0; j < c; j++) {
				trans(j, i) = (*this)(i, j);
			}
		}

		return trans;
	}

	Matrix operator^(char val) const {
		assert(val == 'T' || val == 't');
		return ~(*this);
	}

	Matrix operator^(int val) const {
		int r = getNumRow();
		int c = getNumColumn();

		Matrix pot = *this;

		for (int i = 0; i < r; i++) {
			for (int j = 0; j < c; j++) {
				pot(i, j) = std::pow((*this)(i, j), val);
			}
		}

		return pot;
	}

	Matrix operator*(const Matrix& m) const {
		int r = getNumRow();
		int c = getNumColumn();

		int r_m = m.getNumRow();
		int m_c = m.getNumColumn();

		assert(c == r_m);

		Matrix product(r, m_c);

		for (int i = 0; i < r; i++) {
			for (int j = 0; j < m_c; j++) {
				product(i, j) = 0.0;

				for (int k = 0; k < c; k++) {
					product(i, j) += (*this)(i, k) * m(k, j);
				}
			}
		}

		return product;
	}

	friend ostream& operator<<(ostream& out, const Matrix& m) {
		const static string SPACE(" ");

		int r = m.getNumRow();
		int c = m.getNumColumn();

		out << fixed << setprecision(6);

		for (int i = 0; i < r; i++) {
			out << "|";
			for (int j = 0; j < c; j++) {
				out << SPACE << m(i, j);
			}
			out << SPACE << "|";

			if (i + 1 < r) {
				out << "\n";
			}
		}

		return out;
	}
};

struct Cluster {
	set<int, less<int> > conjunto;

	Cluster() {
		clear();
	}

	virtual ~Cluster() {

	}

	void add(int valor) {
		conjunto.insert(valor);
	}

	void clear() {
		conjunto.clear();
	}

	bool contem(int valor) const {
		return (conjunto.find(valor) != conjunto.end());
	}

	int size() const {
		return conjunto.size();
	}

	int operator&&(const Cluster& c) const {
		int intersection = 0;
		set<int, less<int> >::iterator iter;

		for (iter = conjunto.begin(); iter != conjunto.end(); iter++) {
			if (c.contem(*iter)) {
				intersection++;
			}
		}

		return intersection;
	}

	friend ostream& operator<<(ostream& out, const Cluster& c) {
		set<int, less<int> >::iterator iter = c.conjunto.begin();
		int tam = c.size();

		out << "Elements: {";
		for (int i = 0; i < tam; i++, iter++) {
			out << " " << (*iter);

			if (i + 1 < tam) {
				out << ",";
			}
		}
		out << " }";

		return out;
	}
};

vector<Cluster> readPrioriCluster(int n, const char* filename) {
	vector<Cluster> priori;

	if(var_classe == -1) {
		return priori;
	}

	ifstream soda(filename, ios::in);

	string linha;
	do {
		getline(soda, linha);
	} while (linha.find("RECTANGLE_MATRIX = (") == string::npos);

	char c = '?';
	int ct;
	size_t pri = 0;
	size_t pos;
	vector<string> v;
	string piece;

	for (int i = 0; i < n; i++) {

		linha.clear();

		ct = 0;
		while (!ct) {
			if (!(soda >> c))
				__throw_ios_failure("Arquivo invalido");
			if (c == '(')
				ct++;
		}

		while (ct) {
			if (!(soda >> c))
				__throw_ios_failure("Arquivo invalido");
			if (c == '(')
				ct++;
			if (c == ')')
				ct--;
			if (ct)
				linha.push_back(c);
		}

		ct = 0;
		piece.clear();
		v.clear();
		
		linha.push_back(',');

    for(unsigned int j = 0; j < linha.size(); ++j) if(linha[j] != ' ') {
      if (linha[j] == '(') ct++; else if (linha[j] == ')') ct--;
      if(!ct && linha[j] == ',') {
        v.push_back(piece), piece.clear();
      } else piece.push_back(linha[j]);
    }
    
		piece = v[var_classe - 1];
		assert(sscanf(piece.c_str(),"%u",&pri) == 1);
		assert(pri > 0);

		// for (size_t j = a.prioriCluster.size(); j < priori; j++) {
		// a.prioriCluster.push_back(Cluster(0));
		//
		for (size_t j = priori.size(); j < pri; j++) {
			priori.push_back(Cluster());
		}
		//a.prioriCluster[priori - 1].insert(i);
		priori[pri - 1].add(i);
	}

	soda.close();

	return priori;
}
//===============================================
// PRINTERS FUNCTIONS
//===============================================
#define DISTANCIA 5
#define DIVISA '|'
#define TRACO '-'
void printLinha1(ostream &out) {
	out << DIVISA;
	for (int i = 0; i < 38; i++) {
		out << TRACO;
	}
	for (int i = 6; i < 6 * DISTANCIA; i++) {
		out << TRACO;
	}
	out << DIVISA << '\n';
}

void printLinha2(ostream &out) {
	out << DIVISA;
	for (int i = 0; i < 11; i++) {
		out << TRACO;
	}
	for (int i = 2; i < 2 * DISTANCIA; i++) {
		out << TRACO;
	}
	out << DIVISA;
	for (int i = 0; i < 15; i++) {
		out << TRACO;
	}
	for (int i = 2; i < 2 * DISTANCIA; i++) {
		out << TRACO;
	}
	out << DIVISA;
	for (int i = 0; i < 10; i++) {
		out << TRACO;
	}
	for (int i = 2; i < 2 * DISTANCIA; i++) {
		out << TRACO;
	}
	out << DIVISA << '\n';
}

void printInicializacaoHeader(ostream &out, int inicializacao) {
	out << DIVISA;
	out << " # Initialization: ";
	out << setw(19 + (DISTANCIA - 1) * 6) << left << inicializacao;
	out << DIVISA << '\n';
}
void print2(ostream &out) {
	out << DIVISA;
	out << setw(DISTANCIA) << " ";

	out << "Iteration";

	out << setw(DISTANCIA) << " ";
	out << DIVISA;
	out << setw(DISTANCIA) << " ";

	out << setw(6) << " ";
	out << "J";
	out << setw(6) << " ";

	out << setw(DISTANCIA) << " ";
	out << DIVISA;
	out << setw(DISTANCIA) << " ";

	out << setw(3) << " ";
	out << "CR";
	out << setw(3) << " ";

	out << setw(DISTANCIA) << " ";
	out << DIVISA << '\n';
}
void print3(ostream &out, int iteracao, double J, double CR) {
	out << DIVISA;
	out << setw(DISTANCIA) << " ";

	out << setfill('0');
	out << setw(9) << right << iteracao;
	out << setfill(' ');

	out << setw(DISTANCIA) << " ";
	out << DIVISA;
	out << setw(DISTANCIA) << " ";

	out << scientific << J;

	out << setw(DISTANCIA) << " ";
	out << DIVISA;
	out << setw(DISTANCIA) << " ";

	out << fixed << setprecision(6) << CR;

	out << setw(DISTANCIA) << " ";
	out << DIVISA << '\n';
}
#undef DISTANCIA
#define DISTANCIA 2
//=======================================================
void printLinha1(ostream &out, int clus) {
	out << DIVISA;
	for (int i = 0; i < 43; i++) {
		out << TRACO;
	}
	for (int i = 0; i < 15 * clus; i++) {
		out << TRACO;
	}
	int espaco = DISTANCIA * (3 + clus);
	for (int i = 3 + clus; i < espaco; i++) {
		out << TRACO << TRACO;
	}
	out << DIVISA << '\n';
}

void printLinha2(ostream &out, int clus) {
	out << DIVISA;
	for (int i = 0; i < 9; i++) {
		out << TRACO;
	}
	for (int i = 2; i < 2 * DISTANCIA; i++) {
		out << TRACO;
	}
	out << DIVISA;

	for (int j = 0; j < clus; j++) {
		for (int i = 0; i < 14; i++) {
			out << TRACO;
		}
		for (int i = 2; i < 2 * DISTANCIA; i++) {
			out << TRACO;
		}
		out << DIVISA;
	}

	for (int i = 0; i < 14; i++) {
		out << TRACO;
	}
	for (int i = 2; i < 2 * DISTANCIA; i++) {
		out << TRACO;
	}
	out << DIVISA;

	for (int i = 0; i < 18; i++) {
		out << TRACO;
	}
	for (int i = 2; i < 2 * DISTANCIA; i++) {
		out << TRACO;
	}
	out << DIVISA << '\n';
}

void printRelacaoHeader(ostream &out, int clus) {
	out << DIVISA;
	out << setw(DISTANCIA) << " ";

	out << "Pattern";

	out << setw(DISTANCIA) << " ";
	out << DIVISA;

	for (int i = 0; i < clus; i++) {
		out << setw(DISTANCIA) << " ";
		out << "cluster[";
		out << setfill('0');
		out << setw(3) << (i + 1);
		out << setfill(' ');
		out << "]";
		out << setw(DISTANCIA) << " ";
		out << DIVISA;
	}

	out << setw(DISTANCIA) << " ";

	out << "Hard Cluster";

	out << setw(DISTANCIA) << " ";

	out << DIVISA;
	out << setw(DISTANCIA) << " ";

	out << "A Priori Cluster";

	out << setw(DISTANCIA) << " ";

	out << DIVISA << '\n';
}
void print2(ostream &out, int p, int clus, int select, int priori) {
	out << DIVISA;
	out << setw(DISTANCIA) << " ";

	out << setfill('0');
	out << setw(7) << right << p;
	out << setfill(' ');

	out << setw(DISTANCIA) << " ";

	for (int i = 0; i < clus; i++) {
		out << DIVISA;
		out << setw(DISTANCIA) << " ";

		out << setw(5) << " ";
		if (i == select) {
			out << setw(1) << '1';
		} else {
			out << setw(1) << '0';
		}
		out << setw(6) << " ";
		out << setw(DISTANCIA) << " ";
	}

	out << DIVISA;
	out << setw(DISTANCIA) << " ";

	out << setw(5) << " ";
	out << setfill('0');
	out << setw(3) << (select + 1);
	out << setfill(' ');
	out << setw(4) << " ";

	out << setw(DISTANCIA) << " ";
	out << DIVISA;

	out << setw(DISTANCIA) << " ";

	out << setw(7) << " ";
	out << setfill('0');
	out << setw(3) << (priori + 1);
	out << setfill(' ');
	out << setw(6) << " ";

	out << setw(DISTANCIA) << " ";
	out << DIVISA << '\n';
}
#undef DISTANCIA
#undef DIVISA
#undef TRACO
//===============================================
//===============================================

int readNumeroIndividuos(const char *filename) {
	ifstream soda(filename, ios::in);
	assert(soda != NULL);

	string word;
	int individuos;

	// read "indiv_nb"
	do {
		soda >> word;
	} while (word.compare("indiv_nb"));

	// read "="
	soda >> word;
	assert(word.compare("=") == 0);

	// read the individual number
	soda >> individuos;
	assert(individuos> 0);

	// close file
	soda.close();

	return individuos;
}

#define prox(a) (a + 1) % 2
#define comb(a) ( (a) * ( (a) - 1 ) / 2 )
struct Algorithm {
	vector<string> files;
	vector<Matrix> R;
	Matrix I, M;
	int m, q, s;
	int n;
	int c, p;

	int r; // maximal number of iterations

	Matrix U[2], W;
	int pos;
	double beta;

	const vector<Cluster> prioriCluster;
	vector<Cluster> cluster;

	int init;

	int initialization;
	int iteration;

	int melhorInicializacao;
	double melhorJ;
	double melhorCR;
	Matrix melhorU;
	Matrix melhorW;
	vector<Cluster> melhorCluster;

	ostream& out;

	Algorithm(vector<string> _files, int _c, int _m, int _q, int _init, int _r,
			ostream& saida) :
		n(readNumeroIndividuos(_files[0].c_str())), prioriCluster(
				readPrioriCluster(n, _files[0].c_str())), out(saida) {
		files = _files;
		m = _m;
		q = _q;
		s = files.size();

		c = _c;

		r = _r;

		R.resize(s);
		for (int i = 0; i < s; i++) {
			printf("lendo tabela [%d]\n",i+1); fflush(stdout);
			R[i] = Matrix(files[i].c_str());
		}
		printf("leu TUDO\n"); fflush(stdout);
		I = Matrix::getIndentityMatrix(n);
		M = Matrix::getAntiIndentityMatrix(n);

		U[1] = Matrix(c, n);
		W = Matrix(c, s);

		init = _init;

		cluster.resize(c);

		initialization = 0;
		melhorJ = DBL_MAX;

		printf("contrutor OK!\n"); fflush(stdout);
	}

	virtual ~Algorithm() {

	}

	void atulizarResultados(double j, double cr) {
		if (cmp(j, melhorJ) < 0) {
			melhorInicializacao = initialization;
			melhorJ = j;
			melhorCR = cr;
			melhorU = U[pos];
			melhorW = W;
			melhorCluster = cluster;
		}
	}

	void begin() {
		initialization++;
		iteration = 0;

		U[0] = Matrix::getRandInitialU(c, n);
		pos = 0;
		beta = 0;

		for (int i = 0; i < c; i++) {
			for (int j = 0; j < s; j++) {
				W(i, j) = 1.0 / s;
			}
		}
	}

	bool loop() {
		iteration++;

		Matrix R_[c];
		for (int i = 0; i < c; i++) {
			R_[i] = Matrix(n, n);
		}
		Matrix UpowM = U[pos] ^ m;
		Matrix WpowQ = W ^ q;

		// TODO
		// equation [5]
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				for (int cc = 0; cc < c; cc++) {
					R_[cc](i, j) = 0.0;

					for (int ss = 0; ss < s; ss++) {
						R_[cc](i, j) += (WpowQ(cc, ss) * R[ss](i, j));
					}
				}
			}
		}

		// equation [8]
		Matrix Rbeta_[c];
		for (int i = 0; i < c; i++) {
			Rbeta_[i] = R_[i] + (M - I) * beta;
		}

		// equation [6]
		Matrix v[c];
		for (int i = 0; i < c; i++) {
			double div = 0.0;
			for (int j = 0; j < n; j++) {
				div += UpowM(i, j);
			}
			assert(cmp(div) != 0);
			v[i] = ((UpowM.getRow(i)) ^ 'T') / div;
		}

		// equation [5]
		Matrix d(c, n);
		Matrix A, B;
		for (int i = 0; i < c; i++) {
			A = Rbeta_[i] * v[i];
			B = ((v[i] ^ 'T') * (Rbeta_[i] * v[i])) * 0.5;

			assert(B.getNumRow() == 1 && B.getNumColumn() == 1);

			for (int j = 0; j < n; j++) {
				d(i, j) = A(j, 0) - B(0, 0);
			}
		}

		//----------------------------------------------------------------------------
		// TODO
		double delta = 0.0;
		for (int i = 0; i < c; i++) {
			for (int j = 0; j < n; j++) {
				if (cmp(d(i, j)) < 0) {
					double norm =
							Matrix::getEuclidianNorm(v[i], I.getColumn(j));
					assert(cmp(norm) != 0);

					delta = max(delta, -2 * d(i, j) / norm);
				}
			}
		}

		beta += delta;
		delta /= 2;

		if (cmp(delta) > 0) {
			for (int i = 0; i < c; i++) {
				for (int j = 0; j < n; j++) {
					d(i, j) += delta * Matrix::getEuclidianNorm(v[i],
							I.getColumn(j));
				}
			}
		}
		//----------------------------------------------------------------------------

		// equation [7]
		for (int i = 0; i < c; i++) {
			for (int j = 0; j < n; j++) {
				if (cmp(d(i, j)) > 0) {
					double val = 0.0;
					for (int k = 0; k < c; k++) {
						if(cmp(d(k,j)) == 0) {
							return false;
						}
						val += std::pow(d(i, j) / d(k, j), 1.0 / (m - 1.0));
					}

					assert(cmp(val)> 0);
					U[prox(pos)](i, j) = 1.0 / val;
				} else {
					U[prox(pos)](i, j) = 0.0;
				}
			}
		}

#ifndef NDEBUG
		//double j1 = calculaJ();
#endif

		pos = prox(pos);

		for (int i = 0; i < c; i++) {
			cluster[i].clear();
		}

		int group;
		double dist;
		for (int i = 0; i < n; i++) {
			dist = DBL_MIN;
			for (int j = 0; j < c; j++) {
				if (cmp(U[pos](j, i), dist) > 0) {
					dist = U[pos](j, i);
					group = j;
				}
			}
			cluster[group].add(i);
		}

		// equation [22]
		Matrix D(c, s);
		for (int cc = 0; cc < c; cc++) {
			for (int ss = 0; ss < s; ss++) {
				D(cc, ss) = 0;
				for (int j = 0; j < n; j++) {
					for (int k = 0; k < n; k++) {
						D(cc, ss)
								+= (UpowM(cc, j) * UpowM(cc, k) * R[ss](j, k));
					}
				}
			}
		}

		for (int cc = 0; cc < c; cc++) {
			for (int ss = 0; ss < s; ss++) {
				double val = 0.0;
				for (int pp = 0; pp < s; pp++) {
					val += std::pow(D(cc, ss) / D(cc, pp), 1.0 / (q - 1));
				}
				//assert( cmp(val) );
				if (cmp(val)) {
					W(cc, ss) = 1.0 / val;
				} else {
					return false;
				}
			}
		}
#ifndef NDEBUG
		//double j2 = calculaJ();
#endif

#ifndef NDEBUG
		//assert(cmp(j1,j2) >= 0);
#endif

		return true;
	}

	bool repeat() {
		Matrix dif = (U[pos] - U[prox(pos)]).abs();
		dif.abs();

		double maior = dif.getMaxElement();

		cout << iteration << ": [" << scientific << maior << "]" << endl;
		return cmp(maior);
	}

	void printMelhorResultado() {
		out << "\n\n# best Initialization: " << melhorInicializacao << "\n";
		out << "# best J: " << scientific << melhorJ << "\n";
		out << "# CR: " << fixed << setprecision(6) << melhorCR << "\n\n";

		out << "# Algorithm Solution:\n";
		for (int i = 0; i < c; i++) {
			out << "--- Cluster [" << (i + 1) << "]  ("
					<< melhorCluster[i].size() << " elements)\n";
			out << melhorCluster[i] << "\n\n";
		}

		out << "\n";

		out << "# Priori Cluster:\n";
		for (int i = 0; i < (int) prioriCluster.size(); i++) {
			out << "--- Cluster [" << (i + 1) << "]  ("
					<< prioriCluster[i].size() << " elements)\n";
			out << prioriCluster[i] << "\n\n";
		}

		out << "\n";

		out << "# Matrix U:\n";
		out << melhorU << "\n";

		out << "\n";

		if(var_classe != -1) {
			out << "# Matrix W:\n";
			out << melhorW << "\n";
			out << "\n";
			imprimirMatriz(out);
		}

		out << "\n----------------------------X----------------------------\n";
	}

	void process() {
		execute();
	}

	void execute() {
		out << "Results:\n\n";

		out << setw(30) << right << "File(s): " << files[0] << "\n";
		for (int i = 1; i < (int) files.size(); i++) {
			out << setw(30) << right << "         " << files[i] << "\n";
		}

		out << setw(30) << right << "Clusters: " << c << "\n";
		out << setw(30) << right << "Fuzzification parameter: " << m << "\n";
		out << setw(30) << right << "Initialization number: " << init << "\n";
		out << setprecision(6) << "\n";

		double J, CR;
		while (initialization < init) {

			out.flush();
			printf("init %d\n",initialization); fflush(stdout);

			begin();

			printLinha1(out);
			printInicializacaoHeader(out, initialization);
			printLinha2(out);

			print2(out);

			do {
				if (!loop()) {
					break;
				}

				printf("iter %d\n",iteration); fflush(stdout);

				J = calculaJ();
				CR = calculaCR();

				printLinha2(out);
				print3(out, iteration, J, CR);
			} while (repeat() && iteration < r);

			atulizarResultados(J, CR);

			printLinha1(out);
			out << "\n";
		}

		printMelhorResultado();
	}

	double calculaJ() {
		double J = 0.0;

		Matrix UpowM = U[pos] ^ m;
		Matrix WpowQ = W ^ q;

		double num, den;
		for (int i = 0; i < c; i++) {
			num = 0.0;
			for (int j = 0; j < n; j++) {
				for (int k = 0; k < n; k++) {
					double val = 0.0;
					for (int ss = 0; ss < s; ss++) {
						val += (WpowQ(i, ss) * R[ss](j, k));
					}
					num += (UpowM(i, j) * UpowM(i, k) * val);
				}
			}

			den = 0.0;
			for (int j = 0; j < n; j++) {
				den += UpowM(i, j);
			}
			den *= 2;

			assert(cmp(den));
			J += num / den;
		}

		return J;
	}

	double calculaCR() {

		if(var_classe == -1) {
			return 0;
		}

		int k = c;
		int p = prioriCluster.size();

		int table[k + 1][p + 1];

		for (int i = 0; i < k; i++) {
			for (int j = 0; j < p; j++) {
				table[i][j] = cluster[i] && prioriCluster[j];
			}
		}

		for (int i = 0; i < k; i++) {
			table[i][p] = cluster[i].size();
		}

		for (int i = 0; i < p; i++) {
			table[k][i] = prioriCluster[i].size();
		}

		const double pot = 1.0 / (comb(n));
		double termo[4] = { 0 };
		double temp[2] = { 0 };

		// termo 1
		for (int i = 0; i < k; i++) {
			for (int j = 0; j < p; j++) {
				termo[1] += comb(table[i][j]);
			}
		}

		// temps
		for (int i = 0; i < k; i++) {
			temp[0] += comb(table[i][p]);
		}
		for (int i = 0; i < p; i++) {
			temp[1] += comb(table[k][i]);
		}

		// termo 2
		termo[2] = pot * (temp[0] * temp[1]);

		// termo 3
		termo[3] = 0.5 * (temp[0] + temp[1]);

		return (termo[1] - termo[2]) / (termo[3] - termo[2]);
	}

	void imprimirMatriz(ostream &out) {
		double calculaErro(int individuos,vector<vector<int> > &table);
		double fMeasure(vector<vector<int> > &table);

		int k = (int) melhorCluster.size();
		int p = prioriCluster.size();

		//if (k != p) {
		//	return;
		//}

		vector<vector<int> > table;
		table.resize(k + 1);
		for (int i = 0; i <= k; i++) {
			table[i].resize(p + 1);
		}

		// table[i][j] = cluster[i] && prioriCluster[j]
		for (int i = 0; i < k; i++) {
			for (int j = 0; j < p; j++) {
				table[i][j] = melhorCluster[i] && prioriCluster[j];
			}
		}

		// table[i][p] = melhorCluster[i].tamanho
		for (int i = 0; i < k; i++) {
			table[i][p] = melhorCluster[i].size();
		}

		// table[k][i] = prioriCluster[i].tamanho
		for (int i = 0; i < p; i++) {
			table[k][i] = prioriCluster[i].size();
		}

		table[k][p] = this->n;

		//--------------------------------------------------------------------------------
		out << "--- Confusing Matrix ------------------------------\n";

		for (int i = 0; i < p; i++) {
			out << '\t' << i;
		}
		out << '\n';

		for (int i = 0; i < k; i++) {
			out << i << '\'';

			for (int j = 0; j < p; j++) {
				out << '\t' << table[i][j];
			}
			out << '\n';
		}

		//--------------------------------------------------------------------------------
		out << "\n# error: " << fixed << setprecision(2) << (100 * calculaErro(n, table)) << "%\n";
		out << "# F measure: " << fixed << setprecision(6) << fMeasure(table) << "\n";
		//--------------------------------------------------------------------------------

	// more:
	Array< Array<double> > PRIORI;
	PRIORI.resize(this->n);
	for (size_t i = 0; i < this->n; i++) {
		PRIORI[i].resize(p);
		for(int k=0;k<p;k++) {
		  PRIORI[i][k] = prioriCluster[k].contem(i) ? 1.0 : 0.0;
		}
	}
	
	// U[1] = Matrix(c, n);
	
	Array< Array<double> > MINE;
	MINE.resize(this->n);
	for (size_t i = 0; i < this->n; i++) {
		MINE[i].resize(k);
		for(int p=0;p<k;p++) {
		  MINE[i][p] = melhorU(p,i);
		}
	}
	
	double campello = Validation::fuzzy_rand_index_campello(PRIORI,MINE);
	double hullermeier = Validation::fuzzy_rand_index_hullermeier(PRIORI,MINE);
	
	// frigui-campello
	out << "# fuzzy_rand_index_campello: " << fixed << setprecision(6) << campello << "\n";
	// hullermeier
  out << "# fuzzy_rand_index_hullermeier: " << fixed << setprecision(6) << hullermeier << "\n";

	}
};

double calculaErro(int individuos, vector<vector<int> > &table) {

	//MinCostFlow &mcf = Algoritmo::mcf;

	int N = table.size() - 1;
	int M = table[0].size() - 1;

	int source = 0, sink = N + M + 1;

	mcf.init(N + M + 2);

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			mcf.add_edge(i + 1, j + N + 1, 1, -table[i][j]);
		}
	}

	for (int i = 0; i < N; i++) {
		mcf.add_edge(source, i + 1, 1, 0);
	}

	if (N == M) {
		for (int i = 0; i < M; i++) {
			mcf.add_edge(i + N + 1, sink, N, 0);
		}
	} else {
		for (int i = 0; i < M; i++) {
			mcf.add_edge(i + N + 1, sink, N, 0);
		}
	}

	pair<int, int> ret = mcf.solve(source, sink);
	assert(ret.first == N);

	//dbg(ret.second);

	return 1 + double(ret.second) / individuos;
}

double fMeasure(vector<vector<int> > &table) {
	double F = 0;

	int K = table.size() - 1;
	int P = table[0].size() - 1;

	for(int j = 0; j < P; j++) {
		double vmax = 0, rappel, precision;
		for(int i = 0; i < K; i++) if(table[i][j] != 0) {
			rappel = double(table[i][j]) / table[K][j];
			precision = double(table[i][j]) / table[i][P];
			vmax = max(vmax, 2 * rappel * precision / (rappel + precision));
			//printf("rappel = %d/%d, precision = %d/%d\n",table[i][j],table[K][j],table[i][j],table[i][P]);
		}
		F  += vmax * table[K][j];
		//printf(" %lf * %d\n",vmax, table[K][j]);
	}
	//printf("--> indi = %d\n",table[K][P]);
	return F / table[K][P];
}

#undef prox
#undef comb

void process() {
	int num_cluster;
	int m_parameter;
	int q_parameter;
	int num_initialization;
	int num_max_iteration;
	vector<string> in;
	char out[1000];

	freopen("config.txt", "r", stdin);
	char linha[1000];

	gets(linha);
	sscanf(linha, "%d", &num_cluster);
	gets(linha);
	sscanf(linha, "%d", &m_parameter);
	gets(linha);
	sscanf(linha, "%d", &q_parameter);
	gets(linha);
	sscanf(linha, "%d", &num_initialization);
	gets(linha);
	sscanf(linha, "%d", &var_classe);
	gets(linha);
	sscanf(linha, "%d", &num_max_iteration);

	gets(linha);
	istringstream str(linha);
	string tmp;
	while (str >> tmp && tmp.compare(":")) {
		if (tmp.compare(",")) {
			in.push_back(tmp);
			cerr << (tmp) << endl;
		}
	}

	gets(linha);
	sscanf(linha, "%s", out);

	ofstream saida(out, ios::out);
	Algorithm algorithm(in, num_cluster, m_parameter, q_parameter,
			num_initialization, num_max_iteration, saida);

	printf("processando...\n"); fflush(stdout);

	algorithm.process();
}

int main() {
	srand(time(NULL));
	process();
	return 0;
}
