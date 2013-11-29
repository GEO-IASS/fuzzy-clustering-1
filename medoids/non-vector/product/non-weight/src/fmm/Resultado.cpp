/*
 * Resultado.cpp
 *
 *  Created on: 13/07/2009
 *      Author: Filipe
 */

#include "Resultado.h"

Resultado::Resultado(size_t k, size_t p, size_t t, size_t n, size_t m, int id, double J, double CR) {
	for (size_t i = 0; i < k; i++) {
		this->cluster.push_back(Cluster(p));
	}
	this->m = m;
	this->U.resize(n);
	for (size_t i = 0; i < n; i++) {
		this->U[i].resize(k);
	}
	init(id, J, CR);
}

Resultado::~Resultado() {
}

void Resultado::init(size_t id, double J, double CR) {
	this->id = id;
	this->J = J;
	this->CR = CR;
}

void Resultado::clear() {
	for (size_t i = 0; i < this->cluster.size(); i++) {
		this->cluster[i].clear();
	}
}

void Resultado::srand(const Repositorio& repositorio) {
	const Array<Tabela> &tabela = repositorio.tabela;
	size_t individuos = tabela[0].n;
	size_t prototipos = this->cluster[0].prototipo.size();

	if (this->cluster.size() * prototipos > individuos) {
		__throw_invalid_argument("Quantidade invalida de prototipos");
	}
	map<int, int> repete;
	for (size_t i = 0; i < this->cluster.size(); i++) {
		this->cluster[i].prototipo.clear();
		for (size_t j = 0; j < prototipos; j++) {
			int novoPrototipo;
			do {
				novoPrototipo = (int) (((double) individuos * rand()) / (RAND_MAX+1.0));
				novoPrototipo %= individuos;
			} while (repete[novoPrototipo]++);
			this->cluster[i].prototipo.push_back(novoPrototipo);
		}
	}

	// inicializa matriz U
	this->atualizarU(repositorio);

	// inicializa J
	this->J = DBL_MAX;
	this->atualizaJ(repositorio);

	// inicializa CR
	this->CR = 0;
	this->atualizaCR(repositorio);
}

void Resultado::atualizaCluster(const Repositorio& repositorio) {
	const Array<Tabela> &tabela = repositorio.tabela;
	size_t individuos = tabela[0].n;
	for (size_t i = 0; i < this->cluster.size(); i++) {
		vector<pair<double, int> > v;
		for (size_t j = 0; j < individuos; j++) {
			double d = 0;
			for (size_t k = 0; k < individuos; k++) {
				for (size_t t = 0; t < tabela.size(); t++) {
					d += (tabela[t](k, j) * pow(this->U[k][i], m));
				}
			}
			v.push_back(mp(d,j));
		}
		sort(v.begin(), v.end(), less<pair<double, int> > ());
		for (size_t j = 0; j < this->cluster[i].prototipo.size(); j++) {
			this->cluster[i].prototipo[j] = v[j].second;
		}
	}
}

void Resultado::atualizarU(const Repositorio& repositorio) {
	const Array<Tabela> &tabela = repositorio.tabela;
	size_t individuos = tabela[0].n;
	double dist[this->cluster.size()];
	for (size_t i = 0; i < individuos; i++) {
		for (size_t j = 0; j < this->cluster.size(); j++) {
			dist[j] = 0;
			for (size_t t = 0; t < tabela.size(); t++) {
				dist[j] += this->cluster[j].distancia(i, tabela[t]);
			}
		}

		vector<size_t> v;
		for (size_t j = 0; j < this->cluster.size(); j++) {
			if (cmp(dist[j]) == 0) {
				v.push_back(j);
			}
		}

		if (v.size()) {
			for (size_t j = 0; j < this->cluster.size(); j++) {
				this->U[i][j] = 0;
			}
			for (size_t j = 0; j < v.size(); j++) {
				this->U[i][v[j]] = (1.0 / v.size());
			}
		} else {
			for (size_t j = 0; j < this->cluster.size(); j++) {
				this->U[i][j] = 0;
				for (size_t k = 0; k < this->cluster.size(); k++) {
					this->U[i][j] += pow(dist[j] / dist[k], 1.0 / (m - 1));
				}
				assert(cmp(this->U[i][j]));
				this->U[i][j] = 1.0 / this->U[i][j];
			}
		}
	}

	this->clear();
	for (size_t i = 0; i < individuos; i++) {
		vector<pair<double, int> > v;
		for (size_t j = 0; j < this->cluster.size(); j++) {
			v.push_back(mp(this->U[i][j],j));
		}
		sort(v.begin(), v.end(), greater<pair<double, int> > ());
		this->cluster[v[0].second].insert(i);
	}
}

double Resultado::atualizaJ(const Repositorio& repositorio) {
	const Array<Tabela> &tabela = repositorio.tabela;
	size_t individuos = tabela[0].n;
	double novoJ = 0;
	for (size_t i = 0; i < this->cluster.size(); i++) {
		for (size_t j = 0; j < individuos; j++) {
			double dist = 0;
			for (size_t t = 0; t < tabela.size(); t++) {
				dist += (this->cluster[i].distancia(j, tabela[t]));
			}
			novoJ += ((dist * pow(this->U[j][i], m)));
		}
	}

	assert(cmp(novoJ,this->J) <= 0);
	return (this->J = novoJ);
}

double Resultado::atualizaCR(const Repositorio& repositorio) {
	if (repositorio.rotulado) {
		this->CR = repositorio.dados.calculaCR(this->cluster);
	}
	return this->CR;
}

Cluster& Resultado::operator[](size_t __n) {
	return this->cluster[__n];
}

const Cluster& Resultado::operator[](size_t __n) const {
	return this->cluster[__n];
}

double& Resultado::operator()(size_t __n, size_t __m) {
	return this->U[__n][__m];
}

const double& Resultado::operator()(size_t __n, size_t __m) const {
	return this->U[__n][__m];
}

bool Resultado::operator<(const Resultado& r) const {
	return (cmp(this->J, r.J) < 0);
}

Resultado::operator string() const {
	os out;
	_dbg(out,id);
	_dbg(out,J);
	_dbg(out,CR);
	_dbg(out,cluster.size());
	for (size_t i = 0; i < cluster.size(); i++) {
		_dbg(out,i);
		_print(out,cluster[i]);
	}
	return out.str();
}

ostream& operator<<(ostream& out, const Resultado& r) {
	return out << (string(r));
}
