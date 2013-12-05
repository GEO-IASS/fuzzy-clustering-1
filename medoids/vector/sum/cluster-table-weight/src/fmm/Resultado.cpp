/*
 * Resultado.cpp
 *
 *  Created on: 13/07/2009
 *      Author: Filipe
 */

#include "Resultado.h"

Resultado::Resultado(size_t k, size_t t, size_t n, size_t m,size_t q, size_t s, int id, double J, double CR) {
	for (size_t i = 0; i < k; i++) {
		this->cluster.push_back(Cluster(t,q));
	}
	this->m = m;
	this->s = s;
	this->U.resize(n);
	for (size_t i = 0; i < n; i++) {
		this->U[i].resize(k);
	}
	this->coeficiente.resize(k);
	for (size_t i = 0; i < k; i++) {
		this->coeficiente[i].resize(t);
	}
	init(id, J, CR);

	overallPrototype.resize(t);
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
	size_t Q = this->cluster[0].prototipo[0].size();
	for (size_t i = 0; i < this->cluster.size(); i++) {
		for (size_t t = 0; t < tabela.size(); t++) {
			map<int,int> repete;
			for(size_t q = 0; q < Q; q++) {
				int novoPrototipo = -1;
				do {
					novoPrototipo = (int) (((double) individuos * rand()) / (RAND_MAX+1.0));
					novoPrototipo %= individuos;
				} while(repete[novoPrototipo]++);
				this->cluster[i].prototipo[t][q] = novoPrototipo;
			}
		}
	}

	for (size_t i = 0; i < this->coeficiente.size(); i++) {
		for (size_t j = 0; j < this->coeficiente[i].size(); j++) {
			this->coeficiente[i][j] = 1.0;
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
	size_t Q = this->cluster[0].prototipo[0].size();

	vector<pair<double, int> > v;
	for (size_t i = 0; i < this->cluster.size(); i++) {
		for (size_t t = 0; t < tabela.size(); t++) {
			v.clear();
			for (size_t j = 0; j < individuos; j++) {
				double d = 0;
				//for (Cluster::iterator iter = this->cluster[i].begin(); iter != this->cluster[i].end(); iter++) {
					//int k = *iter;
					// TODO
				for (size_t k = 0; k < individuos; k++) {
					d += (tabela[t](k, j) * pow(this->U[k][i], m) * pow(this->coeficiente[i][t],s));
				}
				v.push_back(mp(d,j));
			}
			sort(v.begin(), v.end(), less<pair<double, int> > ());
			for(size_t q = 0; q < Q; q++) {
				this->cluster[i].prototipo[t][q] = v[q].second;
			}
		}
	}
}

void Resultado::atualizaCoeficiente(const Repositorio& repositorio) {

	const Array<Tabela> &tabela = repositorio.tabela;
	size_t individuos = tabela[0].n;
	size_t Q = this->cluster[0].prototipo[0].size();
	double num = 1, den[tabela.size()];

	for(size_t i = 0; i < this->cluster.size(); i++) {
		num = 1;
		for(size_t j = 0; j < tabela.size(); j++) {
			den[j] = 0;
			for(size_t k = 0; k < individuos; k++) {
				//den[j] += cluster[i].distancia(k,tabela[j]) * pow(this->U[k][i],m) * this->coeficiente[i][j];
				for(size_t q = 0; q < Q; q++) {
					den[j] += tabela[j](k,cluster[i].prototipo[j][q]) * pow(this->U[k][i], m);
				}
			}
			//num *= den[j];
		}
		
		for (size_t j = 0; j < tabela.size(); j++) {
			double tot = 0;
			for (size_t h = 0; h < tabela.size(); h++) {
				tot += pow( den[j] / den[h], 1.0 / (s-1));
			}
			assert(tot);
			this->coeficiente[i][j] = 1.0 / tot;
		}
		/*num = std::pow(num, 1.0 / tabela.size());
		for (size_t j = 0; j < tabela.size(); j++) {
			assert(den[i]);
			this->coeficiente[i][j] = num / den[j];
			dbg(this->coeficiente[i][j]);
		}
		print("---------------------");*/
	}

	//print7("OK");

}

void Resultado::atualizarU(const Repositorio& repositorio) {
	const Array<Tabela> &tabela = repositorio.tabela;
	size_t individuos = tabela[0].n;
	size_t Q = this->cluster[0].prototipo[0].size();
	double dist[this->cluster.size()];
	for (size_t i = 0; i < individuos; i++) {
		for (size_t j = 0; j < this->cluster.size(); j++) {
			dist[j] = 0;
			for (size_t t = 0; t < tabela.size(); t++) {
				//dist[j] += (this->cluster[j].distancia(i, tabela[t]) * this->coeficiente[j][t]); // TODO
				for(size_t q = 0; q < Q; q++) {
					// TODO
					dist[j] += tabela[t](i, this->cluster[j].prototipo[t][q]) * pow(this->coeficiente[j][t],s);
				}
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
	double novoJ = calculaJ(repositorio);
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

void Resultado::atualizaOverallPrototype(const Repositorio& repositorio) {
	const Array<Tabela> &tabela = repositorio.tabela;
	size_t individuos = tabela[0].n;
	for (size_t t = 0; t < tabela.size(); t++) {
		int prototipo = -1;
		double dmin = 1e100;
		for(size_t p = 0; p < individuos; p++) {
			double dist = 0;
			for (size_t i = 0; i < this->cluster.size(); i++) {
				for (size_t j = 0; j < individuos; j++) {
					// TODO
					dist += tabela[t](j, p) * pow(this->U[j][i], m) * pow(this->coeficiente[i][t],s);
				}
			}
			if (cmp(dist, dmin) < 0) {
				dmin = dist;
				prototipo = p;
			}
		}
		assert(prototipo != -1);
		this->overallPrototype[t] = prototipo;
	}
}

double Resultado::calculaT(const Repositorio& repositorio) {
	atualizaOverallPrototype(repositorio);
	return this->calculaT(repositorio, this->overallPrototype);
}

double Resultado::calculaT(const Repositorio& repositorio, Array<int> prototipo) {
	double T = 0;
	for (size_t i = 0; i < this->cluster.size(); i++) {
		T += this->calculaT(repositorio, i, prototipo);
	}
	return T;
}

double Resultado::calculaT(const Repositorio& repositorio, int k, Array<int> prototipo) {
	const Array<Tabela> &tabela = repositorio.tabela;
	double T = 0;
	for (size_t t = 0; t < tabela.size(); t++) {
		T += this->calculaT(repositorio, k, t, prototipo);
	}
	return T;
}

double Resultado::calculaT(const Repositorio& repositorio, int k, int t, Array<int> prototipo) {
	const Array<Tabela> &tabela = repositorio.tabela;
	size_t individuos = tabela[0].n;
	double T = 0;
	for(size_t i = 0; i < individuos; i++) {
		// TODO
		T += tabela[t](i,prototipo[t]) * pow(this->U[i][k],m) * pow(this->coeficiente[k][t],s);
	}
	return T;
}

double Resultado::calculaJ(const Repositorio& repositorio) {
	double J = 0;
	for (size_t i = 0; i < this->cluster.size(); i++) {
		J += this->calculaJ(repositorio, i);
	}
	return J;
}

double Resultado::calculaJ(const Repositorio& repositorio, int k) {
	const Array<Tabela> &tabela = repositorio.tabela;
	double J = 0;
	for (size_t t = 0; t < tabela.size(); t++) {
		J += this->calculaJ(repositorio, k, t);
	}
	return J;
}

double Resultado::calculaJ(const Repositorio& repositorio, int k, int t) {
	const Array<Tabela> &tabela = repositorio.tabela;
	size_t individuos = tabela[0].n;
	size_t Q = this->cluster[0].prototipo[0].size();
	double J = 0;
	for(size_t i = 0; i < individuos; i++) {
		for(size_t q = 0; q < Q; q++) {
			// TODO
			J += tabela[t](i,this->cluster[k].prototipo[t][q]) * pow(this->U[i][k],m) * pow(this->coeficiente[k][t],s);
		}
	}
	return J;
}

double Resultado::calculaQ(const Repositorio& repositorio) {
	return (1.0 - (this->calculaJ(repositorio) / this->calculaT(repositorio)));
}
double Resultado::calculaQ(const Repositorio& repositorio, int k) {
	this->atualizaOverallPrototype(repositorio);
	return (1.0 - (this->calculaJ(repositorio, k) / this->calculaT(repositorio, k, this->overallPrototype)));
}
double Resultado::calculaQ(const Repositorio& repositorio, int k, int t) {
	this->atualizaOverallPrototype(repositorio);
	return (1.0 - (this->calculaJ(repositorio, k, t) / this->calculaT(repositorio, k, t, this->overallPrototype)));
}
