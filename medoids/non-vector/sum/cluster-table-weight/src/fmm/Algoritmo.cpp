/*
 * Algoritmo.cpp
 *
 *  Created on: 12/07/2009
 *      Author: Filipe
 */

#include "Algoritmo.h"
#include "validation.h"

MinCostFlow Algoritmo::mcf;
string Algoritmo::saida = "";

Algoritmo::Algoritmo(size_t inicializacoes, size_t clusters, size_t m, size_t q, size_t s, size_t limite, const Repositorio& repositorio, ostream& out) :
	inicializacoes(inicializacoes), clusters(clusters), repositorio(repositorio), out(out) {
	melhor = atual = Resultado(clusters, repositorio.tabela.size(), repositorio.tabela[0].n, m, q, s, 0);
	this->individuos = this->repositorio.tabela[0].n;
	this->limite = limite;
	niter = 0;
}

Algoritmo::~Algoritmo() {
}

void Algoritmo::executar() {


	this->out << "Resultados:\n\n";

	this->out << "# Numero de inicializacoes: " << this->inicializacoes << "\n";
	this->out << "# Numero de individuos: " << this->individuos << "\n";
	this->out << "# Numero de clusters: " << this->clusters << "\n";
	this->out << "# " << repositorio << "\n";

	this->out << "\n\n";

	for (size_t i = 1; i <= this->inicializacoes; i++) {
		this->atual.init(i);
		while (this->inicializacao())
			;

		// ------------------------------------------------
		size_t iteracao = 1;
		bool repete;
		// ------------------------------------------------

		// ------------------------------------------------
		// - IMPRESSAO ------------------------------------
		printLinha1(this->out);
		printInicializacaoHeader(this->out, i);
		printLinha2(this->out);

		print2(this->out);

		printLinha2(this->out);
		print3(this->out, iteracao++, this->atual.J, this->atual.CR);
		// ------------------------------------------------

		do {
			this->etapa1();
			repete = this->etapa2();

			// ------------------------------------------------
			// - IMPRESSAO ------------------------------------
			if (repete) {
				printLinha2(this->out);
				print3(this->out, iteracao++, this->atual.J, this->atual.CR);
			}
			// ------------------------------------------------

		} while (repete && iteracao <= limite);

		niter += iteracao;

		// ------------------------------------------------
		// - IMPRESSAO ------------------------------------
		printLinha1(this->out);
		this->out << "\n";
		// ------------------------------------------------
	}

	// ------------------------------------------------
	// - IMPRESSAO ------------------------------------
	this->out << "\n\n===================================================";
	this->out << "===================================================\n\n";

	// parte 1
	this->out << "# Melhor inicializacao: " << this->melhor.id << "\n";
	this->out << "# Melhor J: " << scientific << this->melhor.J << "\n";
	this->out << "# CR: " << scientific << this->melhor.CR << "\n";

	this->out << "\n";

	// parte 2
	size_t k = this->melhor.cluster.size();
	size_t n = this->repositorio.tabela[0].n;

	size_t c[n], p[n];
	for (size_t i = 0; i < this->melhor.cluster.size(); i++) {
		tr( (this->melhor.cluster[i]) ,iter) {
			c[*iter] = i;
		}
	}
	if (this->repositorio.rotulado) {
		for (size_t i = 0; i < this->repositorio.dados.prioriCluster.size(); i++) {
			tr( (this->repositorio.dados.prioriCluster[i]) ,iter) {
				p[*iter] = i;
			}
		}
	} else {
		memset(p, -1, sizeof(p));
	}

	printLinha1(this->out, k);
	printRelacaoHeader(this->out, k);

	for (size_t i = 0; i < n; i++) {
		printLinha2(this->out, k);
		print2(this->out, i, k, this->melhor.U[i], c[i], p[i]);
	}
	printLinha1(this->out, k);
	this->out << "\n";

	// parte 3
	this->out << "# Solucao:\n";
	for (size_t i = 0; i < this->clusters; i++) {
		this->out << "Cluster[" << (i + 1) << "] :\n";
		this->out << this->melhor[i] << "\n\n";
	}

	if (this->repositorio.rotulado) {
		this->out << "# Priori Cluster:\n";
		for (size_t i = 0; i < this->repositorio.dados.prioriCluster.size(); i++) {
			this->out << "Cluster[" << (i + 1) << "] :\n";
			this->out << this->repositorio.dados.prioriCluster[i] << "\n\n";
		}
	}

	this->out << "# Coeficientes clusters-tabelas:\n";
	this->out << fixed << setprecision(6);
	for (size_t i = 0; i < this->melhor.coeficiente.size(); i++) {
		this->out << "Cluster #" << (i+1) << "\n";
		for(size_t j = 0; j < this->melhor.coeficiente[i].size(); j++) {
			this->out << "\ttabela[" << j << "] : " << this->melhor.coeficiente[i][j] << "\n";
		}
	}
	this->out << "\n";

	// ------------------------------------------------
	// confusing matrix
	if(Dados::var_classe > 0) {
		this->imprimirMatriz(out);
	}

	// ------------------------------------------------
	// access file
	this->printAcessFile(c);

	// ------------------------------------------------
	// interpretation index
	this->out << "\n";
	this->printGlobalInertia();
	this->out << "\n\n";
	this->printWithinClusterInertia();
	this->out << "\n\n";
	this->printGeneralIndex();
	this->out << "\n";
	this->print_estranho();
	this->out << "\n\n";
}

bool Algoritmo::inicializacao() {
	this->atual.clear();
	this->atual.srand(this->repositorio);
	for (size_t i = 0; i < this->atual.cluster.size(); i++) {
		if (this->atual[i].size() == 0) {
			return true;
		}
	}
	print7("~~ ######################################### ~~");
	return false;
}

void Algoritmo::etapa1() {
	this->atual.atualizaCluster(this->repositorio);
	this->atual.atualizaCoeficiente(this->repositorio);
}

bool Algoritmo::etapa2() {

	this->atual.atualizarU(this->repositorio);
	double dif = this->atual.J;
	dif -= this->atual.atualizaJ(this->repositorio);
	this->atual.atualizaCR(this->repositorio);
	if (this->atual < this->melhor) {
		this->melhor = this->atual;
	}
	dbg7(atual.J);
	return (cmp(dif) != 0);
}

Algoritmo::operator string() const {
	os out;
	_dbg(out,this->inicializacoes);
	_dbg(out,this->clusters);
	_print(out,this->repositorio);
	return out.str();
}

void Algoritmo::printAcessFile(size_t c[]) {
	// ---------------------------------------------------------------------
	// arquivo para ser usado pelo acess da microsoft
	string nome_arquivo = Algoritmo::saida;
	for (int i = nome_arquivo.size() - 1; i >= 0; i--) {
		if (nome_arquivo[i] == '.') {
			nome_arquivo.erase(i, nome_arquivo.size() - i);
			break;
		}
	}
	nome_arquivo += "_acess.txt";

	ofstream extra(nome_arquivo.c_str(), ios::out);
	for (size_t i = 0; i < individuos; i++) {
		extra << (i + 1) << " " << (c[i] + 1) << "\n";
	}
	extra.close();
	// ---------------------------------------------------------------------
}

void Algoritmo::imprimirMatriz(ostream &out) {
	double calculaErro(int individuos, vector<vector<int> > &table);
	double fMeasure(vector<vector<int> > &table);

	//int k = (int) melhorCluster.size();
	//int p = dados.getNumPrioriCluster();
	const Array<Cluster>& C = this->melhor.cluster;
	const Array<Cluster>& P = this->repositorio.dados.prioriCluster;

	int k = C.size();
	int p = P.size();

	vector<vector<int> > table;
	table.resize(k + 1);
	for (int i = 0; i <= k; i++) {
		table[i].resize(p + 1);
	}

	// table[i][j] = cluster[i] && prioriCluster[j]
	for (int i = 0; i < k; i++) {
		for (int j = 0; j < p; j++) {
			table[i][j] = (C[i] && P[j]).size();
		}
	}

	// table[i][p] = melhorCluster[i].tamanho
	for (int i = 0; i < k; i++) {
		table[i][p] = C[i].size();
	}

	// table[k][i] = prioriCluster[i].tamanho
	for (int i = 0; i < p; i++) {
		table[k][i] = P[i].size();
	}

	table[k][p] = this->individuos;

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

	error1 = calculaErro(individuos, table);
	fmed = fMeasure(table);

	//--------------------------------------------------------------------------------
	out << "\n# error: " << fixed << setprecision(2) << (100 * error1) << "%\n";
	out << "# F measure: " << fixed << setprecision(6) << fmed << "\n";
	//--------------------------------------------------------------------------------

	// more:
	Array< Array<double> > PRIORI;
	PRIORI.resize(individuos);
	for (size_t i = 0; i < individuos; i++) {
		PRIORI[i].resize(p);
		for(int k=0;k<p;k++) {
		  PRIORI[i][k] = P[k].count(i) ? 1.0 : 0.0;
		}
	}
	double campello = Validation::fuzzy_rand_index_campello(PRIORI,this->melhor.U);
	double hullermeier = Validation::fuzzy_rand_index_hullermeier(PRIORI,this->melhor.U);
	
	// frigui-campello
	out << "# fuzzy_rand_index_campello: " << fixed << setprecision(6) << campello << "\n";
	// hullermeier
  out << "# fuzzy_rand_index_hullermeier: " << fixed << setprecision(6) << hullermeier << "\n";

}

double calculaErro(int individuos, vector<vector<int> > &table) {

	MinCostFlow &mcf = Algoritmo::mcf;

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
			mcf.add_edge(i + N + 1, sink, 1, 0);
		}
	} else {
		for (int i = 0; i < M; i++) {
			mcf.add_edge(i + N + 1, sink, N, 0);
		}
	}

	pair<int, int> ret = mcf.solve(source, sink);
	assert(ret.first == N);

	dbg7(ret.first);
	dbg7(ret.second);

	return 1 + double(ret.second) / individuos;
}

double fMeasure(vector<vector<int> > &table) {
	double F = 0;
	int K = table.size() - 1;
	int P = table[0].size() - 1;

	for (int j = 0; j < P; j++) {
		double vmax = 0, rappel, precision;
		for (int i = 0; i < K; i++)
			if (table[i][j] != 0) {
				rappel = double(table[i][j]) / table[K][j];
				precision = double(table[i][j]) / table[i][P];
				vmax = max(vmax, 2 * rappel * precision / (rappel + precision));
			}
		F += vmax * table[K][j];
	}

	return F / table[K][P];
}


ostream& operator<<(ostream& out, Algoritmo& a) {
	return out << (string(a));
}

void Algoritmo::printGlobalInertia() {
	this->melhor.atualizaOverallPrototype(this->repositorio);
	Array<int> & prototipo = this->melhor.overallPrototype;
	double T = this->melhor.calculaT(this->repositorio, prototipo);

	this->out << fixed << setprecision(1);

	this->out
			<< "/------------------------------  Global Inertia - T ------------------------------/"
			<< "\n";
	this->out << "\n";
	this->out << "# Data set global inertia = " << T << "\n";
	this->out << "\n";

	// ------------------------------------------------------------
	// (item a)
	this->out << "# The global inertia for each cluster" << "\n";
	this->out << "\n";
	this->out << "===========================" << "\n";
	this->out << "CLUSTER               VALUE" << "\n";
	this->out << "---------------------------" << "\n";
	for (size_t k = 0; k < this->melhor.cluster.size(); k++) {
		out << setw(7) << (k + 1);
		out << setw(20) << this->melhor.calculaT(this->repositorio, k,
				prototipo);
		out << "\n";
	}
	this->out << "===========================" << "\n";
	this->out << setw(7) << "Sum";
	this->out << setw(20) << T;
	this->out << "\n";
	// ------------------------------------------------------------

	this->out << "\n";

	// ------------------------------------------------------------
	// (item b)
	this->out << "# The global inertia for each table" << "\n";
	this->out << "\n";
	this->out << "=========================" << "\n";
	this->out << "TABLE               VALUE" << "\n";
	this->out << "-------------------------" << "\n";
	for (size_t t = 0; t < this->repositorio.tabela.size(); t++) {
		double dist = 0;
		for (size_t k = 0; k < this->melhor.cluster.size(); k++) {
			dist += this->melhor.calculaT(this->repositorio, k, t, prototipo);
		}
		this->out << setw(5) << (t + 1);
		this->out << setw(20) << dist;
		this->out << "\n";
	}
	this->out << "=========================" << "\n";
	this->out << setw(5) << "Sum";
	this->out << setw(20) << T;
	this->out << "\n";
	// ------------------------------------------------------------

	this->out << "\n";

	// ------------------------------------------------------------
	// (item c)
	this->out << "# The global inertia for each table and cluster" << "\n";
	this->out << "\n";

	this->out << "============";
	for (size_t k = 0; k < this->melhor.cluster.size(); k++) {
		this->out << "====================";
	}
	this->out << "\n";

	this->out << "            ";
	for (size_t k = 0; k < this->melhor.cluster.size(); k++) {
		this->out << setw(16) << "CLUSTER";
		this->out << setw(4) << (k + 1);
	}
	this->out << "\n";

	this->out << "------------";
	for (size_t k = 0; k < this->melhor.cluster.size(); k++) {
		this->out << "--------------------";
	}
	this->out << "\n";

	for (size_t t = 0; t < this->repositorio.tabela.size(); t++) {
		this->out << "TABLE";
		this->out << setw(7) << (t + 1);
		double T_t = 0;
		for (size_t k = 0; k < this->melhor.cluster.size(); k++) {
			double T_kt = this->melhor.calculaT(this->repositorio, k, t,
					prototipo);
			T_t += T_kt;
			this->out << setw(20) << T_kt;
		}
		this->out << " | " << T_t << "\n";
	}

	this->out << "============";
	for (size_t k = 0; k < this->melhor.cluster.size(); k++) {
		this->out << "====================";
	}
	this->out << "\n";

	this->out << "            ";
	for (size_t k = 0; k < this->melhor.cluster.size(); k++) {
		this->out << setw(20) << this->melhor.calculaT(this->repositorio, k,
				prototipo);
	}
	this->out << "\n";
	// ------------------------------------------------------------
}

void Algoritmo::printWithinClusterInertia() {
	double J = this->melhor.calculaJ(this->repositorio);

	this->out
			<< "/------------------------------  Within-cluster Inertia - J ------------------------------/"
			<< "\n";
	this->out << "\n";
	this->out << "# Data set within-cluster inertia = " << J << "\n";
	this->out << "\n";

	// ------------------------------------------------------------
	// (item a)
	this->out << "# The within-cluster inertia for each cluster" << "\n";
	this->out << "\n";
	this->out << "===========================" << "\n";
	this->out << "CLUSTER               VALUE" << "\n";
	this->out << "---------------------------" << "\n";
	for (size_t k = 0; k < this->melhor.cluster.size(); k++) {
		this->out << setw(7) << (k + 1);
		this->out << setw(20) << this->melhor.calculaJ(this->repositorio, k);
		this->out << "\n";
	}
	this->out << "===========================" << "\n";
	this->out << setw(7) << "Sum";
	this->out << setw(20) << J;
	this->out << "\n";
	// ------------------------------------------------------------

	this->out << "\n";

	// ------------------------------------------------------------
	// (item b)
	this->out << "# The within-cluster inertia for each table" << "\n";
	this->out << "\n";
	this->out << "=========================" << "\n";
	this->out << "TABLE               VALUE" << "\n";
	this->out << "-------------------------" << "\n";
	for (size_t t = 0; t < this->repositorio.tabela.size(); t++) {
		double dist = 0;
		for (size_t k = 0; k < this->melhor.cluster.size(); k++) {
			dist += this->melhor.calculaJ(this->repositorio, k, t);
		}
		this->out << setw(5) << (t + 1);
		this->out << setw(20) << dist;
		this->out << "\n";
	}
	this->out << "=========================" << "\n";
	this->out << setw(5) << "Sum";
	this->out << setw(20) << J;
	this->out << "\n";
	// ------------------------------------------------------------

	this->out << "\n";

	// ------------------------------------------------------------
	// (item c)
	this->out << "# The global inertia for each table and cluster" << "\n";
	this->out << "\n";

	this->out << "============";
	for (size_t k = 0; k < this->melhor.cluster.size(); k++) {
		this->out << "====================";
	}
	this->out << "\n";

	this->out << "            ";
	for (size_t k = 0; k < this->melhor.cluster.size(); k++) {
		this->out << setw(16) << "CLUSTER";
		this->out << setw(4) << (k + 1);
	}
	this->out << "\n";

	this->out << "------------";
	for (size_t k = 0; k < this->melhor.cluster.size(); k++) {
		this->out << "--------------------";
	}
	this->out << "\n";

	for (size_t t = 0; t < this->repositorio.tabela.size(); t++) {
		this->out << "TABLE";
		this->out << setw(7) << (t + 1);
		double J_t = 0;
		for (size_t k = 0; k < this->melhor.cluster.size(); k++) {
			double J_kt = this->melhor.calculaJ(this->repositorio, k, t);
			J_t += J_kt;
			this->out << setw(20) << J_kt;
		}
		this->out << " | " << J_t << "\n";
	}

	this->out << "============";
	for (size_t k = 0; k < this->melhor.cluster.size(); k++) {
		this->out << "====================";
	}
	this->out << "\n";

	this->out << "            ";
	for (size_t k = 0; k < this->melhor.cluster.size(); k++) {
		this->out << setw(20) << this->melhor.calculaJ(this->repositorio, k);
	}
	this->out << "\n";
	// ------------------------------------------------------------

}

void Algoritmo::printGeneralIndex() {
	double Q = this->melhor.calculaQ(this->repositorio);
	this->melhor.atualizaOverallPrototype(this->repositorio);
	Array<int> & prototipo = this->melhor.overallPrototype;

	this->out
			<< "/------------------------------  General Index  ------------------------------/"
			<< "\n";
	this->out << "# Quality of the partition - Q" << "\n";
	this->out << "(The proportion of inertia explains by the partition): ";
	this->out << fixed << setprecision(6) << Q << "\n";
	this->out << "\n";

	this->out << "\n";

	// ------------------------------------------------------------
	this->out << "# Quality of the table j - Q(j)" << "\n";
	this->out << "\n";
	this->out << "=========================" << "\n";
	this->out << "TABLE               VALUE" << "\n";
	this->out << "-------------------------" << "\n";
	double sum = 0;
	for (size_t t = 0; t < this->repositorio.tabela.size(); t++) {
		double J_t = 0;
		double T_t = 0;
		for (size_t k = 0; k < this->melhor.cluster.size(); k++) {
			J_t += this->melhor.calculaJ(this->repositorio, k, t);
			T_t += this->melhor.calculaT(this->repositorio, k, t, prototipo);
		}
		double Q_t = 1 - (J_t / T_t);
		sum += Q_t;
		this->out << setw(5) << (t + 1);
		this->out << setw(20) << Q_t << " | 1 - (" << J_t << "/" << T_t << ")";
		;
		this->out << "\n";
	}
	this->out << "=========================" << "\n";
	this->out << setw(5) << "Sum";
	this->out << setw(20) << sum;
	this->out << "\n";
	// ------------------------------------------------------------

	this->out << "\n";
}

void Algoritmo::print_estranho() {
	this->out
			<< "/------------------------------  Cluster interpretation indices  ------------------------------/"
			<< "\n";
	this->out << "\n";

	// print J's
	double J = this->melhor.calculaJ(this->repositorio);
	this->out
			<< "# The relative contribution of clusters to the within-cluster inertia - J(i)"
			<< "\n";
	this->out << "\n";
	this->out << "===========================" << "\n";
	this->out << "CLUSTER               VALUE" << "\n";
	this->out << "---------------------------" << "\n";
	double sum = 0;
	for (size_t k = 0; k < this->melhor.cluster.size(); k++) {
		double J_k = this->melhor.calculaJ(this->repositorio, k);
		sum += (J_k / J);
		this->out << setw(7) << (k + 1);
		this->out << setw(20) << (J_k / J) << "  |  (" << J_k << "/" << J
				<< ")";
		this->out << "\n";
	}
	this->out << "===========================" << "\n";
	this->out << setw(7) << "Sum";
	this->out << setw(20) << sum;
	this->out << "\n";

	this->out << "\n";

	// print T's
	double T = this->melhor.calculaT(this->repositorio);
	this->melhor.atualizaOverallPrototype(this->repositorio);
	Array<int> & prototipo = this->melhor.overallPrototype;

	this->out
			<< "# The proportion of the global inertia explained by clusters - T(i)"
			<< "\n";
	this->out << "\n";
	this->out << "===========================" << "\n";
	this->out << "CLUSTER               VALUE" << "\n";
	this->out << "---------------------------" << "\n";
	sum = 0;
	for (size_t k = 0; k < this->melhor.cluster.size(); k++) {
		double T_k = this->melhor.calculaT(this->repositorio, k, prototipo);
		sum += (T_k / T);
		this->out << setw(7) << (k + 1);
		this->out << setw(20) << (T_k / T) << "  |  (" << T_k << "/" << T
				<< ")";
		this->out << "\n";
	}
	this->out << "===========================" << "\n";
	this->out << setw(7) << "Sum";
	this->out << setw(20) << sum;
	this->out << "\n";

	this->out << "\n";

	// ------------------------------------------------------------
	this->out << "# Quality of the cluster i - Q(i)" << "\n";
	this->out << "\n";
	this->out << "===========================" << "\n";
	this->out << "CLUSTER               VALUE" << "\n";
	this->out << "---------------------------" << "\n";
	sum = 0;
	for (size_t k = 0; k < this->melhor.cluster.size(); k++) {
		double Q_k = this->melhor.calculaQ(this->repositorio, k);
		sum += Q_k;
		this->out << setw(7) << (k + 1);
		this->out << setw(20) << (Q_k) << " | 1 - (" << (this->melhor.calculaJ(
				this->repositorio, k)) << "/" << (this->melhor.calculaT(
				this->repositorio, k, prototipo)) << ")";
		this->out << "\n";
	}
	this->out << "===========================" << "\n";
	this->out << setw(7) << "Sum";
	this->out << setw(20) << sum;
	this->out << "\n";
	// ------------------------------------------------------------

	this->out << "\n";

	// ------------------------------------------------------------
	this->out << "# Quality of the cluster i for the table j - Q(i)(j) "
			<< "\n";
	this->out << "\n";

	this->out << "============";
	for (size_t k = 0; k < this->melhor.cluster.size(); k++) {
		this->out << "====================";
	}
	this->out << "\n";

	this->out << "            ";
	for (size_t k = 0; k < this->melhor.cluster.size(); k++) {
		this->out << setw(16) << "CLUSTER";
		this->out << setw(4) << (k + 1);
	}
	this->out << "\n";

	this->out << "------------";
	for (size_t k = 0; k < this->melhor.cluster.size(); k++) {
		this->out << "--------------------";
	}
	this->out << "\n";

	for (size_t t = 0; t < this->repositorio.tabela.size(); t++) {
		this->out << "TABLE";
		this->out << setw(7) << (t + 1);
		for (size_t k = 0; k < this->melhor.cluster.size(); k++) {
			this->out << setw(20) << this->melhor.calculaQ(this->repositorio,
					k, t);
		}
		this->out << "\n";
	}

	this->out << "============";
	for (size_t k = 0; k < this->melhor.cluster.size(); k++) {
		this->out << "====================";
	}
	this->out << "\n";
	// ------------------------------------------------------------
}
