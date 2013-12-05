/*
 * Main.cpp
 *
 *  Created on: 04/08/2009
 *      Author: Filipe
 */

/*************************************************************************/
//Vari�veis Globais
int numeroClassesPriori; //N�mero de classes a priori
int numeroClasses; //N�mero de classes na parti��o
double m = 2; //Par�metro de nebulosidade (maior que 1)
int limiteIteracao = 350; //N�mero m�ximo de itera��es
double epsilon = 1e-10; //Diferen�a m�nima entre os crit�rios W
int numeroVariaveis; //N�mero de vari�veis
int totalPadroes; //N�mero total de padr�es
int numeroInicializacoes = 60; //N�mero de inicializa��es da parti��o para cada conjunto
int numeroMonteCarlo = 60; //N�mero de ciclos do Monte Carlo
int numeroClassesEscolhidas = 0; //Numero de classes que foram selecionadas para a penalizacao
/*************************************************************************/

#include "fmm/Algoritmo.h"

vector<string> v;
string arq, saida;
int numCluster;
int numInicializacao;
int var_classe;
int parametro_m;
int parametro_q;
int parametro_s;
int numIteracoes;

void read(const char *nome) {
	ifstream in(nome, ios::in);

	string s;
	bool input = false, output = false;
	while (in >> s) {
		if (s.find("(numCluster)") != string::npos) {
			in >> s;
			sscanf(s.c_str(), "%d", &numCluster);
			input = output = false;
		} else if (s.find("(numInicializacao)") != string::npos) {
			in >> s;
			sscanf(s.c_str(), "%d", &numInicializacao);
			input = output = false;
		} else if (s.find("(variavel_classe)") != string::npos) {
			in >> s;
			sscanf(s.c_str(), "%d", &var_classe);
			Dados::var_classe = var_classe;
			input = output = false;
		}
		else if (s.find("(parametro_m)") != string::npos) {
			in >> s;
			sscanf(s.c_str(), "%d", &parametro_m);
			input = output = false;
		}
		else if (s.find("(parametro_q)") != string::npos) {
			in >> s;
			sscanf(s.c_str(), "%d", &parametro_q);
			input = output = false;
		}
		else if (s.find("(parametro_s)") != string::npos) {
			in >> s;
			sscanf(s.c_str(), "%d", &parametro_s);
			input = output = false;
		}
		else if (s.find("(numIteracoes)") != string::npos) {
			in >> s;
			sscanf(s.c_str(), "%d", &numIteracoes);
			input = output = false;
		} else if (s.find("(input)") != string::npos) {
			input = true;
			output = false;
		} else if (s.find("(output)") != string::npos) {
			input = false;
			output = true;
		} else if (input) {
			v.push_back(s);
		} else if (output) {
			saida = s;
			Algoritmo::saida = s;
		} else {
			assert(0);
		}
	}

	dbg7(numCluster);
	dbg7(numInicializacao);
	dbg7(numIteracoes);
	dbg7(parametro_m);
	dbg7(parametro_q);
	dbg7(parametro_s);
	dbg7(var_classe);

	print7("input:");
	fr(i,0,(int)v.size()) print7((i+1) << " : " << v[i]);
	print7("output:");
	print7(saida);

}

void wait(int seconds) {
	clock_t endwait;
	endwait = clock() + seconds * CLOCKS_PER_SEC;
	while (clock() < endwait);
}

int main(int argc, char **argv) {

	cout << "Digite o nome do arquivo de configuracao: " << endl;
	cin >> arq;

	read(arq.c_str());

	srand(time(NULL));

	try {
		Repositorio repositorio(v);
		ofstream out(saida.c_str(), ios::out);
		Algoritmo algoritmo(numInicializacao, numCluster, parametro_m,parametro_q,parametro_s, numIteracoes, repositorio, out);
		algoritmo.executar();
	} catch (exception &error) {
		cout << error.what() << endl;
	} catch (...) {
		cout << "eita esse erro foi muito do mau!!\n";
	}

	print7("ending...");
	cout.flush();
	wait(2);

	return 0;
}
