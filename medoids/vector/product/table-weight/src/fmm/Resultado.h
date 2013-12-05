/*
 * Resultado.h
 *
 *  Created on: 13/07/2009
 *      Author: Filipe
 */

#ifndef RESULTADO_H_
#define RESULTADO_H_

#include "Includes.h"
#include "Cluster.h"
#include "Array.h"
#include "Repositorio.h"

struct Resultado {
	int id;
	double J, CR;
	Array<Cluster> cluster;
	Array< Array<double> > U;
	Array< double > coeficiente;
	int m;

	Array<int> overallPrototype;

	Resultado(size_t = 0, size_t = 0, size_t = 0, size_t = 2, size_t = 0, int = 0, double = DBL_MAX, double = 0);
	virtual ~Resultado();

	void init(size_t = 0, double = DBL_MAX, double = 0);
	void clear();
	void srand(const Repositorio&);

	void atualizaCluster(const Repositorio&);
	void atualizaCoeficiente(const Repositorio&);
	void atualizarU(const Repositorio&);
	double atualizaJ(const Repositorio&);
	double atualizaCR(const Repositorio&);

	Cluster& operator[](size_t);
	const Cluster& operator[](size_t) const;

	double& operator()(size_t, size_t);
	const double& operator()(size_t, size_t) const;

	bool operator<(const Resultado&) const;
	operator string() const;

	friend ostream& operator<<(ostream&, const Resultado&);

	// ------------------------------------------------
	// interpretation functions
	void atualizaOverallPrototype(const Repositorio& repositorio);

	double calculaT(const Repositorio& repositorio);
	double calculaT(const Repositorio& repositorio, Array<int> prototipo);
	double calculaT(const Repositorio& repositorio, int k, Array<int> prototipo);
	double calculaT(const Repositorio& repositorio, int k, int t, Array<int> prototipo);

	double calculaJ(const Repositorio& repositorio);
	double calculaJ(const Repositorio& repositorio, int k);
	double calculaJ(const Repositorio& repositorio, int k, int t);

	double calculaQ(const Repositorio& repositorio);
	double calculaQ(const Repositorio& repositorio, int k);
	double calculaQ(const Repositorio& repositorio, int k, int t);
	// ------------------------------------------------
};

#endif /* RESULTADO_H_ */
