/*
 * Cluster.cpp
 *
 *  Created on: 12/07/2009
 *      Author: Filipe
 */

#include "Cluster.h"

Cluster::Cluster(size_t n, size_t q) {
	prototipo.resize(n);
	for(size_t i = 0; i < n; i++) {
		prototipo[i].resize(q);
	}
}

Cluster::~Cluster() {
}

Cluster Cluster::operator&&(const Cluster& c) const {
	Cluster intersecao(0,0);
	tr ((*this), iter) {
		if (c.count(*iter)) {
			intersecao.insert(*iter);
		}
	}
	return intersecao;
}

Cluster::operator string() const {
	os out;
	if(this->prototipo.size()) {
		out << "Prototipos(s): {";
		for(size_t j = 0; j < this->prototipo.size(); j++) {
			out << " ";
			if(j) out << ",";
			out << this->prototipo[j];
		}
		out << " }\n";
	}
	out << "Elemento(s): {";
	for (Cluster::iterator iter = begin(); iter != end(); iter++) {
		if (iter != begin()) {
			out << ",";
		}
		out << " " << (*iter);
	}
	out << " }";
	return out.str();
}

ostream& operator<<(ostream& out, const Cluster& c) {
	return out << (string(c));
}
