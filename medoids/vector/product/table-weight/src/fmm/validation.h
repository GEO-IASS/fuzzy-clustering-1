#ifndef VALIDATION_H_
#define VALIDATION_H_

#include "normalizedmetrics.h"

namespace Validation {

  inline unsigned int comb(unsigned int x) {
    return x * (x - 1) / 2;
  }

  inline Array< Array<double> > psi(Array< Array<double> > U) {
    unsigned int N = U.size(), C = U[0].size();
    Array< Array<double> > R;R.resize(N);
    for(unsigned int i = 0; i < N; ++i) {
      R[i].resize(N);
      for(unsigned int j = 0; j < N; ++j) {
        R[i][j]=0;
      }
    }
    for(unsigned int j = 0; j < N; ++j) {
      for(unsigned int k = 0; k < N; ++k) {
        for(unsigned int i = 0; i < C; ++i) {
          R[j][k] += U[j][i] * U[k][i];
        }
      }
    }
    return R;
  }

  inline double fuzzy_rand_index_campello(Array< Array<double> > U1, Array< Array<double> > U2) {
    U1 = psi(U1);
    U2 = psi(U2);
    unsigned int N = U1.size();
    double N_SS = 0, N_SD = 0, N_DS = 0, N_DD = 0;
    for(unsigned int j = 0; j < N; ++j) {
      for(unsigned int k = 0; k < N; ++k) {
        N_SS += U1[j][k]*U2[j][k];
        N_SD += U1[j][k]*(1-U2[j][k]);
        N_DS += (1-U1[j][k])*U2[j][k];
        N_DD += (1-U1[j][k])*(1-U2[j][k]);
      }
    }
    if(N_SS + N_SD + N_DS + N_DD) {
      return (N_SS + N_DD) / (N_SS + N_SD + N_DS + N_DD);
    } else {
      return 0.0;
    }
  }

  inline double fuzzy_rand_index_hullermeier(Array< Array<double> > U1, Array< Array<double> > U2, const Metrics::NormalizedMetrics& E = Metrics::L1) {
    unsigned int N = U1.size();
    double FRHR = 0;
    for(unsigned int i = 0; i < N; i++) {
      for(unsigned int j = i + 1; j < N; j++) {      
        FRHR += fabs(E(U1[i],U1[j]) - E(U2[i],U2[j]));
      }
    }
    FRHR /= comb(N);
    return 1.0 - FRHR;
  }

};

#endif
