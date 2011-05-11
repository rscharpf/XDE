#ifndef INDEX_V2_H
#define INDEX_V2_H



inline int qg2index(const int q,const int g,const int Q,const int G) {
  int index = g * Q + q;
  
  return index;
}


inline int qq2index(const int q1,const int q2,const int Q) {
  int index = q1 * Q + q2;

  return index;
}



inline void makeSigma(std::vector<std::vector<double> > &Sigma,const int Q,
		      const double gamma2,const double *tau2Rho,
		      const double *a,const double *sigma2g,
		      const double *rho) {
  Sigma.resize(Q);
  int q;
  for (q = 0; q < Q; q++) {
    Sigma[q].resize(Q);
    Sigma[q][q] = gamma2 * tau2Rho[q];
    Sigma[q][q] *= exp(a[q] * log(sigma2g[q]));
  }

  int q1,q2;
  for (q1 = 0; q1 < Q; q1++)
    for (q2 = q1 + 1; q2 < Q; q2++) {
      int k = qq2index(q1,q2,Q);
      Sigma[q1][q2] = gamma2 * rho[k];
      Sigma[q1][q2] *= sqrt(tau2Rho[q1] * tau2Rho[q2]);
      Sigma[q1][q2] *= exp(0.5 * (a[q1] * log(sigma2g[q1]) + 
				  a[q2] * log(sigma2g[q2])));

      Sigma[q2][q1] = Sigma[q1][q2];
    }

  return;
}


#endif
