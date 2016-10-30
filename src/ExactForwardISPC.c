extern uniform int8 *uniform seq_data;
extern uniform int8 **uniform seq_locus;
extern uniform int num_seqs;
extern uniform int seq_size;

export void dExactForward_ISPC_st(uniform int t, uniform int L, uniform int N,
                                  uniform double *uniform Pi,
                                  uniform double *uniform mu,
                                  uniform double *uniform rho,
                                  uniform double *uniform alpha) {
  double theta;
  int8 donor_hap, recipient_hap, H;
  uniform int l=0;

  uniform double *uniform f, *uniform fold;
  f    = uniform new uniform double[N];
  fold = uniform new uniform double[N];

  // Locus zero setup
  foreach(recipient = 0 ... N) {
    recipient_hap = (seq_locus[0][recipient/8] >> recipient%8) & 1;
    fold[recipient] = 0.0;

    for(int donor=0; donor<N; ++donor) {
      donor_hap = (seq_locus[0][donor/8] >> donor%8) & 1;
      H = (recipient_hap ^ donor_hap) & 1;
      theta = (H * mu[0]
                 + (1-H) * (1.0 - mu[0]));

      fold[recipient] += alpha[N*donor + recipient] = theta*Pi[N*donor + recipient];

      alpha[N*donor + recipient] = log(alpha[N*donor + recipient]);
    }

    fold[recipient] = -log(fold[recipient]*rho[0]);
  }

  while(l<t) {
    ++l;
    foreach(recipient = 0 ... N) {
      recipient_hap = (seq_locus[l][recipient/8] >> recipient%8) & 1;
      f[recipient] = 0.0;

      for(int donor=0; donor<N; ++donor) {
        donor_hap = (seq_locus[l][donor/8] >> donor%8) & 1;
        H = (recipient_hap ^ donor_hap) & 1;
        theta = (H * mu[l]
                   + (1-H) * (1.0 - mu[l]));

        f[recipient] += alpha[N*donor + recipient] = (theta * Pi[N*donor + recipient]
                                                        + theta * (1.0-rho[l-1]) * exp(alpha[N*donor + recipient] + fold[recipient]));

        alpha[N*donor + recipient] = log(alpha[N*donor + recipient]) - fold[recipient];
      }

      fold[recipient] = -(log(f[recipient] * rho[l]) - fold[recipient]);
    }
  }

  delete[] f;
  delete[] fold;
}
