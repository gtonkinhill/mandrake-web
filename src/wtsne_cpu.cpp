/*
 ============================================================================
 Author      : Zhirong Yang
 Copyright   : Copyright by Zhirong Yang. All rights are reserved.
 Description : Stochastic Optimization for t-SNE
 ============================================================================
 */

// Modified by John Lees
#include <iterator>
#include "wtsne.hpp"

std::vector<double> wtsne(const std::vector<int> &I,
                          const std::vector<int> &J,
                          std::vector<double> &dists,
                          const int nsamples,
                          const double perplexity,
                          const int maxIter, const int nRepuSamp,
                          const double eta0, const bool bInit,
                          const int n_threads, const unsigned int seed) {
  // Check input
  std::vector<double> Y, P;
  std::vector<double> weights(nsamples, 1.0);
  std::tie(Y, P) =
      wtsne_init<double>(I, J, dists, weights, perplexity, n_threads, seed);
  long long nn = weights.size();
  long long ne = P.size();

  // Set up random number generation
  std::random_device rd;
  std::mt19937 gen(rd());
  std::discrete_distribution<> dist_de(P.begin(), P.end());
  std::discrete_distribution<> dist_dn(weights.begin(), weights.end());

  // SNE algorithm
  const double nsq = nn * (nn - 1);
  double Eq = 1.0;
  for (long long iter = 0; iter < maxIter; iter++)
  {
    double eta = eta0 * (1 - (double)iter / maxIter);
    eta = MAX(eta, eta0 * 1e-4);
    double c = 1.0 / (Eq * nsq);

    double qsum = 0;
    long long qcount = 0;

    double attrCoef = (bInit && iter < maxIter / 10) ? 8 : 2;
    double repuCoef = 2 * c / nRepuSamp * nsq;
#pragma omp parallel for reduction(+ : qsum, qcount) num_threads(n_threads)
    for (int worker = 0; worker < n_threads; worker++) {
      std::vector<double> dY(DIM);
      std::vector<double> Yk_read(DIM);
      std::vector<double> Yl_read(DIM);

      // long long e = gsl_ran_discrete(gsl_r_ne, gsl_de) % ne;
      long long e = dist_de(gen) % ne;
      // fprintf(stderr, "e: %d", e);
      long long i = I[e];
      long long j = J[e];

      for (long long r = 0; r < nRepuSamp + 1; r++)
      {
        // fprintf(stderr, "r: %d", r);
        // fflush(stderr);
        long long k, l, k2;
        if (r == 0)
        {
          k = i;
          l = j;
        }
        else
        {
          k = dist_dn(gen) % nn;
          l = dist_dn(gen) % nn;
        }
        if (k == l)
          continue;

        long long lk = k * DIM;
        long long ll = l * DIM;
        double dist2 = 0.0;
        for (long long d = 0; d < DIM; d++)
        {
#pragma omp atomic read
          Yk_read[d] = Y[d + lk];
#pragma omp atomic read
          Yl_read[d] = Y[d + ll];
          dY[d] = Yk_read[d] - Yl_read[d];
          dist2 += dY[d] * dY[d];
        }
        double q = 1.0 / (1 + dist2);

        double g;
        if (r == 0)
          g = -attrCoef * q;
        else
          g = repuCoef * q * q;

        bool overwrite = false;
        for (long long d = 0; d < DIM; d++)
        {
          double gain = eta * g * dY[d];
          double Yk_read_end, Yl_read_end;
#pragma omp atomic capture
          Yk_read_end = Y[d + lk] += gain;
#pragma omp atomic capture
          Yl_read_end = Y[d + ll] -= gain;
          if (Yk_read_end != Yk_read[d] + gain || Yl_read_end != Yl_read[d] - gain) {
            overwrite = true;
            break;
          }
        }
        if (!overwrite)
        {
          qsum += q;
          qcount++;
        } else {
          // Find another neighbour
          for (int d = 0; d < DIM; d++) {
#pragma atomic write
            Y[d + lk] = Yk_read[d];
#pragma atomic write
            Y[d + ll] = Yl_read[d];
          }
          r--;
        }
      }
    }
    Eq = (Eq * nsq + qsum) / (nsq + qcount);
    update_progress(iter, maxIter, eta, Eq);
  }
  std::cerr << std::endl
            << "Optimizing done" << std::endl;

  return Y;
}
