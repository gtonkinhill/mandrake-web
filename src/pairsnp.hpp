#include <iostream>
#include <stdio.h>
#include <string>
#include <zlib.h>

#include "kseq.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <boost/dynamic_bitset.hpp>

struct SparseDist
{
    std::vector<int> rows;
    std::vector<int> cols;
    std::vector<double> distances;
    std::vector<std::string> seq_names;
};

template <typename T>
std::vector<T> combine_vectors(const std::vector<std::vector<T> > &vec,
                               const size_t len)
{
    std::vector<T> all(len);
    auto all_it = all.begin();
    for (size_t i = 0; i < vec.size(); ++i)
    {
        std::copy(vec[i].cbegin(), vec[i].cend(), all_it);
        all_it += vec[i].size();
    }
    return all;
}

KSEQ_INIT(gzFile, gzread)

SparseDist pairsnp(std::string fasta, int n_threads, int dist, int knn)
{
    // open filename and initialise kseq
    int l;
    gzFile fp = gzopen(fasta.c_str(), "r");
    kseq_t *seq = kseq_init(fp);

    size_t n_seqs = 0;
    size_t seq_length = 0;

    // initialise bitmaps
    std::vector<std::string> seq_names;
    std::vector<boost::dynamic_bitset<> > A_snps;
    std::vector<boost::dynamic_bitset<> > C_snps;
    std::vector<boost::dynamic_bitset<> > G_snps;
    std::vector<boost::dynamic_bitset<> > T_snps;

    while (true)
    {
        l = kseq_read(seq);

        if (l == -1) // end of file
            break;
        if (l == -2)
        {
            throw std::runtime_error("Error reading FASTA!");
        }
        if (l == -3)
        {
            throw std::runtime_error("Error reading FASTA!");
        }

        // check sequence length
        if ((n_seqs > 0) && (seq->seq.l != seq_length))
        {
            throw std::runtime_error(
                "Error reading FASTA, variable sequence lengths!");
        }
        seq_length = seq->seq.l;

        seq_names.push_back(seq->name.s);
        boost::dynamic_bitset<> As(seq_length);
        boost::dynamic_bitset<> Cs(seq_length);
        boost::dynamic_bitset<> Gs(seq_length);
        boost::dynamic_bitset<> Ts(seq_length);

        for (size_t j = 0; j < seq_length; j++)
        {

            seq->seq.s[j] = std::toupper(seq->seq.s[j]);

            switch (seq->seq.s[j])
            {
            case 'A':
                As[j] = 1;
                break;
            case 'C':
                Cs[j] = 1;
                break;
            case 'G':
                Gs[j] = 1;
                break;
            case 'T':
                Ts[j] = 1;
                break;

            // M = A or C
            case 'M':
                As[j] = 1;
                Cs[j] = 1;
                break;

            // R = A or G
            case 'R':
                As[j] = 1;
                Gs[j] = 1;
                break;

            // W = A or T
            case 'W':
                As[j] = 1;
                Ts[j] = 1;
                break;

            // S = C or G
            case 'S':
                Cs[j] = 1;
                Gs[j] = 1;
                break;

            // Y = C or T
            case 'Y':
                Cs[j] = 1;
                Ts[j] = 1;
                break;

            // K = G or T
            case 'K':
                Gs[j] = 1;
                Ts[j] = 1;
                break;

            // V = A,C or G
            case 'V':
                As[j] = 1;
                Cs[j] = 1;
                Gs[j] = 1;
                break;

            // H = A,C or T
            case 'H':
                As[j] = 1;
                Cs[j] = 1;
                Ts[j] = 1;
                break;

            // D = A,G or T
            case 'D':
                As[j] = 1;
                Gs[j] = 1;
                Ts[j] = 1;
                break;

            // B = C,G or T
            case 'B':
                Cs[j] = 1;
                Gs[j] = 1;
                Ts[j] = 1;
                break;

            // N = A,C,G or T
            default:
                As[j] = 1;
                Cs[j] = 1;
                Gs[j] = 1;
                Ts[j] = 1;
                break;
            }
        }
        A_snps.push_back(As);
        C_snps.push_back(Cs);
        G_snps.push_back(Gs);
        T_snps.push_back(Ts);

        n_seqs++;
    }
    kseq_destroy(seq);
    gzclose(fp);

    std::vector<std::vector<int> > rows(n_seqs);
    std::vector<std::vector<int> > cols(n_seqs);
    std::vector<std::vector<double> > distances(n_seqs);
    uint64_t len = 0;

#pragma omp parallel for schedule(dynamic, 5) reduction(+ \
                                                        : len) num_threads(n_threads)
    for (uint64_t i = 0; i < n_seqs; i++)
    {

        std::vector<int> comp_snps(n_seqs);
        boost::dynamic_bitset<> res(seq_length);

        for (uint64_t j = 0; j < n_seqs; j++)
        {

            res = A_snps[i] & A_snps[j];
            res |= C_snps[i] & C_snps[j];
            res |= G_snps[i] & G_snps[j];
            res |= T_snps[i] & T_snps[j];

            comp_snps[j] = seq_length - res.count();
        }

        // if using knn find the distance needed
        if (knn >= 0)
        {
            std::vector<int> s_comp = comp_snps;
            std::sort(s_comp.begin(), s_comp.end());
            dist = s_comp[knn + 1];
        }

        // output distances
        for (size_t j = 0; j < n_seqs; j++)
        {
            if ((dist == -1) || (comp_snps[j] <= dist))
            {
                rows[i].push_back(i);
                cols[i].push_back(j);
                distances[i].push_back(comp_snps[j]);
            }
        }
        len += distances[i].size();

        if (i % 10 == 0)
        {
            fprintf(stdout, "Calculating distances: %.1lf%%\n", (double)i / n_seqs * 100);
            fflush(stdout);
        }
    }

    // Combine the lists from each thread
    SparseDist results;
    results.rows = std::vector<int>(combine_vectors(rows, len));
    results.cols = std::vector<int>(combine_vectors(cols, len));
    results.distances = std::vector<double>(combine_vectors(distances, len));
    results.seq_names = seq_names;

    return results;
}
