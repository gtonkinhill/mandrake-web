#include <stdio.h>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <boost/dynamic_bitset.hpp>

SparseDist pairgene(std::string rtab, int n_threads, int dist, int knn)
{
    // first get number of genes
    int ngenes = 0;
    std::string line;
    std::ifstream gf(rtab);

    while (std::getline(gf, line))
        ++ngenes;
    ngenes -= 1;

    // open filename
    size_t n_seqs = 0;

    // initialise bitmaps
    std::vector<std::string> seq_names;
    std::vector<boost::dynamic_bitset<>> allgenes;

    SparseDist results;

    gf.clear();
    gf.seekg(0);

    int linenum=0;
    while (getline(gf, line))
    {
        std::istringstream iss(line);
        std::string token;
        std::vector<std::string> parts;
        while (std::getline(iss, token, '\t'))
        { // but we can specify a different one
            parts.push_back(token);
        }
        if (linenum == 0)
        {
            n_seqs = parts.size() - 1;
            for (size_t j = 0; j < n_seqs; j++)
            {
                boost::dynamic_bitset<> genes(ngenes);
                allgenes.push_back(genes);
                results.seq_names.push_back(parts[j+1]);
            }
        }
        else
        {
            
            for (size_t j = 0; j < n_seqs; j++)
            {
                char character = parts[j+1][0];
                if (character == '1')
                {
                    allgenes[j][linenum-1] = 1;
                }
                else
                {
                    allgenes[j][linenum-1] = 0;
                }
            }
        }
        linenum+=1;
    }
    gf.close();

#pragma omp parallel for ordered shared(allgenes, ngenes, n_seqs, seq_names, dist, knn) default(none) schedule(static, 1) num_threads(n_threads)
    for (size_t i = 0; i < n_seqs; i++)
    {

        std::vector<int> comp_snps(n_seqs);
        boost::dynamic_bitset<> res(ngenes);

        size_t start = 0;
        if (knn < 0)
        {
            start = i + 1;
        }
        else
        {
            start = 0;
        }

        for (size_t j = start; j < n_seqs; j++)
        {
            res = allgenes[i] & allgenes[j];
            comp_snps[j] = ngenes - res.count();
        }

        // if using knn find the distance needed
        if (knn >= 0)
        {
            std::vector<int> s_comp = comp_snps;
            std::sort(s_comp.begin(), s_comp.end());
            dist = s_comp[knn + 1];
            start = 0;
        }
        else
        {
            start = i + 1;
        }

// output distances
#pragma omp critical
        for (size_t j = 0; j < n_seqs; j++)
        {
            if ((dist == -1) || (comp_snps[j] <= dist))
            {
                results.rows.push_back(i);
                results.cols.push_back(j);
                results.distances.push_back(comp_snps[j]);
            }
        }
    }

    return results;
}