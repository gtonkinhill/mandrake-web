/*
 * linear_regression.cpp
 * Regresses k-mer lengths and matches
 *
 */
#include <iostream>
#include <math.h>
#include <algorithm>
#include <limits>

#include "reference.hpp"
#include "Eigen/Dense"

const float core_upper = 0;
const float accessory_upper = 0;

std::tuple<float, float> fit_slope(const Eigen::MatrixXf kmers,
                                   const Eigen::VectorXf dists,
                                   Reference *r1,
                                   Reference *r2)
{

  // Store core/accessory in dists, truncating at zero
  float core_dist = 0, accessory_dist = 0;
  static const double tolerance = (5.0 / (r1->sketchsize64() * 64));
  try
  {
    std::vector<Eigen::Index> truncation;
    for(Eigen::Index i=0; i<dists.size(); ++i){
      if(dists(i)>=tolerance){
        truncation.push_back(i);
      } else {
        break;
      }
    }

    Eigen::VectorXf ldists(dists.size());
    for(Eigen::Index i=0; i<dists.size(); ++i){
      ldists(i) = log(dists(i));
    }

    Eigen::MatrixXf slopes;
    if (truncation.size() > 0)
    {
      slopes = kmers(truncation,Eigen::all).bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(ldists(truncation));
    }
    else
    {
      slopes = kmers.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(ldists);
    }

    if (slopes(1) < core_upper)
    {
      core_dist = 1 - exp(slopes(1));
    }
    if (slopes(0) < accessory_upper)
    {
      accessory_dist = 1 - exp(slopes(0));
    }

  }
  catch (const std::exception &e)
  {
    std::cerr << e.what() << std::endl;
    std::cerr << "Fitting k-mer gradient failed, for:" << r1->name() << "vs." << r2->name() << std::endl;
    std::cerr << dists << std::endl;
    std::cerr << std::endl
              << "Check for low quality genomes" << std::endl;
    exit(1);
  }
  return (std::make_tuple(core_dist, accessory_dist));
}


std::tuple<float, float> regress_kmers(Reference *r1,
                                       Reference *r2,
                                       const Eigen::MatrixXf &kmers,
                                       const RandomMC &random)
{
  // Vector of points
  Eigen::VectorXf dists(kmers.rows());
  for (unsigned int i = 0; i < dists.size(); ++i)
  {
    dists(i) = r1->jaccard_dist(*r2, (int)kmers(i, 1), random);
    // std::cout << "jaccard dist A:" << dists(i) << std::endl;
  }
  return (fit_slope(kmers, dists, r1, r2));
}


std::tuple<float, float> regress_kmers(Reference *r1,
                                       Reference *r2,
                                       const Eigen::MatrixXf &kmers,
                                       const std::vector<double> &random)
{
  // Vector of points
  Eigen::VectorXf dists(kmers.rows());
  for (unsigned int i = 0; i < dists.size(); ++i)
  {
    dists(i) = r1->jaccard_dist(*r2, (int)kmers(i, 1), random[i]);
    // std::cout << "jaccard dist B:" << dists(i) << std::endl;
  }
  return (fit_slope(kmers, dists, r1, r2));
}