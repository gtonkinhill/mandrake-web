#include <emscripten.h>
#include <emscripten/val.h>
#include <emscripten/bind.h>

#include "pairsnp.hpp"
#include "pairgene.hpp"
#include "wtsne_cpu.cpp"

using namespace emscripten;

EMSCRIPTEN_BINDINGS(module) {
  function("pairsnp", &pairsnp);
  function("pairgene", &pairgene);
  function("wtsne", &wtsne);

  value_object<SparseDist>("SparseDist")
  .field("rows", &SparseDist::rows) // Need to register the array type
  .field("cols", &SparseDist::cols) // Need to register the array type
  .field("distances", &SparseDist::distances) // Need to register the array type
  .field("seq_names", &SparseDist::seq_names) // Need to register the array type
  ;

  register_vector<int>("vector<int>");
  register_vector<uint64_t>("vector<uint64_t>");
  register_vector<double>("vector<double>");
  register_vector<std::string>("vector<std::string>");
}

