#include <utility>
#include <vector>
#include <string>
#include <cstdio>
#include <iostream>
#include "src/main.hxx"

using namespace std;




template <class G, class K, class V>
void getModularity(const G& x, const LouvainResult<K>& a, V M) {
  auto fc = [&](auto u) { return a.membership[u]; };
  return modularity(x, fc, M, V(1));
}


template <class G>
void runLouvain(const G& x, int repeat) {
  using K = typename G::key_type;
  using V = typename G::edge_value_type;
  vector<K> *init = nullptr;
  random_device dev;
  default_random_engine rnd(dev());
  int retries  = 5;
  V resolution = V(1);
  V tolerance  = V(1e-2);
  V passTolerance = V(0);
  V toleranceDeclineFactor = V(10);
  auto M = edgeWeight(x)/2;
  auto Q = modularity(x, M, 1.0f);
  printf("[%01.6f modularity] noop\n", Q);

  // Get last pass community memberships (static).
  LouvainResult<K> al = louvainSeq(x, init, {repeat, resolution, tolerance, passTolerance, toleranceDeclineFactor});
  printf("[%1.1e batch_size; %d batch_count; %09.3f ms; %04d iters.; %03d passes; %01.9f modularity] louvainSeqLast\n", 0.0, 0, al.time, al.iterations, al.passes, getModularity(x, al, M), tolerance, toleranceDeclineFactor);
  // Get first pass community memberships (static).
  LouvainResult<K> af = louvainSeq(x, init, {repeat, resolution, tolerance, passTolerance, toleranceDeclineFactor, 500, 1});
  printf("[%1.1e batch_size; %d batch_count; %09.3f ms; %04d iters.; %03d passes; %01.9f modularity] louvainSeqLast\n", 0.0, 0, af.time, af.iterations, af.passes, getModularity(x, af, M), tolerance, toleranceDeclineFactor);
  // Batch of additions only (dynamic).
  for (int batchSize=500, i=0; batchSize<=100000; batchSize*=i&1? 5:2, ++i) {
    for (int batchCount=1; batchCount<=5; ++batchCount) {
      auto y  = duplicate(x);
      auto fe = [&](auto u, auto v, auto w) { y.addEdge(u, v, w); };
      for (int i=0; i<batchSize; ++i)
        retry([&]() { return addRandomEdge(y, rnd, x.span(), V(1), fe); }, retries);
      y.correct();
      LouvainResult<K> bl = louvainSeq(y, &al.membership, {repeat, resolution, tolerance, passTolerance, toleranceDeclineFactor});
      printf("[%1.1e batch_size; %d batch_count; %09.3f ms; %04d iters.; %03d passes; %01.9f modularity] louvainSeqDynamicLast\n", double(batchSize), batchCount, bl.time, bl.iterations, bl.passes, getModularity(y, bl, M), tolerance, toleranceDeclineFactor);
      LouvainResult<K> bf = louvainSeq(y, &af.membership, {repeat, resolution, tolerance, passTolerance, toleranceDeclineFactor});
      printf("[%1.1e batch_size; %d batch_count; %09.3f ms; %04d iters.; %03d passes; %01.9f modularity] louvainSeqDynamicLast\n", double(batchSize), batchCount, bf.time, bf.iterations, bf.passes, getModularity(y, bf, M), tolerance, toleranceDeclineFactor);
    }
  }
  // Batch of deletions only (dynamic).
  for (int batchSize=500, i=0; batchSize<=100000; batchSize*=i&1? 5:2, ++i) {
    for (int batchCount=1; batchCount<=5; ++batchCount) {
      auto y  = duplicate(x);
      auto fe = [&](auto u, auto v) { y.removeEdge(u, v); };
      for (int i=0; i<batchSize; ++i)
        retry([&]() { return removeRandomEdge(y, rnd, x.span(), fe); }, retries);
      y.correct();
      LouvainResult<K> bl = louvainSeq(y, &al.membership, {repeat, resolution, tolerance, passTolerance, toleranceDeclineFactor});
      printf("[%1.1e batch_size; %d batch_count; %09.3f ms; %04d iters.; %03d passes; %01.9f modularity] louvainSeqDynamicLast\n", double(-batchSize), batchCount, bl.time, bl.iterations, bl.passes, getModularity(y, bl, M), tolerance, toleranceDeclineFactor);
      LouvainResult<K> bf = louvainSeq(y, &af.membership, {repeat, resolution, tolerance, passTolerance, toleranceDeclineFactor});
      printf("[%1.1e batch_size; %d batch_count; %09.3f ms; %04d iters.; %03d passes; %01.9f modularity] louvainSeqDynamicLast\n", double(-batchSize), batchCount, bf.time, bf.iterations, bf.passes, getModularity(y, bf, M), tolerance, toleranceDeclineFactor);
    }
  }
}


int main(int argc, char **argv) {
  using K = int;
  using V = float;
  char *file = argv[1];
  int repeat = argc>2? stoi(argv[2]) : 5;
  OutDiGraph<K, None, V> x; V w = 1;
  printf("Loading graph %s ...\n", file);
  readMtxW(x, file); println(x);
  auto y  = symmetricize(x); print(y); printf(" (symmetricize)\n");
  auto fl = [](auto u) { return true; };
  // selfLoopU(y, w, fl); print(y); printf(" (selfLoopAllVertices)\n");
  runLouvain(y, repeat);
  printf("\n");
  return 0;
}
