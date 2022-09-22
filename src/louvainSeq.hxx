#pragma once
#include <vector>
#include "_main.hxx"
#include "properties.hxx"
#include "duplicate.hxx"
#include "modularity.hxx"
#include "louvain.hxx"

using std::vector;








template <class G, class K, class V=float>
auto louvainSeq(const G& x, const vector<K>* q=nullptr, LouvainOptions<V> o={}) {
  V   R = o.resolution;
  V   D = o.passTolerance;
  int L = o.maxIterations, l = 0;
  int P = o.maxPasses, p = 0;
  K   S = x.span();
  V   M = edgeWeight(x)/2;
  vector<K> vcom(S), vcs, a(S);
  vector<V> vtot(S), ctot(S), vcout(S);
  float t = measureDurationMarked([&](auto mark) {
    V E  = o.tolerance;
    V Q0 = modularity(x, M, R);
    G y  = duplicate(x);
    fillValueU(vcom, K());
    fillValueU(vtot, V());
    fillValueU(ctot, V());
    mark([&]() {
      louvainVertexWeights(vtot, y);
      if (q) copyValues(*q, vcom);
      if (q) louvainCommunityWeights(ctot, y, vcom, vtot);
      else   louvainInitialize(vcom, ctot, y, vtot);
      copyValues(vcom, a);
      for (l=0, p=0; p<P;) {
        l += louvainMove(vcom, ctot, vcs, vcout, y, vtot, M, R, E, L);
        y  = louvainAggregate(y, vcom); ++p;
        louvainLookupCommunities(a, vcom);
        V Q = modularity(y, M, R);
        if (Q-Q0<=D) break;
        fillValueU(vcom, K());
        fillValueU(vtot, V());
        fillValueU(ctot, V());
        louvainVertexWeights(vtot, y);
        louvainInitialize(vcom, ctot, y, vtot);
        E /= o.tolerenceDeclineFactor;
        Q0 = Q;
      }
    });
  }, o.repeat);
  return LouvainResult<K>(a, l, p, t);
}
