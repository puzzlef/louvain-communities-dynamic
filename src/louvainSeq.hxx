#pragma once
#include <utility>
#include <vector>
#include "_main.hxx"
#include "properties.hxx"
#include "duplicate.hxx"
#include "modularity.hxx"
#include "louvain.hxx"

using std::pair;
using std::vector;




// LOUVAIN-SEQ
// -----------

template <class G, class K, class V, class FA>
auto louvainSeq(const G& x, const vector<K>* q, const LouvainOptions<V>& o, FA fa) {
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
        if (p==0) l += louvainMove(vcom, ctot, vcs, vcout, y, vtot, M, R, E, L, fa);
        else      l += louvainMove(vcom, ctot, vcs, vcout, y, vtot, M, R, E, L);
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




// LOUVAIN-SEQ-STATIC
// ------------------

template <class G, class K, class V=float>
inline auto louvainSeqStatic(const G& x, const vector<K>* q=nullptr, const LouvainOptions<V>& o={}) {
  auto fa = [](auto u) { return true; };
  return louvainSeq(x, q, o, fa);
}




// LOUVAIN-SEQ-DYNAMIC
// -------------------

template <class G, class K, class V=float>
inline auto louvainSeqDynamic(const G& x, const vector<pair<K, K>>& deletions, const vector<pair<K, K>>& insertions, const vector<K>* q, const LouvainOptions<V>& o={}) {
  K S = x.span();
  V R = o.resolution;
  V M = edgeWeight(x)/2;
  const vector<K>& vcom = *q;
  vector<V> vtot(S), ctot(S);
  louvainVertexWeights(vtot, x);
  louvainCommunityWeights(ctot, x, vcom, vtot);
  auto vaff = louvainAffectedVertices(x, deletions, insertions, vcom, vtot, ctot, M, R);
  auto fa   = [&](auto u) { return vaff[u]==true; };
  return louvainSeq(x, q, o, fa);
}
