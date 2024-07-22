Comparing *static* vs *dynamic* approaches of the [Louvain algorithm] for
[community detection].

[Louvain] is an algorithm for **detecting communities in graphs**. *Community*
*detection* helps us understand the *natural divisions in a network* in an
**unsupervised manner**. It is used in **e-commerce** for *customer*
*segmentation* and *advertising*, in **communication networks** for *multicast*
*routing* and setting up of *mobile networks*, and in **healthcare** for
*epidemic causality*, setting up *health programmes*, and *fraud detection* is
hospitals. *Community detection* is an **NP-hard** problem, but heuristics exist
to solve it (such as this). **Louvain algorithm** is an **agglomerative-hierarchical**
community detection method that **greedily optimizes** for **modularity**
(**iteratively**).

**Modularity** is a score that measures *relative density of edges inside* vs
*outside* communities. Its value lies between `−0.5` (*non-modular clustering*)
and `1.0` (*fully modular clustering*). Optimizing for modularity *theoretically*
results in the best possible grouping of nodes in a graph.

Given an *undirected weighted graph*, all vertices are first considered to be
*their own communities*. In the **first phase**, each vertex greedily decides to
move to the community of one of its neighbors which gives greatest increase in
modularity. If moving to no neighbor's community leads to an increase in
modularity, the vertex chooses to stay with its own community. This is done
sequentially for all the vertices. If the total change in modularity is more
than a certain threshold, this phase is repeated. Once this **local-moving**
**phase** is complete, all vertices have formed their first hierarchy of
communities. The **next phase** is called the **aggregation phase**, where all
the *vertices belonging to a community* are *collapsed* into a single
**super-vertex**, such that edges between communities are represented as edges
between respective super-vertices (edge weights are combined), and edges within
each community are represented as self-loops in respective super-vertices
(again, edge weights are combined). Together, the local-moving and the
aggregation phases constitute a **pass**. This super-vertex graph is then used
as input for the next pass. This process continues until the increase in
modularity is below a certain threshold. As a result from each pass, we have a
*hierarchy of community memberships* for each vertex as a **dendrogram**. We
generally consider the *top-level hierarchy* as the *final result* of community
detection process.

*Louvain* algorithm is a hierarchical algorithm, and thus has two different
tolerance parameters: `tolerance` and `passTolerance`. **tolerance** defines the
minimum amount of increase in modularity expected, until the local-moving phase
of the algorithm is considered to have converged. We compare the increase in
modularity in each iteration of the local-moving phase to see if it is below
`tolerance`. **passTolerance** defines the minimum amount of increase in
modularity expected, until the entire algorithm is considered to have converged.
We compare the increase in modularity across all iterations of the local-moving
phase in the current pass to see if it is below `passTolerance`. `passTolerance`
is normally set to `0` (we want to maximize our modularity gain), but the same
thing does not apply for `tolerance`. Adjusting values of `tolerance` between
each pass have been observed to impact the runtime of the algorithm, without
significantly affecting the modularity of obtained communities. In this
experiment, we compare the performance of *three different types* of **dynamic**
**Louvain** with respect to the *static* version.

**Naive dynamic**:
- We start with previous community membership of each vertex (instead of each vertex its own community).

**Delta screening**:
- All edge batches are undirected, and sorted by source vertex-id.
- For edge additions across communities with source vertex `i` and highest modularity changing edge vertex `j*`,
  `i`'s neighbors and `j*`'s community is marked as affected.
- For edge deletions within the same community `i` and `j`,
  `i`'s neighbors and `j`'s community is marked as affected.

**Frontier-based**:
- All source and destination vertices are marked as affected for insertions and deletions.
- For edge additions across communities with source vertex `i` and destination vertex `j`,
  `i` is marked as affected.
- For edge deletions within the same community `i` and `j`,
  `i` is marked as affected.
- Vertices whose communities change in local-moving phase have their neighbors marked as affected.

The input data used for the experiments is available from the
[SuiteSparse Matrix Collection]. These experiments are done with guidance
from [Prof. Kishore Kothapalli] and [Prof. Dip Sankar Banerjee].

<br>


### Comparing various Naive-dynamic approaches

In this experiment ([approaches-naive]), we compare the performance of *two*
*different types* of **naive dynamic Louvain** with respect to the *static*
version. The **last** approach (`louvainSeqDynamicLast`) considers the community
membership of each vertex *after* the Louvain algorithm has *converged*
(community membership from the "last" pass) and then performs the Louvain
algorithm upon the new (updated) graph. This is *similar* to naive dynamic
approaches with other algorithms. On the other hand, the **first** approach
(`louvainSeqDynamicFirst`) considers the community membership of each vertex
right after the *first pass* of the Louvain algorithm (this is the first
community membership hierarchy) and then performs the Louvain algorithm upon the
updated graph. With this approach, we allow the affected vertices to choose
their community membership from the first pass itself, which which to my
intuition would lead to better communities.

First, we compute the community membership of each vertex using the static
Louvain algorithm (`louvainSeqLast`). We also run the static Louvain algorithm
for only one pass (`louvainSeqFirst`). We then generate *batches* of *insertions*
*(+)* and *deletions (-)* of edges of sizes 500, 1000, 5000, ... 100000. For each
batch size, we generate *five* different batches for the purpose of *averaging*.
Each batch of edges (insertion / deletion) is generated randomly such that the
selection of each vertex (as endpoint) is *equally probable*. We choose the
Louvain *parameters* as `resolution = 1.0`, `tolerance = 1e-2` (for local-moving
phase) with *tolerance* decreasing after every pass by a factor of
`toleranceDeclineFactor = 10`, and a `passTolerance = 0.0` (when passes stop).
In addition we limit the maximum number of iterations in a single local-moving
phase with `maxIterations = 500`, and limit the maximum number of passes with
`maxPasses = 500`. We run the Louvain algorithm until convergence (or until the
maximum limits are exceeded), and measure the **time** **taken** for the
*computation* (performed 5 times for averaging), the **modularity score**, the
**total number of iterations** (in the *local-moving* *phase*), and the number
of **passes**. This is repeated for *seventeen* different graphs.

From the results, we make make the following observations. The performance of
dynamic approaches upon a batch of deletions appears to *increase* with *increasing*
batch size*. This makes sense since, as the graph keeps getting smaller, the
computation would complete *sooner*. Next, the `first` naive dynamic approach is
found to be *significantly slower* (~0.3x speedup) than the `last` approach.
However, the `first` approach is *still faster* than the static approach upto a
batch size of `50000`. On the other hand, the `last` approach is *faster* than the
static approach for all batch sizes. A similar behavior is observed with the
total number of iterations. The `first` approach seems to have a *slightly higher*
modularity with respect to the `last` approach. Since the modularity between the
two dynamic approaches are almost the same, the **last** approach is clearly the
**best choice**.

[approaches-naive]: https://github.com/puzzlef/louvain-communities-dynamic/tree/approaches-naive

<br>


### Comparision with Static approach

First ([compare-static], [main]), we compute the community membership of each vertex
using the static Louvain algorithm. We then generate *batches* of *insertions*
*(+)* and *deletions (-)* of edges of sizes 500, 1000, 5000, ... 100000. For
each batch size, we generate *five* different batches for the purpose of
*averaging*. Each batch of edges (insertion / deletion) is generated randomly
such that the selection of each vertex (as endpoint) is *equally probable*. We
choose the Louvain *parameters* as `resolution = 1.0`, `tolerance = 1e-2` (for
local-moving phase) with *tolerance* decreasing after every pass by a factor of
`toleranceDeclineFactor = 10`, and a `passTolerance = 0.0` (when passes stop).
In addition we limit the maximum number of iterations in a single local-moving
phase with `maxIterations = 500`, and limit the maximum number of passes with
`maxPasses = 500`. We run the Louvain algorithm until convergence (or until the
maximum limits are exceeded), and measure the **time taken** for the
*computation* (performed 5 times for averaging), the **modularity score**, the
**total number of iterations** (in the *local-moving* *phase*), and the number
of **passes**. This is repeated for *seventeen* different graphs.

From the results, we make make the following observations. The frontier-based
dynamic approach converges the fastest, which obtaining communities with only
slightly lower modularity than other approaches. We also observe that
delta-screening based dynamic Louvain algorithm has the same performance as that
of the naive dynamic approach. Therefore, **frontier-based dynamic Louvain**
would be the **best choice**. All outputs are saved in a [gist]. Some [charts]
are also included below, generated from [sheets].

[![](https://i.imgur.com/MVFks1a.png)][sheetp]
[![](https://i.imgur.com/du5qe69.png)][sheetp]
[![](https://i.imgur.com/AulWzoj.png)][sheetp]
[![](https://i.imgur.com/bxpYvsF.png)][sheetp]

[compare-static]: https://github.com/puzzlef/louvain-communities-dynamic/tree/compare-static
[main]: https://github.com/puzzlef/louvain-communities-dynamic

<br>

<br>
<br>


## References

- [Fast unfolding of communities in large networks; Vincent D. Blondel et al. (2008)](https://arxiv.org/abs/0803.0476)
- [Community Detection on the GPU; Md. Naim et al. (2017)](https://arxiv.org/abs/1305.2006)
- [Scalable Static and Dynamic Community Detection Using Grappolo; Mahantesh Halappanavar et al. (2017)](https://ieeexplore.ieee.org/document/8091047)
- [From Louvain to Leiden: guaranteeing well-connected communities; V.A. Traag et al. (2019)](https://www.nature.com/articles/s41598-019-41695-z)
- [Community detection in networks: A multidisciplinary review; Muhammad Aqib Javed et al. (2018)](https://www.sciencedirect.com/science/article/abs/pii/S1084804518300560)
- [Adaptive parallel Louvain community detection on a multicore platform; Mahmood Fazlali et al. (2017)](https://www.sciencedirect.com/science/article/abs/pii/S014193311630240X)
- [A signal processing perspective to community detection in dynamic networks; Selin Aviyente (2021)](https://www.sciencedirect.com/science/article/abs/pii/S1051200421002311)
- [Incremental Community Detection in Distributed Dynamic Graph; Tariq Abughofa et al. (2021)](https://ieeexplore.ieee.org/abstract/document/9564329)
- [Community Detection in Dynamic Attributed Graphs; Gonzalo A. Bello et al. (2016)](https://link.springer.com/chapter/10.1007/978-3-319-49586-6_22)
- [Dynamic Community Detection Using Nonnegative Matrix Factorization; Feng Gao et al. (2017)](https://ieeexplore.ieee.org/abstract/document/8327706)
- [Improving fuzzy C-mean-based community detection in social networks using dynamic parallelism; Mahmoud Al-Ayyoub et al. (2019)](https://www.sciencedirect.com/science/article/abs/pii/S0045790617315641)
- [Distributed Louvain Algorithm for Graph Community Detection; Sayan Ghosh et al. (2018)](https://ieeexplore.ieee.org/abstract/document/8425242)
- [Evaluating community detection algorithms for progressively evolving graphs; Remy Cazabet et al. (2020)](https://academic.oup.com/comnet/article-abstract/8/6/cnaa027/6161494?login=false)
- [Scalable distributed Louvain algorithm for community detection in large graphs; Naw Safrin Sattar et al. (2022)](https://link.springer.com/article/10.1007/s11227-021-04224-2)
- [Community Detection in Evolving Networks; Tejas Puranik et al. (2017)](https://ieeexplore.ieee.org/abstract/document/9069125)
- [Incremental local community identification in dynamic social networks; Mansoureh Takaffoli et al. (2013)](https://ieeexplore.ieee.org/abstract/document/6785692)
- [Dynamic community detection based on the Matthew effect; Zejun Sun et al. (2022)](https://www.sciencedirect.com/science/article/abs/pii/S0378437122002564)
- [Community detection in social networks; Punam Bedi et al. (2016)](https://wires.onlinelibrary.wiley.com/doi/abs/10.1002/widm.1178)
- [Long range community detection; Thomas Aynaud et al. (2010)](https://inria.hal.science/inria-00531750/)
- [Dynamic Community Detection Algorithm Based on Incremental Identification; Xiaoming Li et al. (2015)](https://ieeexplore.ieee.org/abstract/document/7395763)
- [An Improved Louvain Algorithm for Community Detection; Jicun Zhang et al. (2021)](https://onlinelibrary.wiley.com/doi/10.1155/2021/1485592)
- [DyPerm: Maximizing Permanence for Dynamic Community Detection; Prerna Agarwal et al. (2018)](https://link.springer.com/chapter/10.1007/978-3-319-93034-3_35)
- [A fast and efficient incremental approach toward dynamic community detection; Neda Zarayeneh et al. (2019)](https://dl.acm.org/doi/abs/10.1145/3341161.3342877)
- [Scalable Community Detection with the Louvain Algorithm; Xinyu Que et al. (2015)](https://ieeexplore.ieee.org/abstract/document/7161493)
- [DynaMo: Dynamic Community Detection by Incrementally Maximizing Modularity; Di Zhuang et al. (2019)](https://ieeexplore.ieee.org/abstract/document/8890861)
- [C-Blondel: An Efficient Louvain-Based Dynamic Community Detection Algorithm; Mahsa Seifikar et al. (2020)](https://ieeexplore.ieee.org/abstract/document/8982190)
- [Scalable static and dynamic community detection using Grappolo; Mahantesh Halappanavar et al. (2017)](https://ieeexplore.ieee.org/abstract/document/8091047)
- [Dynamic community detection in evolving networks using locality modularity optimization; Mário Cordeiro et al. (2016)](https://link.springer.com/article/10.1007/s13278-016-0325-1)
- [A Dynamic Modularity Based Community Detection Algorithm for Large-scale Networks: DSLM; Riza Aktunc et al. (2015)](https://dl.acm.org/doi/abs/10.1145/2808797.2808822)
- [An Analysis of the Dynamic Community Detection Algorithms in Complex Networks; Dhananjay Kumar Singh et al. (2020)](https://ieeexplore.ieee.org/abstract/document/9067224)
- [Dynamic Community Detection Decouples Multiple Time Scale Behavior of Complex Chemical Systems; Neda Zarayeneh et al. (2022)](https://pubs.acs.org/doi/full/10.1021/acs.jctc.2c00454)
- [Targeted revision: A learning-based approach for incremental community detection in dynamic networks; Jiaxing Shang et al. (2016)](https://www.sciencedirect.com/science/article/abs/pii/S0378437115008080)
- [DynComm R Package -- Dynamic Community Detection for Evolving Networks; Rui Portocarrero Sarmento et al. (2019)](https://arxiv.org/abs/1905.01498)
- [CS224W: Machine Learning with Graphs | Louvain Algorithm; Jure Leskovec (2021)](https://www.youtube.com/watch?v=0zuiLBOIcsw)
- [The University of Florida Sparse Matrix Collection; Timothy A. Davis et al. (2011)](https://doi.org/10.1145/2049662.2049663)

<br>
<br>


[![](https://i.imgur.com/UGB0g2L.jpg)](https://www.youtube.com/watch?v=pIF3wOet-zw)<br>
[![ORG](https://img.shields.io/badge/org-puzzlef-green?logo=Org)](https://puzzlef.github.io)
[![DOI](https://zenodo.org/badge/538336155.svg)](https://zenodo.org/badge/latestdoi/538336155)


[Prof. Dip Sankar Banerjee]: https://sites.google.com/site/dipsankarban/
[Prof. Kishore Kothapalli]: https://faculty.iiit.ac.in/~kkishore/
[SuiteSparse Matrix Collection]: https://sparse.tamu.edu
[Louvain algorithm]: https://en.wikipedia.org/wiki/Louvain_method
[community detection]: https://en.wikipedia.org/wiki/Community_search
[Louvain]: https://en.wikipedia.org/wiki/Louvain_method
[gist]: https://gist.github.com/wolfram77/de2c1e1c8f6efb7f4053a122b688c7a7
[charts]: https://imgur.com/a/xcoVmDw
[sheets]: https://docs.google.com/spreadsheets/d/1U8cdi0Y9i-Pl0SkUdSsuM87Khw0Jn6SCM7u6qlVp--Y/edit?usp=sharing
[sheetp]: https://docs.google.com/spreadsheets/d/e/2PACX-1vQjZGdwSy2Cd5mojCQgvvl0P5eaZBjJwhIBqcX-RN5MfeFVl9V2o9G4XCjR2_qaH_mpeBSm7eMIhSeQ/pubhtml
