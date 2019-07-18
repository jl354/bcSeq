#ifndef __TRIE_H__
#define __TRIE_H__

#include <Rcpp.h>
#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "../alignments/alignment.hpp"
#include "Node.h"
#include "NodePool.h"

// typedef tuple<int, int, double, double> res_t;
typedef tuple<int, int, shared_ptr<SA>, double> res_t;

using namespace std;

struct state_t {
  const Pool& pool;
  int cur;
  int seqIdx;
  int readIdx;
  int misMatch;
  double pen;
  const string &seq;
  const vector<double> &cErr;
  vector<res_t> &results;
  shared_ptr<SA> align;

  bool match() { return seq[seqIdx] == pool[cur].baseType; }

  double phred() { return cErr[seqIdx]; }

  bool leaf() { return pool[cur].isLeaf(); }

  int child(int i) { return pool[cur].children[i]; }

  template <typename T, typename... Args> void add(Args &&... args) {
    align = make_alignment<T>(align, args...);
  }

  bool end() { return static_cast<std::size_t>(seqIdx) == seq.size() - 1; }

  void addResult() {
    results.emplace_back(readIdx, pool[cur].leaf, align, 0.0);
  }
};

struct Trie {
  Pool pool;
  int root;
  map<string, double> tMat;
  int corLength;
  array<double, 4> _penalties;
  double pen_max;
  mutex mut;
  std::vector<res_t> results;

  double getMatEle(const string &read, const int seqIdx);

  void addSeq(const string &seq, const int idx);

  void hammingSearch(state_t state);

  void editSearch(state_t state, bool ind);

  void indel(state_t state, bool left);

  void extend(state_t state, bool left);

Trie() : pool(), root(pool.createNode('N', 0)), results() {}

  Trie(double gap, double ext, double max)
  : pool(), root(pool.createNode('N', 0)), _penalties{{gap, ext, gap, ext}}, pen_max(max), results() {}

  Trie(double gap0, double ext0, double gap1, double ext1, double max)
  : pool(), root(pool.createNode('N', 0)), _penalties{{gap0, ext0, gap1, ext1}}, pen_max(max), results() {}

  ~Trie() {}

  void add_results(std::vector<res_t>& results);
  
  array<double, 4> &penalties() { return _penalties; }

  double max() { return pen_max; }

  bool bounded() {
    return pen_max != 0 &&
           (_penalties[0] || _penalties[1] || _penalties[2] || _penalties[3]);
  }

  void fromLibrary(const vector<string> &library);

  bool setTMat(const Rcpp::StringVector tMatSeq,
               const Rcpp::NumericVector tMatProb);

  void count(vector<res_t> &results, vector<double> &countTable);

  void count(vector<res_t> &results, vector<double> &countTable,
             std::ostream &out);

  template <typename... Args> void hamming(Args &&... args) {
    for (auto i = 0; i < 4; i++) {
      if(pool[root].children[i] != -1){
        auto state = state_t{pool, pool[root].children[i], args..., shared_ptr<SA>()};
        hammingSearch(state);
      }
    }
  }

  template <typename... Args> void edit(Args &&... args) {    
    for (auto i = 0; i < 4; i++) {
      if(pool[root].children[i] != -1){
        auto state = state_t{pool, pool[root].children[i], args..., shared_ptr<SA>()};
        editSearch(state, false);
      }
    }
  }
};

#endif
