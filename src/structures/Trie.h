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

// typedef tuple<int, int, double, double> res_t;
typedef tuple<int, int, shared_ptr<SA>, double> res_t;

using namespace std;

struct state_t {
  Node *cur;
  int seqIdx;
  int readIdx;
  int misMatch;
  double pen;
  const string &seq;
  const vector<double> &cErr;
  vector<res_t> &results;
  shared_ptr<SA> align;

  auto match() -> bool { return seq[seqIdx] == cur->getBase(); }

  auto phred() -> double { return cErr[seqIdx]; }

  auto leaf() -> bool { return cur->isLeaf(); }

  template <typename T> auto child(T &&t) -> Node * { return cur->getChild(t); }

  template <typename T, typename... Args> auto add(Args &&... args) -> void {
    align = make_alignment<T>(align, args...);
  }

  auto end() -> bool { return seqIdx == seq.size() - 1; }

  auto addResult() -> void {
    results.emplace_back(readIdx, dynamic_cast<Leaf *>(cur)->idx, align, 0.0);
  }
};

class Trie {
private:
  unique_ptr<Node> root;
  map<string, double> tMat;
  int corLength;
  array<double, 4> _penalties;
  double pen_max;

  mutex mut;

  auto getMatEle(const string &read, const int seqIdx) -> double;

  auto addSeq(const string &seq, const int idx) -> void;

  auto hammingSearch(state_t state) -> void;

  auto editSearch(state_t state, bool ind) -> void;

  auto indel(state_t state, bool left) -> void;

  auto extend(state_t state, bool left) -> void;
  // auto lev2Leaf(state_t state) -> void;

public:
  Trie() : root(new Node()) {}

  Trie(double gap, double ext, double max)
      : root(new Node()), _penalties{{gap, ext, gap, ext}}, pen_max(max) {}

  Trie(double gap0, double ext0, double gap1, double ext1, double max)
      : root(new Node()), _penalties{{gap0, ext0, gap1, ext1}}, pen_max(max) {}

  ~Trie() {}

  auto penalties() -> array<double, 4> & { return _penalties; }

  auto max() -> double { return pen_max; }

  auto bounded() -> bool {
    return pen_max != 0 &&
           (_penalties[0] || _penalties[1] || _penalties[2] || _penalties[3]);
  }

  auto fromLibrary(const vector<string> &library) -> void;

  auto setTMat(const Rcpp::StringVector tMatSeq,
               const Rcpp::NumericVector tMatProb) -> bool;

  auto count(vector<res_t> &results, vector<double> &countTable) -> void;

  auto count(vector<res_t> &results, vector<double> &countTable,
             std::ostream &out) -> void;

  template <typename... Args> auto hamming(Args &&... args) -> void {
    for (auto i = 0; i < 4; i++) {
      auto state = state_t{root->getChild(i), args..., shared_ptr<SA>()};
      hammingSearch(state);
    }
  }

  template <typename... Args> auto edit(Args &&... args) -> void {
    for (auto i = 0; i < 4; i++) {
      auto state = state_t{root->getChild(i), args..., shared_ptr<SA>()};
      editSearch(state, false);
    }
  }
};

#endif
