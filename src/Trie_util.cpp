// Implementation of utility methods for the Trie class.
#include "structures/Trie.h"

// Method for adding a sequence to the Trie data structure.
auto Trie::addSeq(const string &seq, const int idx) -> void {
  if (seq.empty()) {
    return;
  }

  auto cur = root;
  for (auto i = seq.begin(); i < seq.end(); ++i) {
    auto child = pool[cur].getChild(*i);

    if (child != -1) {
      cur = child;
    } else {
      cur = pool.add_child(cur, *i, i + 1 == seq.end() ? idx : -1);
    }
  }
}

// Method for creating a Trie data structure from a library.
auto Trie::fromLibrary(const vector<string> &library) -> void {
  for (auto x = library.begin(); x < library.end(); x++) {
    addSeq(*x, x - library.begin());
  }

  Rcpp::Rcout << "Trie(" << pool.count << ") created" << std::endl;
}

// Method for creating count table for barcodes in the library
// from the alignment results.
auto Trie::count(vector<res_t> &results, vector<double> &countTable) -> void {
  // sort by read then probability
  sort(results.begin(), results.end(), [](const res_t &x, const res_t &y) {
    return get<0>(x) < get<0>(y) ||
           (get<0>(x) == get<0>(y) && get<2>(x)->value() < get<2>(y)->value());
  });

  // lock + count maximum per read
  lock_guard<mutex> lock(mut);

  auto last = results.begin();
  for (auto x = results.begin(); x < results.end(); x++) {
    if (x == results.end() - 1 || get<0>(*(x + 1)) != get<0>(*last)) {
      countTable[get<1>(*x)] += 1;
      last = x + 1;
    }
  }
}


void Trie::add_results(std::vector<res_t>& res) {
  lock_guard<mutex> lock(mut);
  std::move(res.begin(), res.end(), back_inserter(results));
}


auto Trie::count(vector<res_t> &results, vector<double> &countTable,
                 std::ostream &out) -> void {
  // sort by read then probability
  sort(results.begin(), results.end(), [](const res_t &x, const res_t &y) {
    return get<0>(x) < get<0>(y) ||
           (get<0>(x) == get<0>(y) && get<2>(x)->value() < get<2>(y)->value());
  });

  // lock + count maximum per read
  lock_guard<mutex> lock(mut);

  auto last = results.begin();
  for (auto x = results.begin(); x < results.end(); x++) {
    if (x == results.end() - 1 || get<0>(*(x + 1)) != get<0>(*last)) {
      countTable[get<1>(*x)] += 1;
      out << get<0>(*x) << "," << get<1>(*x) << "\n";
      last = x + 1;
    }
  }
}

// Method for set the transformation matrix used in the probability model.
// This is advanced usage, notification for caution usage will be printed
// to the console for users.
auto Trie::setTMat(const Rcpp::StringVector tMatSeq,
                   const Rcpp::NumericVector tMatProb) -> bool {
  // ask user to provide the tMat; if not provided will use default 1/3
  int tMatSize = tMatSeq.size();
  if (tMatSize != tMatProb.size()) {
    Rcpp::Rcout << "Error: tMatSeq and tMatProb must be same length" << endl;
    return false;
  }

  if (tMatSize > 0) {
    corLength = ((string)tMatSeq[0]).length();
    string MatSeq;
    double MatProb;

    for (int i = 0; i < tMatSize; i++) {
      MatSeq = tMatSeq[i];
      if (MatSeq.length() != (std::size_t)corLength) {
        Rcpp::Rcout << "Error: All tMatSeq entries must be the same length"
                    << endl;
        return false;
      }

      MatProb = tMatProb[i];
      tMat.insert(pair<string, double>(MatSeq, MatProb));
    }
  }

  return true;
}

// Method for obtaining transfermation matrix elements.
auto Trie::getMatEle(const string &read, const int seqIdx) -> double {
  if (seqIdx - corLength + 1 > 0) {
    auto sub = read.substr(seqIdx - corLength + 1, corLength);
    auto val = tMat.find(sub);

    if (val != tMat.end()) {
      return val->second;
    }
  }

  return 0.333333;
}
