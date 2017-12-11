// The main function for conducting the alignment and assigning mapping
// probability
#include "helper/count.h"
#include "helper/fastqCheck.h"
#include "structures/Trie.h"
#include <cmath>
#include <thread>

using namespace Rcpp;
using namespace std;

// sort by barcode
bool by_bc(const res_t &x, const res_t &y) { return get<1>(x) < get<1>(y); }

// sort by pr
bool by_pr(const res_t &x, const res_t &y) { return get<3>(x) < get<3>(y); }

auto default_transform(double max, double pr, double penalty) -> double {
  return pr * (1.0 - log(2.0) + log1p(max / (max + penalty)));
}

/*
*  Takes maximum over a given read's alignments
*  to a given barcode, replaces the current set of results
*  with set of maxima
*/

void clean(Trie &trie, vector<res_t>::iterator start, vector<res_t> &results) {
  auto keep = vector<res_t>();

  auto m = trie.max();
  auto pen = trie.penalties();

  // fitness function
  for (auto x = start; x < results.end(); x++) {
    get<3>(*x) =
        default_transform(m, get<2>(*x)->value(), get<2>(*x)->penalty(pen));
  }

  // sort by barcode
  sort(start, results.end(), by_bc);

  // take max over each barcode for this read
  auto last = start;
  for (auto i = start; i < results.end(); i++) {
    if (i + 1 == results.end() || get<1>(*(i + 1)) != get<1>(*last)) {
      keep.push_back(*max_element(last, i + 1, by_pr));
      last = i + 1;
    }
  }

  auto end = move(keep.begin(), keep.end(), start);
  results.erase(end, results.end());
}

/*
*  find maxima per mapped barcode for a given read
*/
auto extract(vector<res_t>::iterator start, vector<res_t>::iterator end,
             vector<res_t> &keep, double m, array<double, 4> &pen,
             Function tForm) -> void {
  NumericVector val, pens, out;
  for (auto x = start; x < end; x++) {
    val.push_back(get<2>(*x)->value());
    pens.push_back(get<2>(*x)->penalty(pen));
  }
  out = tForm(m, val, pens);

  // take max over barcodes for this read
  auto mx = max_element(out.begin(), out.end());
  keep.push_back(*(start + (mx - out.begin())));
  get<3>(keep.back()) = *mx;
}

void clean(Trie &trie, vector<res_t> &results, Function tForm) {
  auto keep = vector<res_t>();

  auto m = trie.max();
  auto pen = trie.penalties();

  sort(results.begin(), results.end(),
       [](const res_t &x, const res_t &y) { return get<0>(x) < get<0>(y); });

  // loop over all results
  // whenever we reach a read boundary, find barcode boundaries
  // and take maximum over the barcode
  auto last = results.begin();
  for (auto x = results.begin(); x < results.end(); x++) {
    if (x + 1 == results.end() || get<0>(*(x + 1)) != get<0>(*last)) {
      sort(last, x + 1, by_bc);

      // find barcode boundaries
      auto last2 = last;
      for (auto y = last2; y < x + 1; y++) {
        if (y + 1 == x + 1 || get<1>(*(y + 1)) != get<1>(*last2)) {
          extract(last2, y + 1, keep, m, pen, tForm);
          last2 = y + 1;
        }
      }

      last = x + 1;
    }
  }

  results = keep;
}

void user_alignment(Trie &trie, vector<string> &sequence,
                    vector<string> &phredScore, int misMatch,
                    vector<double> &countTable, int L, int R,
                    vector<res_t> &all_res, mutex &mut, bool count_only) {
  auto cError = vector<double>();
  auto results = std::vector<res_t>();

  for (int i = L; i < R; i++) {
    // auto s = results.size(); //warning from clang variable not used
    phred2err(cError, phredScore[i]);
    trie.edit(0, i, misMatch, 0.0, sequence[i], cError, results);
  }

  // must collect all results so we can call R code on it
  lock_guard<mutex> l(mut);
  move(results.begin(), results.end(), back_inserter(all_res));
}

void alignment(Trie &trie, vector<string> &sequence, vector<string> &phredScore,
               int misMatch, vector<double> &countTable, int L, int R,
               vector<res_t> &all_res, mutex &mut, bool count_only,
               ostream &table_stream, bool detail_info) {
  auto cError = vector<double>();
  auto results = std::vector<res_t>();

  for (int i = L; i < R; i++) {
    auto s = results.size();
    phred2err(cError, phredScore[i]);
    trie.edit(0, i, misMatch, 0.0, sequence[i], cError, results);
    if (trie.bounded()) {
      clean(trie, results.begin() + s, results);
    }
  }

  if (count_only)
    trie.count(results, countTable);
  if (detail_info)
    trie.count(results, countTable, table_stream);

  if (!count_only) {
    lock_guard<mutex> l(mut);
    move(results.begin(), results.end(), back_inserter(all_res));
  }
}

void alignmentH(Trie &trie, vector<string> &sequence,
                vector<string> &phredScore, int misMatch,
                vector<double> &countTable, int L, int R,
                vector<res_t> &all_res, mutex &mut, bool count_only,
                ostream &table_stream, bool detail_info) {
  auto cError = vector<double>();
  auto results = std::vector<res_t>();

  for (int i = L; i < R; i++) {
    phred2err(cError, phredScore[i]);
    trie.hamming(0, i, misMatch, 0.0, sequence[i], cError, results);
  }

  if (count_only)
    trie.count(results, countTable);
  if (detail_info)
    trie.count(results, countTable, table_stream);

  if (!count_only) {
    for (auto &x : results) {
      get<3>(x) = get<2>(x)->value();
    }

    lock_guard<mutex> l(mut);
    move(results.begin(), results.end(), back_inserter(all_res));
  }
}

//[[Rcpp::export]]
SEXP CRISPR_matching(String sampleFile, String libFile, String outFile,
                     int misMatch, Rcpp::StringVector tMatSeq,
                     Rcpp::NumericVector tMatProb, int numThread, bool hamming,
                     bool count_only, double gap_left, double ext_left,
                     double gap_right, double ext_right, double pen_max,
                     bool detail_info) {
  // read in fastq file and store sequence and phred score as vector string.
  auto sequence = vector<string>(), sequence_ids = vector<string>(),
       phredScore = vector<string>(), library = vector<string>(),
       library_ids = vector<string>();

  auto threads = std::vector<thread>();
  auto match = hamming ? alignmentH : alignment;

  Trie trie(gap_left, ext_left, gap_right, ext_right, pen_max);

  // read samples, library, and setup trie's tMat
  if (!(readSamples(sampleFile, sequence, sequence_ids, phredScore) &&
        trie.setTMat(tMatSeq, tMatProb) &&
        readLibrary(library, library_ids, libFile))) {
    return R_NilValue;
  }

  auto countTable = vector<double>(library.size());

  // create trie with library elements
  trie.fromLibrary(library);

  mutex mut;
  vector<res_t> all_res;

  std::string table_out = outFile;
  std::ofstream table_stream(table_out + ".txt");

  // run alignment
  try {

    int nr_item = sequence.size();
    int nr_item_per_thread = ceil(nr_item * 1.0 / numThread);

    Rcpp::Rcout << "Running"
                << (hamming ? " hamming search" : " levenshtein search")
                << " with " << nr_item_per_thread << " sequences per thread in "
                << numThread << " threads" << endl;

    for (int i = 1; i < numThread; ++i) {
      int R = min((i + 1) * nr_item_per_thread, nr_item);
      threads.emplace_back(match, std::ref(trie), std::ref(sequence),
                           std::ref(phredScore), misMatch, std::ref(countTable),
                           i * nr_item_per_thread, R, std::ref(all_res),
                           std::ref(mut), count_only, std::ref(table_stream),
                           detail_info);
    }

    int R = min(nr_item_per_thread, nr_item);
    match(trie, sequence, phredScore, misMatch, countTable, 0, R, all_res, mut,
          count_only, table_stream, detail_info);

    for (auto &t : threads) {
      t.join();
    }

  } catch (Rcpp::exception &e) {
    Rcpp::Rcout << e.what() << endl;
    return R_NilValue;
  }

  if (detail_info)
    trie.count(all_res, countTable, table_stream);

  Rcpp::Rcout << "Compiling results\n";

  count2CSV(countTable, library, outFile);

  if (!count_only) {
    Rcpp::Rcout << "Generating dataframe\n";
    IntegerVector readIdx;
    IntegerVector bcIdx;
    auto prob = vector<double>();

    for (auto &x : all_res) {
      readIdx.push_back(std::get<0>(x));
      bcIdx.push_back(std::get<1>(x));
      prob.push_back(std::get<3>(x));
    }

    return List::create(List::create(Named("reads") = sequence_ids,
                                     Named("barcodes") = library_ids),
                        List::create(Named("i") = readIdx, Named("j") = bcIdx,
                                     Named("x") = prob,
                                     Named("index1") = false));
  }

  return List::create(Named("reads") = sequence_ids,
                      Named("barcodes") = library_ids);
}

//[[Rcpp::export]]
SEXP CRISPR_user_matching(String sampleFile, String libFile, String outFile,
                          int misMatch, Rcpp::StringVector tMatSeq,
                          Rcpp::NumericVector tMatProb, int numThread,
                          bool count_only, double gap_left, double ext_left,
                          double gap_right, double ext_right, double pen_max,
                          Function tForm) {
  // read in fastq file and store sequence and phred score as vector string.
  auto sequence = vector<string>(), sequence_ids = vector<string>(),
       phredScore = vector<string>(), library = vector<string>(),
       library_ids = vector<string>();

  auto threads = std::vector<thread>();

  Trie trie(gap_left, ext_left, gap_right, ext_right, pen_max);

  // read samples, library, and setup trie's tMat
  if (!(readSamples(sampleFile, sequence, sequence_ids, phredScore) &&
        trie.setTMat(tMatSeq, tMatProb) &&
        readLibrary(library, library_ids, libFile))) {
    return R_NilValue;
  }

  auto countTable = vector<double>(library.size());

  // create trie with library elements
  trie.fromLibrary(library);

  mutex mut;
  vector<res_t> all_res;

  // run alignment
  try {

    int nr_item = sequence.size();
    int nr_item_per_thread = ceil(nr_item * 1.0 / numThread);

    Rcpp::Rcout << "Running"
                << " levenshtein search with " << nr_item_per_thread
                << " sequences per thread in " << numThread << " threads"
                << endl;

    for (int i = 1; i < numThread; ++i) {
      int R = min((i + 1) * nr_item_per_thread, nr_item);
      threads.emplace_back(user_alignment, std::ref(trie), std::ref(sequence),
                           std::ref(phredScore), misMatch, std::ref(countTable),
                           i * nr_item_per_thread, R, std::ref(all_res),
                           std::ref(mut), count_only);
    }

    int R = min(nr_item_per_thread, nr_item);
    user_alignment(trie, sequence, phredScore, misMatch, countTable, 0, R,
                   all_res, mut, count_only);

    for (auto &t : threads) {
      t.join();
    }

  } catch (Rcpp::exception &e) {
    Rcpp::Rcout << e.what() << endl;
    return R_NilValue;
  }

  clean(trie, all_res, tForm);

  trie.count(all_res, countTable);

  Rcpp::Rcout << "Compiling results\n";

  count2CSV(countTable, library, outFile);

  if (!count_only) {
    Rcpp::Rcout << "Generating dataframe\n";
    IntegerVector readIdx;
    IntegerVector bcIdx;
    auto prob = vector<double>();

    for (auto &x : all_res) {
      readIdx.push_back(std::get<0>(x));
      bcIdx.push_back(std::get<1>(x));
      prob.push_back(std::get<3>(x));
    }

    return List::create(List::create(Named("reads") = sequence_ids,
                                     Named("barcodes") = library_ids),
                        List::create(Named("i") = readIdx, Named("j") = bcIdx,
                                     Named("x") = prob,
                                     Named("index1") = false));
  }

  return List::create(Named("reads") = sequence_ids,
                      Named("barcodes") = library_ids);
}

//[[Rcpp::export]]
SEXP CRISPR_matching_DNAString(
    Rcpp::StringVector readSeq, Rcpp::StringVector readSeq_ids,
    Rcpp::StringVector readPhred, Rcpp::StringVector libSeq,
    Rcpp::StringVector libSeq_ids, String outFile, int misMatch,
    Rcpp::StringVector tMatSeq, Rcpp::NumericVector tMatProb, int numThread,
    bool hamming, bool count_only, double gap_left, double ext_left,
    double gap_right, double ext_right, double pen_max, bool detail_info) {
  // read in fastq file and store sequence and phred score as vector string.
  auto sequence = vector<string>(), sequence_ids = vector<string>(),
       phredScore = vector<string>(), library = vector<string>(),
       library_ids = vector<string>();

  // read in reads and library barcode from input
  for (int i = 0; i < readSeq.size(); i++) {
    sequence[i] = readSeq[i];
    sequence_ids[i] = readSeq_ids[i];
    phredScore[i] = readPhred[i];
  }

  for (int i = 0; i < libSeq.size(); i++) {
    library[i] = libSeq[i];
    libSeq_ids[i] = libSeq_ids[i];
  }

  auto threads = std::vector<thread>();
  auto match = hamming ? alignmentH : alignment;

  Trie trie(gap_left, ext_left, gap_right, ext_right, pen_max);

  // read samples, library, and setup trie's tMat
  // if( !(readSamples(sampleFile, sequence, sequence_ids, phredScore) &&
  //      trie.setTMat(tMatSeq, tMatProb)                             &&
  //      readLibrary(library, library_ids, libFile)                    ) )
  if (trie.setTMat(tMatSeq, tMatProb)) {
    return R_NilValue;
  }

  auto countTable = vector<double>(library.size());

  // create trie with library elements
  trie.fromLibrary(library);

  mutex mut;
  vector<res_t> all_res;

  std::string table_out = outFile;
  std::ofstream table_stream(table_out + ".txt");

  // run alignment
  try {

    int nr_item = sequence.size();
    int nr_item_per_thread = ceil(nr_item * 1.0 / numThread);

    Rcpp::Rcout << "Running"
                << (hamming ? " hamming search" : " levenshtein search")
                << " with " << nr_item_per_thread << " sequences per thread in "
                << numThread << " threads" << endl;

    for (int i = 1; i < numThread; ++i) {
      int R = min((i + 1) * nr_item_per_thread, nr_item);
      threads.emplace_back(match, std::ref(trie), std::ref(sequence),
                           std::ref(phredScore), misMatch, std::ref(countTable),
                           i * nr_item_per_thread, R, std::ref(all_res),
                           std::ref(mut), count_only, std::ref(table_stream),
                           detail_info);
    }

    int R = min(nr_item_per_thread, nr_item);
    match(trie, sequence, phredScore, misMatch, countTable, 0, R, all_res, mut,
          count_only, table_stream, detail_info);

    for (auto &t : threads) {
      t.join();
    }

  } catch (Rcpp::exception &e) {
    Rcpp::Rcout << e.what() << endl;
    return R_NilValue;
  }

  if (detail_info)
    trie.count(all_res, countTable, table_stream);

  Rcpp::Rcout << "Compiling results\n";

  count2CSV(countTable, library, outFile);

  if (!count_only) {
    Rcpp::Rcout << "Generating dataframe\n";
    IntegerVector readIdx;
    IntegerVector bcIdx;
    auto prob = vector<double>();

    for (auto &x : all_res) {
      readIdx.push_back(std::get<0>(x));
      bcIdx.push_back(std::get<1>(x));
      prob.push_back(std::get<3>(x));
    }

    return List::create(List::create(Named("reads") = sequence_ids,
                                     Named("barcodes") = library_ids),
                        List::create(Named("i") = readIdx, Named("j") = bcIdx,
                                     Named("x") = prob,
                                     Named("index1") = false));
  }

  return List::create(Named("reads") = sequence_ids,
                      Named("barcodes") = library_ids);
}

//[[Rcpp::export]]
SEXP CRISPR_user_matching_DNAString(
    Rcpp::StringVector readSeq, Rcpp::StringVector readSeq_ids,
    Rcpp::StringVector readPhred, Rcpp::StringVector libSeq,
    Rcpp::StringVector libSeq_ids, String outFile, int misMatch,
    Rcpp::StringVector tMatSeq, Rcpp::NumericVector tMatProb, int numThread,
    bool count_only, double gap_left, double ext_left, double gap_right,
    double ext_right, double pen_max, Function tForm) {
  // read in fastq file and store sequence and phred score as vector string.
  auto sequence = vector<string>(), sequence_ids = vector<string>(),
       phredScore = vector<string>(), library = vector<string>(),
       library_ids = vector<string>();

  for (int i = 0; i < readSeq.size(); i++) {
    sequence[i] = readSeq[i];
    sequence_ids[i] = readSeq_ids[i];
    phredScore[i] = readPhred[i];
  }

  for (int i = 0; i < libSeq.size(); i++) {
    library[i] = libSeq[i];
    libSeq_ids[i] = libSeq_ids[i];
  }

  auto threads = std::vector<thread>();

  Trie trie(gap_left, ext_left, gap_right, ext_right, pen_max);

  // read samples, library, and setup trie's tMat
  // if( !(readSamples(sampleFile, sequence, sequence_ids, phredScore) &&
  //      trie.setTMat(tMatSeq, tMatProb)                             &&
  //      readLibrary(library, library_ids, libFile)                    ) )
  if (trie.setTMat(tMatSeq, tMatProb)) {
    return R_NilValue;
  }

  auto countTable = vector<double>(library.size());

  // create trie with library elements
  trie.fromLibrary(library);

  mutex mut;
  vector<res_t> all_res;

  // run alignment
  try {

    int nr_item = sequence.size();
    int nr_item_per_thread = ceil(nr_item * 1.0 / numThread);

    Rcpp::Rcout << "Running"
                << " levenshtein search with " << nr_item_per_thread
                << " sequences per thread in " << numThread << " threads"
                << endl;

    for (int i = 1; i < numThread; ++i) {
      int R = min((i + 1) * nr_item_per_thread, nr_item);
      threads.emplace_back(user_alignment, std::ref(trie), std::ref(sequence),
                           std::ref(phredScore), misMatch, std::ref(countTable),
                           i * nr_item_per_thread, R, std::ref(all_res),
                           std::ref(mut), count_only);
    }

    int R = min(nr_item_per_thread, nr_item);
    user_alignment(trie, sequence, phredScore, misMatch, countTable, 0, R,
                   all_res, mut, count_only);

    for (auto &t : threads) {
      t.join();
    }

  } catch (Rcpp::exception &e) {
    Rcpp::Rcout << e.what() << endl;
    return R_NilValue;
  }

  clean(trie, all_res, tForm);

  trie.count(all_res, countTable);

  Rcpp::Rcout << "Compiling results\n";

  count2CSV(countTable, library, outFile);

  if (!count_only) {
    Rcpp::Rcout << "Generating dataframe\n";
    IntegerVector readIdx;
    IntegerVector bcIdx;
    auto prob = vector<double>();

    for (auto &x : all_res) {
      readIdx.push_back(std::get<0>(x));
      bcIdx.push_back(std::get<1>(x));
      prob.push_back(std::get<3>(x));
    }

    return List::create(List::create(Named("reads") = sequence_ids,
                                     Named("barcodes") = library_ids),
                        List::create(Named("i") = readIdx, Named("j") = bcIdx,
                                     Named("x") = prob,
                                     Named("index1") = false));
  }

  return List::create(Named("reads") = sequence_ids,
                      Named("barcodes") = library_ids);
}
