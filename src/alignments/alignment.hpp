#ifndef _ALIGNMENT_HPP
#define _ALIGNMENT_HPP

#include <Rcpp.h>
#include <array>
#include <iostream>
#include <memory>

#include "SA.hpp"

/*
*  Match + MisMatch
*/
class Match : public AlignmentBase<Match> {
  bool match;
  double phred;
  double mult;

public:
  auto value() -> double;

  auto penalty(std::array<double, 4> &pen) -> double;

  Match(std::shared_ptr<SA> &last, bool match, double phred, double mult)
      : AlignmentBase<Match>(last), match(match), phred(phred), mult(mult) {}

  Match(bool match, double phred, double mult)
      : AlignmentBase<Match>(), match(match), phred(phred), mult(mult) {}

  ~Match() {}
};

class Insertion : public AlignmentBase<Insertion> {
public:
  auto value() -> double;

  auto penalty(std::array<double, 4> &pen) -> double;

  Insertion(std::shared_ptr<SA> &last) : AlignmentBase<Insertion>(last) {}

  Insertion() {}
};

class I_Extension : public AlignmentBase<I_Extension> {
public:
  auto value() -> double;

  auto penalty(std::array<double, 4> &pen) -> double;

  I_Extension(std::shared_ptr<SA> &last) : AlignmentBase<I_Extension>(last) {}

  ~I_Extension() {}
};

class Deletion : public AlignmentBase<Deletion> {
public:
  auto value() -> double;

  auto penalty(std::array<double, 4> &pen) -> double;

  Deletion(std::shared_ptr<SA> &last) : AlignmentBase<Deletion>(last) {}

  Deletion() {}

  ~Deletion() {}
};

class D_Extension : public AlignmentBase<D_Extension> {
public:
  auto value() -> double;

  auto penalty(std::array<double, 4> &pen) -> double;

  D_Extension(std::shared_ptr<SA> &last) : AlignmentBase<D_Extension>(last) {}

  ~D_Extension() {}
};

#endif
