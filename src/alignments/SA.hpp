#ifndef _SA_HPP
#define _SA_HPP

#include <memory>

struct SA;

/*
*  CRTP classes
*  with cached sub-values
*/
template <typename T> class AlignmentBase {
public:
  std::shared_ptr<SA> last;
  double val;
  double pty;
  bool done_val;
  bool done_pen;

  auto value() -> double;

  auto penalty(std::array<double, 4> &pen) -> double;

  AlignmentBase(std::shared_ptr<SA> &last)
      : last(last), done_val(false), done_pen(false) {}

  AlignmentBase() : last(nullptr), done_val(false), done_pen(false) {}

  ~AlignmentBase() {}
};

class Match;
class Insertion;
class Deletion;
class I_Extension;
class D_Extension;

enum class tag_t { Match, Insertion, Deletion, I_Extension, D_Extension, Nil };

struct deleter_t {
  tag_t tag;

  void operator()(void *x);

  deleter_t(tag_t t) : tag(t) {}

  ~deleter_t() {}
};

/*
*  holds alignments w/ erased type,
*  re-types them to call methods
*  Using so we can use standard containers
*/
struct SA {
  std::unique_ptr<void, std::function<void(void *)>> alignment;
  tag_t tag;

public:
  auto value() -> double;

  auto penalty(std::array<double, 4> &pen) -> double;

  template <typename T>
  SA(T *x, tag_t t) : alignment(static_cast<void *>(x), deleter_t(t)), tag(t) {}

  SA() : alignment(nullptr) {}

  SA(SA &&o) : alignment(std::move(o.alignment)), tag(o.tag) {}
};

template <typename T> struct get_tag {};

template <> struct get_tag<Match> { tag_t tag = tag_t::Match; };

template <> struct get_tag<Insertion> { tag_t tag = tag_t::Insertion; };

template <> struct get_tag<Deletion> { tag_t tag = tag_t::Deletion; };

template <> struct get_tag<I_Extension> { tag_t tag = tag_t::I_Extension; };

template <> struct get_tag<D_Extension> { tag_t tag = tag_t::D_Extension; };

template <typename T, typename... Args>
auto make_alignment(Args &&... args) -> std::shared_ptr<SA> {
  return std::make_shared<SA>(new T(args...), get_tag<T>().tag);
}

#endif
