#include "alignments/alignment.hpp"


auto Match::value() -> double
{
  return val = last == nullptr ? ( match ? 1.0 - phred : phred * mult ) :
                      last->value() * ( match ? 1.0 - phred : phred * mult );
}

auto Match::penalty(std::array<double, 4>& pen) -> double
{
  return pty = last == nullptr ? 0 : last->penalty(pen);
}


auto Insertion::value() -> double
{
  return val = last == nullptr ? 1.0 :  last->value();
}

auto Insertion::penalty(std::array<double, 4>& pen) -> double
{
  return pty = last == nullptr ? pen[0] : last->penalty(pen) + pen[0];
}


auto I_Extension::value() -> double
{
  return val = last == nullptr ? 1.0 : last->value();
}

auto I_Extension::penalty(std::array<double, 4>& pen) -> double
{
  return pty = last == nullptr ? pen[1] : last->penalty(pen) + pen[1];
}


auto Deletion::value() -> double
{
  return val = last == nullptr ? 1.0 : last->value();
}

auto Deletion::penalty(std::array<double, 4>& pen) -> double
{
  return pty = last == nullptr ? pen[2] : last->penalty(pen) + pen[2];
}


auto D_Extension::value() -> double
{
  return val = last == nullptr ? 1.0 : last->value();
}

auto D_Extension::penalty(std::array<double, 4>& pen) -> double
{
  return pty = last == nullptr ? pen[3] : last->penalty(pen) + pen[3];
}


template<typename T>
auto AlignmentBase<T>::value() -> double
{
  if( !done_val )
  {
    static_cast<T*>(this)->value();
    done_val = true;
  }

  return val;
}

template<typename T>
auto AlignmentBase<T>::penalty(std::array<double, 4>& pen) -> double
{
  if( !done_pen )
  {
    static_cast<T*>(this)->penalty(pen);
    done_pen = true;
  }

  return pty;
}


void deleter_t::operator()(void* x)
{
  switch(tag)
  {
    case tag_t::Match:
      delete static_cast<Match*>(x);
      break;
    case tag_t::Insertion:
      delete static_cast<Insertion*>(x);
      break;
    case tag_t::Deletion:
      delete static_cast<Deletion*>(x);
      break;
    case tag_t::I_Extension:
      delete static_cast<I_Extension*>(x);
      break;
    case tag_t::D_Extension:
      delete static_cast<D_Extension*>(x);
      break;
    case tag_t::Nil:
      break;
  }
}


auto SA::value() -> double
{
  switch(tag){
    case tag_t::Match:
      return static_cast<AlignmentBase<Match>*>(alignment.get())->value();
      break;
    case tag_t::Insertion:
      return static_cast<AlignmentBase<Insertion>*>(alignment.get())->value();
      break;
    case tag_t::Deletion:
      return static_cast<AlignmentBase<Deletion>*>(alignment.get())->value();
      break;
    case tag_t::I_Extension:
      return static_cast<AlignmentBase<I_Extension>*>(alignment.get())->value();
      break;
    case tag_t::D_Extension:
      return static_cast<AlignmentBase<D_Extension>*>(alignment.get())->value();
      break;
  }
}

auto SA::penalty(std::array<double, 4>& pen) -> double
{
  switch(tag){
    case tag_t::Match:
      return static_cast<AlignmentBase<Match>*>(alignment.get())->penalty(pen);
      break;
    case tag_t::Insertion:
      return static_cast<AlignmentBase<Insertion>*>(alignment.get())->penalty(pen);
      break;
    case tag_t::Deletion:
      return static_cast<AlignmentBase<Deletion>*>(alignment.get())->penalty(pen);
      break;
    case tag_t::I_Extension:
      return static_cast<AlignmentBase<I_Extension>*>(alignment.get())->penalty(pen);
      break;
    case tag_t::D_Extension:
      return static_cast<AlignmentBase<D_Extension>*>(alignment.get())->penalty(pen);
      break;
  }
}
