#include "structures/Trie.h"

auto Trie::editSearch(state_t state, bool ind) -> void {
  // ------- exit conditions --------
  if (state.misMatch < 0 || state.pen > pen_max) {
    return;
  }

  if (state.leaf()) {
    if (state.end()) {
      state.addResult();
    } else if (!ind) {
      indel(state, true);
    }
    return;
  }

  if (state.end()) {
    if (!ind) {
      indel(state, false);
    }
    return;
  }

  // ------- not an exit state ---------

  if (!ind) {
    indel(state, true);
    indel(state, false);
  }

  auto s = state;
  s.misMatch -= !state.match();
  s.add<Match>(s.match(), s.phred(), getMatEle(s.seq, s.seqIdx));
  s.seqIdx += 1;

  for (auto i = 0; i < 4; i++) {
    s.cur = state.child(i);

    if (s.cur != -1) {
      editSearch(s, false);
    }
  }
}

auto Trie::indel(state_t state, bool left) -> void {
  state.misMatch -= 1;

  if (state.misMatch < 0 || state.pen > pen_max) {
    return;
  }

  if (left) // add gap in reference
  {
    state.seqIdx += 1;
    state.pen += _penalties[0];
    state.add<Insertion>();

    if (!state.leaf()) {
      editSearch(state, true);
    }

    extend(state, left);

  } else {
    auto s = state;
    s.pen += _penalties[2];
    s.add<Deletion>();

    for (auto i = 0; i < 4; ++i) // add gap in sequence
    {
      s.cur = state.child(i);
      if (s.cur != -1) {
        if (!state.end()) {
          editSearch(s, true);
        }

        extend(s, left);
      }
    }
  }
}

auto Trie::extend(state_t state, bool left) -> void {
  if (state.pen > pen_max) {
    return;
  }

  if (left) // add gap in reference
  {
    if (state.end()) {
      if (state.leaf()) {
        state.addResult();
      }

      return;
    }

    state.seqIdx += 1;
    state.pen += _penalties[1];
    state.add<I_Extension>();

    if (!state.leaf()) {
      editSearch(state, true);
    }

    extend(state, left);

  } else {
    if (state.leaf()) {
      if (state.end()) {
        state.addResult();
      }
      return;
    }

    auto s = state;
    s.pen += _penalties[3];
    s.add<D_Extension>();

    for (auto i = 0; i < 4; ++i) // add gap in sequence
    {
      s.cur = state.child(i);
      if (s.cur != -1) {
        if (!state.end()) {
          editSearch(s, true);
        }
        extend(s, left);
      }
    }
  }
}
