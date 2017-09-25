#include "structures/Trie.h"
#include <thread>

// Key Trie class method: read search with hamming type error.
auto Trie::hammingSearch(state_t state) -> void
{
    // ---------update and check mismatch bound -------------------

    state.misMatch -= !state.match();

    if (state.misMatch < 0) {
        return;
    }

    state.add<Match>(state.match(), state.phred(), getMatEle(state.seq, state.seqIdx));

    // ---------------- check exit conditions ---------------------

    if (state.leaf() || state.end())
    {
        if (state.leaf() && state.end())
        {
            state.addResult();
        }

        return;
    }

    // ----------------- recursion -------------------------------

    auto s = state;
    s.seqIdx += 1;

    for (auto i = 0; i < 4; i++) {
        auto* child = state.child(i);
        if (child != nullptr) {
            s.cur = child;
            hammingSearch(s);
        }
    }
}
