#ifndef NODE_POOL_HPP
#define NODE_POOL_HPP


#include <cassert>
#include <memory>
#include "Node.h"

// Memory pool for Tries
// all elements constructed in-place in buffer
// buffer size changes only occur while building library
// and threaded portion of algorithm uses this as read-only structure
struct Pool {
  std::unique_ptr<Node[]> buff;
  std::size_t size;
  int count;

  int createNode(char base, int leaf) {
    if(static_cast<std::size_t>(count + 1) >= size) {
      resize();
    }

    auto idx = count;
    ++count;

    buff[idx].baseType = base;
    buff[idx].leaf = leaf;
      
    return idx;
  }

  const Node &operator[](int i) const {
    assert(i != -1);
    assert(static_cast<std::size_t>(i) < size);
    return buff[i];
  }
  
  Pool(const int n = 32)
    : buff(new Node[n]())
    , size(n)
    , count(0) 
  {}

  ~Pool(){}

  int add_child(int parent, char base, int leaf) {
    assert(static_cast<std::size_t>(parent) < size);
    int child = createNode(base, leaf);
    assert(static_cast<std::size_t>(child) < size);
    buff[parent].addChild(base, child);
    return child;
  }
  
private:
  void resize() {
    size *= 2;
    auto *x = new Node[size]();
    std::memcpy(x, buff.get(), count * sizeof(Node));
    buff.reset(x);
  }
};

#endif
