#ifndef __Node_H__
#define __Node_H__

#include <memory>
#include <mutex>
#include <string>
#include <vector>

using namespace std;

struct Node {
  array<int, 4> children;
  int leaf;
  char baseType;

  Node() : children{-1, -1, -1, -1}, leaf(0), baseType('N') {}
  
  Node(const char &base, const int l = 0) : children{-1, -1, -1, -1}, leaf(l), baseType(base) {}

  ~Node(){}
  
  // leaf
  bool isLeaf() const { return leaf != -1; }

  int getChild(char x) const {
    switch (x) {
      case 'A':
        return children[0];
        break;
      case 'C':
        return children[1];
        break;
      case 'T':
        return children[2];
        break;
      case 'G':
        return children[3];
        break;
    }

    return children[0];
  }

  
  void addChild(char x, int idx) {
    switch(x) {
      case 'A':
	children[0] = idx;
	break;
      case 'C':
        children[1] = idx;
        break;
      case 'T':
        children[2] = idx;
	break;
      case 'G':
	children[3] = idx;
	break;
    }
  }
};


#endif
