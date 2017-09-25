#ifndef __Node_H__
#define __Node_H__

#include <memory>
#include <mutex>
#include <string>
#include <vector>

using namespace std;

class Node {
    char baseType;
    bool leaf;
    array<unique_ptr<Node>, 4> children;

public:
    Node()
        : baseType('N')
        , leaf(false)
    {
    }

    Node(const char& base, const bool& l = false)
        : baseType(base)
        , leaf(l)
    {
    }

    Node(Node& o)
        : baseType(o.baseType)
        , leaf(o.leaf)
    {
    }

    virtual ~Node() {}

    // base
    char getBase() { return baseType; }
    // leaf
    bool isLeaf() { return leaf; }

    Node* getChild(char x)
    {
        int i;
        switch (x) {
        case 'A':
            i = 0;
            break;
        case 'C':
            i = 1;
            break;
        case 'T':
            i = 2;
            break;
        case 'G':
            i = 3;
            break;
        }

        return children[i].get();
    }

    Node* getChild(int i)
    {
        return children[i].get();
    }

    void addChild(std::unique_ptr<Node>& child)
    {
        int i;
        switch (child->baseType) {
        case 'A':
            i = 0;
            break;
        case 'C':
            i = 1;
            break;
        case 'T':
            i = 2;
            break;
        case 'G':
            i = 3;
            break;
        }

        children[i] = move(child);
    }
};

class Leaf : public Node {
public:
    const int idx;

    Leaf(char baseType, const int idx)
        : Node(baseType, true)
        , idx(idx)
    {
    }

    ~Leaf() {}
};

template <typename... Args>
static auto make_leaf(Args... args) -> unique_ptr<Node>
{
    return unique_ptr<Node>(dynamic_cast<Node*>(new Leaf(args...)));
}

template <typename... Args>
static auto make_node(Args... args) -> unique_ptr<Node>
{
    return unique_ptr<Node>(new Node(args...));
}

#endif
