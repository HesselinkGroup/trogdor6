/*
 *  AVL.h
 *  rle1d
 *
 *  Created by Paul Hansen on 11/3/09.
 *  Copyright 2009 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifndef _AVL_
#define _AVL_

#include <iostream>
#include <cstdint>

namespace RLE
{

class UniversalLessThan
{
public:
    template<class T, class S>
    bool operator() (const T & lhs, const S & rhs) const
    {
        return lhs < rhs;
    }
};

/**
 * Implementation of an AVL tree with some slightly dangerous iterators.  The
 * iterators allow changing the key of a node, which can invalidate the tree.
 * It is safe to change a node's key as int64 as it does not cross another node's
 * key.
 *
 * @author Paul C. Hansen
 */
template <class T, class Compare = UniversalLessThan>
class AVLTree
{
public:
    class Iterator;
    class ConstIterator;
    
    AVLTree();
    AVLTree(const AVLTree & copyMe);
    AVLTree & operator=(const AVLTree & rhs);
    
    ~AVLTree();
    void clear();
    
    Iterator insert(T key);
    void erase(T key);
    bool count(T key) const;
    void validate() const;
    void print(std::ostream & stream) const;
    std::uint64_t size() const { return mSize; }
    std::uint64_t bytes() const;
    
    class Node
    {
    public:
        Node(T key);
        ~Node();
        const T & key() const { return mKey; }
        T & key() { return mKey; } // unsafe access
        void setKey(T key) { mKey = key; }
        
        int balanceFactor() const;
        std::uint64_t height() const { return mHeight; }
        void setHeight(std::uint64_t height) { mHeight = height; }
        void fixHeight();
        
        Node* left() const { return mLeft; }
        Node* right() const { return mRight; }
        Node* parent() const { return mParent; }
        
        Node* leftmost() const;
        Node* rightmost() const;
        
        Node* next() const;
        Node* previous() const;
        
        void setLeft(Node* n) { mLeft = n; }
        void setRight(Node* n) { mRight = n; }
        void setParent(Node* n) { mParent = n; }
        void unlink() { mLeft = mRight = mParent = 0L; }
        
        void print(std::ostream & str, int indent) const;
    protected:
        T mKey;
        
        std::uint64_t mHeight;
        Node* mLeft;
        Node* mRight;
        Node* mParent;
    };
    
protected:
    void erase(Node* node); // returns 1

public:
    class Iterator
    {
    protected:
        Iterator(AVLTree<T, Compare>* tree, Node* nextNode);
    public:
        Iterator();
        Iterator(AVLTree<T, Compare>* tree);
        Iterator operator++(int unused);
        Iterator & operator++();
        Iterator operator--(int unused);
        Iterator & operator--();
        
        Iterator previousPosition() const;
        Iterator nextPosition() const;
        
        void seek(std::int64_t position);
        
        bool valid() const;
        bool hasNext() const { return mNext != 0L; }
        bool hasPrevious() const { return mPrev != 0L; }
        
        /**
         * "Unsafe" access to the next and previous items: if you change them
         * you may invalidate the tree, and there is nothing to protect against
         * that situation.
         */
        T & next();  // like Java: gets next element and moves iterator
        T & previous(); // like Java: gets previous element and moves iterator
        T & peekNext() const; // like Qt4 Java-style iterators: just access
        T & peekPrevious() const; // just access
        
        void eraseNext();
        void erasePrevious();
        
        bool operator==(const Iterator & rhs) const;
        bool operator!=(const Iterator & rhs) const;
        
        void print(std::ostream & str) const;
        friend class AVLTree<T, Compare>::ConstIterator;
    protected:
        mutable Node* mNext;
        mutable Node* mPrev;
        AVLTree* mTree;
        Compare mLessThan;
        
        friend class AVLTree<T, Compare>;
    };
    
    class ConstIterator
    {
    public:
        ConstIterator();
        ConstIterator(const typename AVLTree<T, Compare>::Iterator & itr);
        ConstIterator(const AVLTree<T, Compare>* tree);
        ConstIterator operator++(int unused);
        ConstIterator & operator++();
        ConstIterator operator--(int unused);
        ConstIterator & operator--();
        
        ConstIterator previousPosition() const;
        ConstIterator nextPosition() const;
        
        void seek(std::int64_t position);
        
        bool valid() const;
        bool hasNext() const { return mNext != 0L; }
        bool hasPrevious() const { return mPrev != 0L; }
        
        const T & next();  // like Java: gets next element and advances
        const T & previous(); // like Java: gets previous element and retreats
        const T & peekNext() const; // just access
        const T & peekPrevious() const; // just access
        
        bool operator==(const ConstIterator & rhs) const;
        bool operator!=(const ConstIterator & rhs) const;
        
        void print(std::ostream & str) const;
    protected:
        const Node* mNext;
        const Node* mPrev;
        const AVLTree* mTree;
        Compare mLessThan;
        
        friend class AVLTree<T, Compare>;
    };
    friend class Iterator;
    friend class ConstIterator;
    
    Iterator begin();
    Iterator end();
    
    ConstIterator begin() const;
    ConstIterator end() const;
    
    // If the key is present, the returned iterator->next() is key.
    // Otherwise, the iterator marks the place the key should be inserted.
    template<class S>
    Iterator location(const S & key);
    
    template<class S>
    ConstIterator location(const S & key) const;
        
protected:
    void validate(Node* node) const;
    
    // Recurse toward root, rotating as necessary.  node must have the correct
    // height already.
    void rebalanceInsert(Node* node);
    void rebalanceErase(Node* node);
    
    // Node's left child becomes parent
    void rotateLL(Node* node);
    
    // Node's right child becomes parent
    void rotateRR(Node* node);
    
    // Node's LR grandchild becomes parent
    void rotateLR(Node* node);
    
    // Node's RL grandchild becomes parent
    void rotateRL(Node* node);
    
    Compare mLessThan;
    std::uint64_t mSize;
    Node* mRoot;
};

template<class T, class Compare>
bool operator==(const AVLTree<T, Compare> & lhs,
    const AVLTree<T, Compare> & rhs);

}; // namespace RLE

#include "AVL-inl.h"

#endif
