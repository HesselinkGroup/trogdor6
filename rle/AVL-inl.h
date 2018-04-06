/*
 *  AVL.cpp
 *  rle1d
 *
 *  Created by Paul Hansen on 11/3/09.
 *  Copyright 2009 Stanford University. All rights reserved.
 *
 *  This file is covered by the MIT license.  See LICENSE.txt.
 */

#ifdef _AVL_
#include <cassert>

namespace RLE
{

template <class T, class Compare>
AVLTree<T, Compare>::
AVLTree() :
    mSize(0),
    mRoot(0L)
{
}

template<class T, class Compare>
AVLTree<T, Compare>::
AVLTree(const AVLTree & copyMe) :
    mSize(0),
    mRoot(0L)
{
    ConstIterator itr;
    for (itr = copyMe.begin(); itr != copyMe.end(); itr++)
        insert(itr.peekNext());
}

template<class T, class Compare>
AVLTree<T,Compare> & AVLTree<T, Compare>::
operator=(const AVLTree & copyMe)
{
    if (&copyMe == this)
        return *this;
    
    delete mRoot;
    mRoot = 0L;
    mSize = 0;
    
    ConstIterator itr;
    for (itr = copyMe.begin(); itr != copyMe.end(); itr++)
        insert(itr.peekNext());
    
    return *this;
}


template <class T, class Compare>
AVLTree<T, Compare>::
~AVLTree()
{
    if (mRoot != 0L)
        delete mRoot;
}

template<class T, class Compare>
void AVLTree<T, Compare>::
clear()
{
    if (mRoot != 0L)
        delete mRoot;
    mRoot = 0L;
    mSize = 0;
}

template <class T, class Compare>
typename AVLTree<T,Compare>::Iterator AVLTree<T, Compare>::
insert(T key)
{
    Node* itrNextNode;
    if (mRoot == 0L)
    {
        mRoot = new Node(key);
        mSize++;
        itrNextNode = mRoot;
    }
    else
    {
        Node* current = mRoot;
        
        bool done = 0;
        while (!done)
        {
            if (mLessThan(key, current->key()))
            {
                if (current->left())
                    current = current->left();
                else
                {
                    itrNextNode = new Node(key);
                    current->setLeft(itrNextNode);
                    current->left()->setParent(current);
                    mSize++;
                    current->fixHeight();
                    rebalanceInsert(current);
                    done = 1;                    
                    //validate();
                }
            }
            else if (mLessThan(current->key(), key))
            {
                if (current->right())
                    current = current->right();
                else
                {
                    itrNextNode = new Node(key);
                    current->setRight(itrNextNode);
                    current->right()->setParent(current);
                    mSize++;
                    current->fixHeight();
                    rebalanceInsert(current);
                    done = 1;
                    //validate();
                }
            }
            else // no insertion for already-present key
            {
                //assert(key == current->key());
                done = 1;
                itrNextNode = current;
            }
        }
    }
    //validate();
        
    return Iterator(this, itrNextNode);
}


template <class T, class Compare>
void AVLTree<T, Compare>::
erase(T key)
{
    if (mRoot == 0L)
        return;
    else
    {
        Node* current = mRoot;
        bool done = 0;
        while (!done)
        {
            if (mLessThan(key, current->key()))
            {
                if (current->left())
                    current = current->left();
                else
                    done = 1; // key not present
            }
            else if (mLessThan(current->key(), key))
            {
                if (current->right())
                    current = current->right();
                else
                    done = 1;
            }
            else // key present: delete current node
            {
                assert(key == current->key());
                erase(current);
                done = 1;
            }
        }
    }
    //validate();
}

template <class T, class Compare>
bool AVLTree<T, Compare>::
count(T key) const
{
    Node* current = mRoot;
    while (current != 0L)
    {
        if (mLessThan(key, current->key()))
            current = current->left();
        else if (mLessThan(current->key(), key))
            current = current->right();
        else
        {
            assert(key == current->key());
            return 1;
        }
    }
    return 0;
}

template <class T, class Compare>
void AVLTree<T, Compare>::
validate() const
{
    if (mRoot != 0L)
        validate(mRoot);
}

template <class T, class Compare>
std::uint64_t AVLTree<T, Compare>::
bytes() const
{
    return sizeof(Node)*size();
}

template <class T, class Compare>
typename AVLTree<T, Compare>::Iterator AVLTree<T, Compare>::
begin()
{
    Iterator itr(this);
    
    Node* n = mRoot;
    if (n != 0L)
    {
        n = n->leftmost();
        itr.mNext = n;
    }
    
    return itr;
}

template <class T, class Compare>
typename AVLTree<T, Compare>::Iterator AVLTree<T, Compare>::
end()
{
    Iterator itr(this);
    
    Node* n = mRoot;
    if (n != 0L)
    {
        n = n->rightmost();
        itr.mPrev = n;
    }
    
    return itr;
}

template <class T, class Compare>
typename AVLTree<T, Compare>::ConstIterator AVLTree<T, Compare>::
begin() const
{
    ConstIterator itr(this);
    
    Node* n = mRoot;
    if (n != 0L)
    {
        n = n->leftmost();
        itr.mNext = n;
    }
    
    return itr;
}

template <class T, class Compare>
typename AVLTree<T, Compare>::ConstIterator AVLTree<T, Compare>::
end() const
{
    ConstIterator itr(this);
    
    Node* n = mRoot;
    if (n != 0L)
    {
        n = n->rightmost();
        itr.mPrev = n;
    }
    
    return itr;
}

template <class T, class Compare>
template<class S>
typename AVLTree<T, Compare>::Iterator AVLTree<T, Compare>::
location(const S & key)
{
    Node* current = mRoot;
    while (current != 0L)
    {
        if (mLessThan(key, current->key()))
        {
            if (current->left())
                current = current->left();
            else
            {
                Iterator iteratorSmaller(this);
                iteratorSmaller.mNext = current;
                iteratorSmaller.mPrev = current->previous();
                return iteratorSmaller;
            }
        }
        else if (mLessThan(current->key(), key))
        {
            if (current->right())
                current = current->right();
            else
            {
                Iterator iteratorLarger(this);
                iteratorLarger.mPrev = current;
                iteratorLarger.mNext = current->next();
                return iteratorLarger;
            }
        }
        else
        {
            //assert(key == current->key());
            //return 1;
            Iterator iteratorFound(this);
            iteratorFound.mNext = current;
            iteratorFound.mPrev = current->previous();
            return iteratorFound;
        }
    }
    //std::cerr << "Function failed to return.\n";
    return Iterator(this);
}

template <class T, class Compare>
template<class S>
typename AVLTree<T, Compare>::ConstIterator AVLTree<T, Compare>::
location(const S & key) const
{
    Node* current = mRoot;
    while (current != 0L)
    {
        if (mLessThan(key, current->key()))
        {
            if (current->left())
                current = current->left();
            else
            {
                ConstIterator iteratorSmaller(this);
                iteratorSmaller.mNext = current;
                iteratorSmaller.mPrev = current->previous();
                return iteratorSmaller;
            }
        }
        else if (mLessThan(current->key(), key))
        {
            if (current->right())
                current = current->right();
            else
            {
                ConstIterator iteratorLarger(this);
                iteratorLarger.mPrev = current;
                iteratorLarger.mNext = current->next();
                return iteratorLarger;
            }
        }
        else
        {
            //assert(key == current->key());
            //return 1;
            ConstIterator iteratorFound(this);
            iteratorFound.mNext = current;
            iteratorFound.mPrev = current->previous();
            return iteratorFound;
        }
    }
    //std::cerr << "Function failed to return.\n";
    return ConstIterator(this);
}



template <class T, class Compare>
void AVLTree<T, Compare>::
validate(AVLTree<T, Compare>::Node* n) const
{
    assert(n != 0L);
    if (n->balanceFactor() < -1 || n->balanceFactor() > 1)
        std::cerr << "Balance error\n";
    
    int lHeight = n->left() ? n->left()->height() : -1;
    int rHeight = n->right() ? n->right()->height() : -1;
    int myHeight = lHeight > rHeight ? lHeight+1 : rHeight+1;
    if (myHeight != n->height())
        std::cerr << "Height error for item " << n->key() << std::endl;
    
    if (n->left() != 0L)
    {
        if (n->left()->parent() != n)
            std::cerr << "Left child " << n->left()->key() <<
                " doesn't recognize parent " << n->key() << std::endl;
        if (n->key() <= n->left()->key())
            std::cerr << "Left child " << n->left()->key() << " out of order with "
                << n->key() << std::endl;
        if (n->height() <= n->left()->height())
            std::cerr << "Height error\n";
        validate(n->left());
    }
    if (n->right() != 0L)
    {
        if (n->right()->parent() != n)
            std::cerr << "Right child " << n->right()->key() <<
                " doesn't recognize parent "<< n->key() << std::endl;
        if (n->key() >= n->right()->key())
            std::cerr << "Right child " << n->right()->key() << " out of order with "
                << n->key() << std::endl;
        if (n->height() <= n->right()->height())
            std::cerr << "Height error\n";
        validate(n->right());
    }
}

// On insertion, rebalancing may stop after the first rotation.
template <class T, class Compare>
void AVLTree<T, Compare>::
rebalanceInsert(Node* node)
{
    assert(node);
    if (node->balanceFactor() == 2)
    {
        assert(node->left());
        if (node->left()->balanceFactor() == -1)
            rotateLR(node);
        else // left balance factor 1 or 0
        {
            assert(node->left()->balanceFactor() == 1);
            rotateLL(node);
        }
    }
    else if (node->balanceFactor() == -2)
    {
        assert(node->right());
        if (node->right()->balanceFactor() == 1)
            rotateRL(node);
        else
        {
            assert(node->right()->balanceFactor() == -1);
            rotateRR(node);
        }
    }
    else if (node->parent())
    {
        node->parent()->fixHeight();
        rebalanceInsert(node->parent());
    }
}

// On deletion, rebalancing may stop after the balance factor becomes +/- 1.
template <class T, class Compare>
void AVLTree<T, Compare>::
rebalanceErase(Node* node)
{
    assert(node);
    while (node != 0L)
    {
        Node* parent = node->parent();
        if (node->balanceFactor() == 2)
        {
            assert(node->left());
            if (node->left()->balanceFactor() == -1)
                rotateLR(node);
            else // left balance factor 1 or 0
                rotateLL(node);
        }
        else if (node->balanceFactor() == -2)
        {
            assert(node->right());
            if (node->right()->balanceFactor() == 1)
                rotateRL(node);
            else
                rotateRR(node);
        }

        // Actually recursion can stop if the node height is +1 or -1, apparently.
        if (parent)
        {
            parent->fixHeight();
            if (abs(parent->balanceFactor()) != 1)
                node = parent;
            else
                node = 0L;
        }
        else
            node = 0L;
    }
}

template <class T, class Compare>
void AVLTree<T, Compare>::
rotateLL(Node* node)
{
    assert(node->left());
    Node* left = node->left();
    
    // 1.  Connection with parent of node
    left->setParent(node->parent());
    if (node->parent())
    {
        if (node->parent()->left() == node)
            node->parent()->setLeft(left);
        else
            node->parent()->setRight(left);
    }
    else
        mRoot = left;
    
    // 2.  LR grandchild moves to right hand side
    node->setLeft(left->right());
    if (left->right())
        left->right()->setParent(node);
    
    // 3.  Left becomes parent of node
    left->setRight(node);
    node->setParent(left);
    
    node->fixHeight();
    left->fixHeight();
    
    //validate(left);
}

template <class T, class Compare>
void AVLTree<T, Compare>::
rotateRR(Node* node)
{
    assert(node->right());
    Node* right = node->right();
    
    // 1.  Connection with parent
    right->setParent(node->parent());
    if (node->parent())
    {
        if (node->parent()->left() == node)
            node->parent()->setLeft(right);
        else
            node->parent()->setRight(right);
    }
    else
        mRoot = right;
    
    // 2.  RL grandchild moves to left hand side
    node->setRight(right->left());
    if (right->left())
        right->left()->setParent(node);
    
    // 3.  Right becomes parent of node
    right->setLeft(node);
    node->setParent(right);
    
    node->fixHeight();
    right->fixHeight();
    
    //validate(right);
}

template <class T, class Compare>
void AVLTree<T, Compare>::
rotateLR(Node* node)
{
    assert(node->left());
    assert(node->left()->right());
    
    Node* left = node->left();
    Node* lr = left->right();
    
    // 1.  Connect to parent
    lr->setParent(node->parent());
    if (node->parent())
    {
        if (node == node->parent()->left())
            node->parent()->setLeft(lr);
        else
            node->parent()->setRight(lr);
    }
    else
        mRoot = lr;
    
    // 2.  LR's left becomes left's right
    left->setRight(lr->left());
    if (lr->left())
        lr->left()->setParent(left);
    
    // 3.  LR's right becomes node's left
    node->setLeft(lr->right());
    if (lr->right())
        lr->right()->setParent(node);
    
    // 4.  Left becomes left child of LR
    left->setParent(lr);
    lr->setLeft(left);
    
    // 5.  LR becomes parent of node
    lr->setRight(node);
    node->setParent(lr);
    
    node->fixHeight();
    left->fixHeight();
    lr->fixHeight();
    
    //validate(lr);
}

template <class T, class Compare>
void AVLTree<T, Compare>::
rotateRL(Node* node)
{
    assert(node->right());
    assert(node->right()->left());
    
    Node* right = node->right();
    Node* rl = right->left();
    
    // 1.  Connect to parent
    rl->setParent(node->parent());
    if (node->parent())
    {
        if (node == node->parent()->left())
            node->parent()->setLeft(rl);
        else
            node->parent()->setRight(rl);
    }
    else
        mRoot = rl;
    
    // 2.  RL's right becomes right's left
    right->setLeft(rl->right());
    if (rl->right())
        rl->right()->setParent(right);
    
    // 3.  RL's left becomes node's right
    node->setRight(rl->left());
    if (rl->left())
        rl->left()->setParent(node);
    
    // 4.  Right becomes right child of RL
    rl->setRight(right);
    right->setParent(rl);
    
    // 5.  RL becomes parent of node
    node->setParent(rl);
    rl->setLeft(node);
    
    node->fixHeight();
    right->fixHeight();
    rl->fixHeight();
    
    //validate(rl);
}

template <class T, class Compare>
void AVLTree<T, Compare>::
print(std::ostream & str) const
{
    if (mRoot)
        mRoot->print(str, 0);
}

template <class T, class Compare>
void AVLTree<T, Compare>::
erase(typename AVLTree<T, Compare>::Node* node)
{
    if (node->left()) // replace with rightmost on left
    {
        Node* replacer = node->left()->rightmost();
        assert(replacer->right() == 0L);
        node->setKey(replacer->key());
        if (replacer == replacer->parent()->right())
            replacer->parent()->setRight(replacer->left());
        else
            replacer->parent()->setLeft(replacer->left());
        
        if (replacer->left())
        {
            replacer->left()->setParent(replacer->parent());
            replacer->setLeft(0L);
        }
        
        if (replacer->parent())
        {
            replacer->parent()->fixHeight();
            rebalanceErase(replacer->parent());
        }
        delete replacer;
    }
    else if (node->right()) // replace with leftmost on right
    {
        Node* replacer = node->right()->leftmost();
        assert(replacer->left() == 0L);
        node->setKey(replacer->key());
        if (replacer == replacer->parent()->left())
            replacer->parent()->setLeft(replacer->right());
        else
            replacer->parent()->setRight(replacer->right());
        
        if (replacer->right())
        {
            replacer->right()->setParent(replacer->parent());
            replacer->setRight(0L);
        }
        
        if (replacer->parent())
        {
            replacer->parent()->fixHeight();
            rebalanceErase(replacer->parent());
        }
        
        delete replacer;
    }
    else // current is a leaf node
    {
        if (node->parent())
        {
            if (node == node->parent()->left())
                node->parent()->setLeft(0L);
            else
                node->parent()->setRight(0L);
        }
        else
            mRoot = 0L;
        
        Node* parent = node->parent();
        if (parent != 0L)
        {
            parent->fixHeight();
            rebalanceErase(parent);
        }
        
        delete node;
    }
    mSize--;
}



#pragma mark *** Node ***


template <class T, class Compare>
AVLTree<T, Compare>::Node::
Node(T key) :
    mKey(key),
    mHeight(0),
    mLeft(0L),
    mRight(0L),
    mParent(0L)
{
}

template <class T, class Compare>
AVLTree<T, Compare>::Node::
~Node()
{
    if (mLeft != 0L)
        delete mLeft;
    if (mRight != 0L)
        delete mRight;
}

template <class T, class Compare>
int AVLTree<T, Compare>::Node::
balanceFactor() const
{
    int leftHeight, rightHeight;
    leftHeight = mLeft ? mLeft->height()+1 : 0;
    rightHeight = mRight ? mRight->height()+1 : 0;
    return (leftHeight - rightHeight);
}

template <class T, class Compare>
void AVLTree<T, Compare>::Node::
fixHeight()
{
    mHeight = 0;
    if (mLeft)
        mHeight = mLeft->height() + 1;
    if (mRight && mRight->height() >= mHeight)
        mHeight = mRight->height() + 1;
}

template <class T, class Compare>
typename AVLTree<T, Compare>::Node* AVLTree<T, Compare>::Node::
leftmost() const
{
    const Node* current = this;
    while (current->left())
        current = current->left();
    return const_cast<Node*>(current);
}

template <class T, class Compare>
typename AVLTree<T, Compare>::Node* AVLTree<T, Compare>::Node::
rightmost() const
{
    const Node* current = this;
    while (current->right())
        current = current->right();
    return const_cast<Node*>(current);
}

template <class T, class Compare>
typename AVLTree<T, Compare>::Node* AVLTree<T, Compare>::Node::
next() const
{
    if (mRight)
        return mRight->leftmost();
    else
    {   
        const Node* current = this;
        
        while (current)
        {
            if (current->parent() && current == current->parent()->left())
                return current->parent();
            else
                current = current->parent();
        }
    }
    return 0L;
}

template <class T, class Compare>
typename AVLTree<T, Compare>::Node* AVLTree<T, Compare>::Node::
previous() const
{
    if (mLeft)
        return mLeft->rightmost();
    else
    {
        const Node* current = this;
        
        while (current)
        {
            if (current->parent() && current == current->parent()->right())
                return current->parent();
            else
                current = current->parent();
        }
    }
    return 0L;
}

template <class T, class Compare>
void AVLTree<T, Compare>::Node::
print(std::ostream & str, int indent) const
{
    if (mRight != 0L)
        mRight->print(str, indent+1);
    for (int nn = 0; nn < indent; nn++)
        str << "\t";
    str << mKey << " (" << mHeight << ") [" << balanceFactor() << "]\n";
    if (mLeft != 0L)
        mLeft->print(str, indent+1);
}

#pragma mark *** Iterator ***

template <class T, class Compare>
AVLTree<T, Compare>::Iterator::
Iterator(AVLTree<T, Compare>* tree, Node* nextNode) :
    mNext(nextNode),
    mPrev(nextNode->previous()),
    mTree(tree)
{
}

template <class T, class Compare>
AVLTree<T, Compare>::Iterator::
Iterator() :
    mNext(0L),
    mPrev(0L),
    mTree(0L)
{
}

template <class T, class Compare>
AVLTree<T, Compare>::Iterator::
Iterator(AVLTree<T, Compare>* tree) :
    mNext(0L),
    mPrev(0L),
    mTree(tree)
{
}

template <class T, class Compare>
typename AVLTree<T, Compare>::Iterator AVLTree<T, Compare>::Iterator::
operator++(int unused)
{
    Iterator returnValue(*this);
    
    if (mNext)
    {
        mPrev = mNext;
        mNext = mNext->next();
    }
    else
    {
        std::cerr << "At end of tree.\n";
    }
    
    return returnValue;
}

template <class T, class Compare>
typename AVLTree<T, Compare>::Iterator & AVLTree<T, Compare>::Iterator::
operator++()
{
    if (mNext)
    {
        mPrev = mNext;
        mNext = mNext->next();
    }
    else
    {
        std::cerr << "At end of tree.\n";
    }
    
    return *this;
}

template <class T, class Compare>
typename AVLTree<T, Compare>::Iterator AVLTree<T, Compare>::Iterator::
operator--(int unused)
{
    Iterator returnValue(*this);
    
    if (mPrev)
    {
        mNext = mPrev;
        mPrev = mPrev->previous();
    }
    else
    {
        std::cerr << "At start of tree.\n";
    }
    
    return returnValue;
}

template <class T, class Compare>
typename AVLTree<T, Compare>::Iterator & AVLTree<T, Compare>::Iterator::
operator--()
{
    if (mPrev)
    {
        mNext = mPrev;
        mPrev = mPrev->previous();
    }
    else
    {
        std::cerr << "At start of tree.\n";
    }
    
    return *this;
}

template<class T, class Compare>
typename AVLTree<T, Compare>::Iterator AVLTree<T, Compare>::Iterator::
nextPosition() const
{
    Iterator itr(*this);
    itr++;
    return itr;
}

template<class T, class Compare>
typename AVLTree<T, Compare>::Iterator AVLTree<T, Compare>::Iterator::
previousPosition() const
{
    Iterator itr(*this);
    itr--;
    return itr;
}

template <class T, class Compare>
void AVLTree<T, Compare>::Iterator::
seek(std::int64_t position)
{
    // while position is greater than next, go right
    // while position is less than or equal to prev, go left
    // uh, and !(b < a) is the same as a <= b.
    
    while (mNext && mLessThan(mNext->key(), position))
        operator++();
    while (mPrev && !mLessThan(mPrev->key(), position))
        operator--();
}

template <class T, class Compare>
bool AVLTree<T, Compare>::Iterator::
valid() const
{
    if (mTree == 0L)
        return 0;
    
    return (mNext != 0L || mPrev != 0L);
}

template <class T, class Compare>
T & AVLTree<T, Compare>::Iterator::
next()
{
    assert(mNext);
    (*this)++;
    return mPrev->key(); // since we just advanced, mPrev is the old mNext
}

template <class T, class Compare>
T & AVLTree<T, Compare>::Iterator::
previous()
{
    assert(mPrev);
    (*this)--;
    return mNext->key(); // since we just backed up, mNext is the old mPrev
}

template <class T, class Compare>
T & AVLTree<T, Compare>::Iterator::
peekNext() const
{
    assert(mNext);
    return mNext->key(); // since we just advanced, mPrev is the old mNext
}

template <class T, class Compare>
T & AVLTree<T, Compare>::Iterator::
peekPrevious() const
{
    assert(mPrev);
    return mPrev->key(); // since we just backed up, mNext is the old mPrev
}


template <class T, class Compare>
void AVLTree<T, Compare>::Iterator::
eraseNext()
{
    assert(mNext);
    assert(mTree);
    
    bool isLeaf = (mNext->right() == 0 && mNext->left() == 0);
    
    // Leaf nodes are deleted, but interior nodes are copied over.
    // I need to be careful about how I update mNext and mPrev, depending on
    // whether I'm deleting a leaf or not.
    if (isLeaf)
    {
        Node* nextNext = mNext->next();
        mTree->erase(mNext);
        mNext = nextNext;
    }
    else
    {
        T key = mNext->key();
        mTree->erase(mNext);
        // mNext is now a valid pointer to a replacement key node
        if (mLessThan(mNext->key(), key))
        {
            mPrev = mNext;
            mNext = mPrev->next();
        }
        else
        {
            mPrev = mNext->previous();
        }
    }
}

template <class T, class Compare>
void AVLTree<T, Compare>::Iterator::
erasePrevious()
{
    assert(mPrev);
    assert(mTree);
    
    bool isLeaf = (mPrev->right() == 0 && mPrev->left() == 0);
    
    // Leaf nodes are deleted, but interior nodes are copied over.
    // I need to be careful about how I update mNext and mPrev, depending on
    // whether I'm deleting a leaf or not.
    if (isLeaf)
    {
        Node* prevPrev = mPrev->previous();
        mTree->erase(mPrev);
        mPrev = prevPrev;
    }
    else
    {
        T key = mPrev->key();
        mTree->erase(mPrev);
        // mNext is now a valid pointer to a replacement key node
        if (mLessThan(key, mPrev->key()))
        {
            mNext = mPrev;
            mPrev = mNext->previous();
        }
        else
        {
            mNext = mPrev->next();
        }
    }
}

template <class T, class Compare>
bool AVLTree<T, Compare>::Iterator::
operator==(const Iterator & rhs) const
{
    return (mPrev == rhs.mPrev && mNext == rhs.mNext && mTree == rhs.mTree);
}



template <class T, class Compare>
bool AVLTree<T, Compare>::Iterator::
operator!=(const Iterator & rhs) const
{
    return !(mPrev == rhs.mPrev && mNext == rhs.mNext && mTree == rhs.mTree);
}


template <class T, class Compare>
void AVLTree<T, Compare>::Iterator::
print(std::ostream & str) const
{
    str << "[";
    if (mPrev)
        str << mPrev->key() << ",";
    else
        str << "_,";
    if (mNext)
        str << mNext->key();
    else
        str << "_";
    str << "]";
}

#pragma mark *** ConstIterator ***

template <class T, class Compare>
AVLTree<T, Compare>::ConstIterator::
ConstIterator() :
    mNext(0L),
    mPrev(0L),
    mTree(0L)
{
}

template <class T, class Compare>
AVLTree<T, Compare>::ConstIterator::
ConstIterator(const AVLTree<T, Compare>* tree) :
    mNext(0L),
    mPrev(0L),
    mTree(tree)
{
}

template <class T, class Compare>
AVLTree<T, Compare>::ConstIterator::
ConstIterator(const AVLTree<T, Compare>::Iterator & itr) :
    mNext(itr.mNext),
    mPrev(itr.mPrev),
    mTree(itr.mTree)
{
}


template <class T, class Compare>
typename AVLTree<T, Compare>::ConstIterator AVLTree<T, Compare>::ConstIterator::
operator++(int unused)
{
    ConstIterator returnValue(*this);
    
    if (mNext)
    {
        mPrev = mNext;
        mNext = const_cast<Node*>(mNext)->next();
    }
    else
    {
        std::cerr << "At end of tree.\n";
    }
    
    return returnValue;
}

template <class T, class Compare>
typename AVLTree<T, Compare>::ConstIterator & AVLTree<T, Compare>::ConstIterator::
operator++()
{
    if (mNext)
    {
        mPrev = mNext;
        mNext = mNext->next();
    }
    else
    {
        std::cerr << "At end of tree.\n";
    }
    
    return *this;
}

template <class T, class Compare>
typename AVLTree<T, Compare>::ConstIterator AVLTree<T, Compare>::ConstIterator::
operator--(int unused)
{
    ConstIterator returnValue(*this);
    
    if (mPrev)
    {
        mNext = mPrev;
        mPrev = const_cast<Node*>(mPrev->previous());
    }
    else
    {
        std::cerr << "At start of tree.\n";
    }
    
    return returnValue;
}

template <class T, class Compare>
typename AVLTree<T, Compare>::ConstIterator & AVLTree<T, Compare>::ConstIterator::
operator--()
{
    if (mPrev)
    {
        mNext = mPrev;
        mPrev = mPrev->previous();
    }
    else
    {
        std::cerr << "At start of tree.\n";
    }
    
    return *this;
}

template<class T, class Compare>
typename AVLTree<T, Compare>::ConstIterator AVLTree<T, Compare>::ConstIterator::
nextPosition() const
{
    ConstIterator itr(*this);
    itr++;
    return itr;
}

template<class T, class Compare>
typename AVLTree<T, Compare>::ConstIterator AVLTree<T, Compare>::ConstIterator::
previousPosition() const
{
    ConstIterator itr(*this);
    itr--;
    return itr;
}

template <class T, class Compare>
void AVLTree<T, Compare>::ConstIterator::
seek(std::int64_t position)
{
    // while position is greater than next, go right
    // while position is less than or equal to prev, go left
    // uh, and !(b < a) is the same as a <= b.
    
    while (mNext && mLessThan(mNext->key(), position))
        operator++();
    while (mPrev && !mLessThan(mPrev->key(), position))
        operator--();
}

template <class T, class Compare>
bool AVLTree<T, Compare>::ConstIterator::
valid() const
{
    if (mTree == 0L)
        return 0;
    
    return (mNext != 0L || mPrev != 0L);
}

template <class T, class Compare>
const T & AVLTree<T, Compare>::ConstIterator::
next()
{
    assert(mNext);
    (*this)++;
    return mPrev->key(); // since we just advanced, mPrev is the old mNext
}

template <class T, class Compare>
const T & AVLTree<T, Compare>::ConstIterator::
previous()
{
    assert(mPrev);
    (*this)--;
    return mNext->key(); // since we just backed up, mNext is the old mPrev
}

template <class T, class Compare>
const T & AVLTree<T, Compare>::ConstIterator::
peekNext() const
{
    assert(mNext);
    return mNext->key(); // since we just advanced, mPrev is the old mNext
}

template <class T, class Compare>
const T & AVLTree<T, Compare>::ConstIterator::
peekPrevious() const
{
    assert(mPrev);
    return mPrev->key(); // since we just backed up, mNext is the old mPrev
}

template <class T, class Compare>
bool AVLTree<T, Compare>::ConstIterator::
operator==(const ConstIterator & rhs) const
{
    return (mPrev == rhs.mPrev && mNext == rhs.mNext && mTree == rhs.mTree);
}

template <class T, class Compare>
bool AVLTree<T, Compare>::ConstIterator::
operator!=(const ConstIterator & rhs) const
{
    return !(mPrev == rhs.mPrev && mNext == rhs.mNext && mTree == rhs.mTree);
}

template <class T, class Compare>
void AVLTree<T, Compare>::ConstIterator::
print(std::ostream & str) const
{
    str << "[";
    if (mPrev)
        str << mPrev->key() << ",";
    else
        str << "_,";
    if (mNext)
        str << mNext->key();
    else
        str << "_";
    str << "]";
}


template<class T, class Compare>
bool operator==(const AVLTree<T, Compare> & lhs,
    const AVLTree<T, Compare> & rhs)
{
    typename AVLTree<T, Compare>::ConstIterator left = lhs.begin();
    typename AVLTree<T, Compare>::ConstIterator right = rhs.begin();
    
    while (left != lhs.end() && right != rhs.end())
    {
        if (!(left.peekNext() == right.peekNext()))
            return 0;
        left++;
        right++;
    }
    if (left != lhs.end() || right != rhs.end())
        return 0;
    
    return 1;
}

}; // namespace RLE

#endif
