#ifndef TREETYPE_H
#define TREETYPE_H

#include <iostream>

using namespace std;

typedef int ItemType;

struct TreeNode
{
    ItemType info;
    TreeNode* left;
    TreeNode* right;
};

class TreeType {// binary search tree.
    public:
        TreeType() {
            root = NULL;
        }
        void CreateTree(ItemType *itemArray, int itemNum);
        void PreorderTraverse(TreeNode *rootNode);
        void InorderTraverse(TreeNode *rootNode);
        void PostorderTraverse(TreeNode *rootNode);
        void InsertItems(ItemType *itemArray, int itemNum);
        bool FindItem(TreeNode *rootNode, ItemType item);
        ItemType MinValue(TreeNode *rootNode);
        ItemType MaxValue(TreeNode *rootNode);
        int NodesNum(TreeNode *rootNode);
        void LeavesNum(TreeNode *rootNode, int &num);
        void PreorderTraverseLeaves(TreeNode *rootNode);
        void InorderTraverseLeaves(TreeNode *rootNode);
        void PostorderTraverseLeaves(TreeNode *rootNode);
        bool DeleteItem(TreeNode *rootNode, ItemType item);

        TreeNode *root;
};


#endif // TREETYPE_H
