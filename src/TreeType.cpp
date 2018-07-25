#include "TreeType.h"

void TreeType::CreateTree(ItemType *itemArray, int itemNum) {
     for(int i=0; i<itemNum; i++) {
         TreeNode *newnode = new TreeNode;
         newnode->info = itemArray[i];
         newnode->left = NULL;
         newnode->right = NULL;
         if(root == NULL) {
             root = newnode;
             std::cout << itemArray[i] << " will be the root" << endl;
         }
         else {
             TreeNode *parent;
             TreeNode *temp = root;
             while(temp != NULL) {
                 parent = temp;
                 if((temp->info) > itemArray[i]) {
                     temp = temp->left;
                 }
                 else {
                     temp = temp->right;
                 }
             }

             if((parent->info) > itemArray[i]) {
                 parent->left = newnode;
                 std::cout << itemArray[i] << " will be " << parent->info << "'s left child" << endl;
             }
             else {
                 parent->right = newnode;
                 std::cout << itemArray[i] << " will be " << parent->info << "'s right child" << endl;
             }
         }
    }
}

 void TreeType::PreorderTraverse(TreeNode *rootNode) {
     if(rootNode!=NULL)
     {
         std::cout << rootNode->info << " ";
         PreorderTraverse(rootNode->left);
         PreorderTraverse(rootNode->right);
     }
 }

 void TreeType::InorderTraverse(TreeNode *rootNode) {
     if(rootNode!=NULL)
     {
         InorderTraverse(rootNode->left);
         std::cout << rootNode->info << " ";
         InorderTraverse(rootNode->right);
     }
 }

 void TreeType::PostorderTraverse(TreeNode *rootNode) {
     if(rootNode!=NULL)
     {
         PostorderTraverse(rootNode->left);
         PostorderTraverse(rootNode->right);
         std::cout << rootNode->info << " ";
     }
 }

void TreeType::InsertItems(ItemType *itemArray, int itemNum) {
     for(int i=0; i<itemNum; i++) {
         TreeNode *newnode = new TreeNode;
         newnode->info = itemArray[i];
         newnode->left = NULL;
         newnode->right = NULL;
         if(root == NULL) {
             root = newnode;
             std::cout << itemArray[i] << " will be the root" << endl;
         }
         else {
             TreeNode *parent;
             TreeNode *temp = root;
             while(temp != NULL) {
                 parent = temp;
                 if((temp->info) > itemArray[i]) {
                     temp = temp->left;
                 }
                 else if((temp->info) == itemArray[i]) {
                    break;
                 }
                 else {
                     temp = temp->right;
                 }
             }

             if((parent->info) > itemArray[i]) {
                 parent->left = newnode;
                 std::cout << itemArray[i] << " will be " << parent->info << "'s left child" << endl;
             }
             else if((parent->info) == itemArray[i]) {
                 std::cout << "there has been a " << itemArray[i] << " in the tree" << endl;
             }
             else {
                 parent->right = newnode;
                 std::cout << itemArray[i] << " will be " << parent->info << "'s right child" << endl;
             }
         }
    }
}

bool TreeType::FindItem(TreeNode *rootNode, ItemType item) {
    if (rootNode == NULL) {
        return  false;
    }
    else {
        if ((rootNode->info) == item) {
            std::cout << "found." << endl;
            if ((rootNode->left) != NULL) {
                std::cout << "its left child is " << rootNode->left->info << " ." << endl;
            }
            if ((rootNode->right) != NULL) {
                std::cout << "its right child is " << rootNode->right->info << " ." << endl;
            }
            return true;
        }
        else {
            if (FindItem(rootNode->left, item)) {
                return true;
            }
            if (FindItem(rootNode->right, item)) {
                return true;
            }
            return false;
        }
    }
}

ItemType TreeType::MinValue(TreeNode *rootNode) {
     if(rootNode!=NULL)
     {
         if ((rootNode->left) != NULL) {
            MinValue(rootNode->left);
         }
         else {
            return rootNode->info;
         }
     }
    return -1;//error
}

ItemType TreeType::MaxValue(TreeNode *rootNode) {
     if(rootNode!=NULL)
     {
         if ((rootNode->right) != NULL) {
            MaxValue(rootNode->right);
         }
         else {
            return rootNode->info;
         }
     }
    return -1;//error
}

int TreeType::NodesNum(TreeNode *rootNode) {
     if(rootNode!=NULL)
     {
        return (NodesNum(rootNode->left) + NodesNum(rootNode->right) + 1);
     }
     else {
        return 0;
     }
}

void TreeType::LeavesNum(TreeNode *rootNode, int &num) {
     if(rootNode!=NULL)
     {
        if (((rootNode->left) == NULL) && ((rootNode->right) == NULL)) {
            num ++;
        }
        else {
            LeavesNum(rootNode->left, num);
            LeavesNum(rootNode->right, num);
        }
     }
}


void TreeType::PreorderTraverseLeaves(TreeNode *rootNode) {
    if(rootNode!=NULL) {
        if ((rootNode->left == NULL) && (rootNode->right == NULL)) {
            std::cout << rootNode->info << " ";
        }
        PreorderTraverseLeaves(rootNode->left);
        PreorderTraverseLeaves(rootNode->right);
    }
}

void TreeType::InorderTraverseLeaves(TreeNode *rootNode) {
    if(rootNode!=NULL) {
        InorderTraverseLeaves(rootNode->left);
        if ((rootNode->left == NULL) && (rootNode->right == NULL)) {
            std::cout << rootNode->info << " ";
        }
        InorderTraverseLeaves(rootNode->right);
    }
}

void TreeType::PostorderTraverseLeaves(TreeNode *rootNode) {
    if(rootNode!=NULL)
    {
        PostorderTraverseLeaves(rootNode->left);
        PostorderTraverseLeaves(rootNode->right);
        if ((rootNode->left == NULL) && (rootNode->right == NULL)) {
            std::cout << rootNode->info << " ";
        }
    }
}

void DeleteTree(TreeNode *rootNode) {
    if (rootNode == NULL) {
    }
    else {
        TreeNode *temp_left = rootNode->left;
        TreeNode *temp_right = rootNode->right;
        delete rootNode;
        if (temp_left != NULL) {
            DeleteTree(temp_left);
        }
        if (temp_right != NULL) {
            DeleteTree(temp_right);
        }
    }
 }

 bool TreeType::DeleteItem(TreeNode *rootNode, ItemType item) {
    if (rootNode == NULL) {
        return  false;
    }
    else {
        if (rootNode->left != NULL) {
            if (rootNode->left->info == item) {
                DeleteTree(rootNode->left);
                rootNode->left = NULL;
                return true;
            }
        }
        if (rootNode->right != NULL) {
            if (rootNode->right->info == item) {
                DeleteTree(rootNode->right);
                rootNode->right = NULL;
                return true;
            }
        }
        if (DeleteItem(rootNode->left, item)) {
            return true;
        }
        if (DeleteItem(rootNode->right, item)) {
            return true;
        }
        else {
            return false;
        }
    }
 }
