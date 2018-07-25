#ifndef CSGTREE_H
#define CSGTREE_H

#include "hitable.h"
#include <iomanip>
using namespace std;

struct SolidStruct
{
    hitable *solid;
    bool hitted;
    float t, t2;
    vec3 normal, normal2;
    int operation;
/*
operation=0: the solid is primitive;
operation=1: union;
operation=2: intersection;
operation=3: difference1;
operation=4: difference2;
when operation is not 0, solid is NULL;
*/
};
typedef SolidStruct ItemType;
struct csgTreeNode
{
    ItemType info;
    csgTreeNode* left;
    csgTreeNode* right;
/*
when operation is 0, left and right are NULL;
*/
};

class csgTree : public hitable
{
    public:
        csgTree() {
            root = NULL;
        }
//        void CreateTree(ItemType *itemArray, int itemNum);
        virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
        csgTreeNode *root;
};

#endif // CSGTREE_H
