#include <random>
#include "Mesh.h"
#include "Primitive.h"
#include "threading.hpp"

struct OctTree : IntersectionAccelerator {
    struct Node {
        BBox box;
        Node* children[8] = {
            nullptr,
        };
        std::vector<Intersectable*> primitives;
        bool isLeaf() const { return children[0] == nullptr; }
    };

    std::vector<Intersectable*> allPrimitives;
    Node* root = nullptr;
    int depth = 0;
    int leafSize = 0;
    int nodes = 0;
    int MAX_DEPTH = 35;
    int MIN_PRIMITIVES = 10;

    void clear(Node* n) {
        if (!n) {
            return;
        }

        for (int c = 0; c < 8; c++) {
            clear(n->children[c]);
            delete n->children[c];
        }
    }

    void clear() {
        clear(root);
        allPrimitives.clear();
    }

    void addPrimitive(Intersectable* prim) override { allPrimitives.push_back(prim); }

    void build(Node* n, int currentDepth = 0) {
        if (currentDepth >= MAX_DEPTH || n->primitives.size() <= MIN_PRIMITIVES) {
            leafSize = std::max(leafSize, int(n->primitives.size()));
            return;
        }

        depth = std::max(depth, currentDepth);
        BBox childBoxes[8];
        n->box.octSplit(childBoxes);

        for (int c = 0; c < 8; c++) {
            Node*& child = n->children[c];
            child = new Node;
            nodes++;
            memset(child->children, 0, sizeof(child->children));
            child->box = childBoxes[c];
            for (int r = 0; r < n->primitives.size(); r++) {
                if (n->primitives[r]->boxIntersect(child->box)) {
                    child->primitives.push_back(n->primitives[r]);
                }
            }
            if (child->primitives.size() == n->primitives.size()) {
                build(child, MAX_DEPTH + 1);
            } else {
                build(child, currentDepth + 1);
            }
        }
        n->primitives.clear();
    }

    void build(Purpose purpose) override {
        const char* treePurpose = "";
        if (purpose == Purpose::Instances) {
            MAX_DEPTH = 5;
            MIN_PRIMITIVES = 4;
            treePurpose = " instances";
        } else if (purpose == Purpose::Mesh) {
            MAX_DEPTH = 35;
            MIN_PRIMITIVES = 20;
            treePurpose = " mesh";
        }

        if (root) {
            clear(root);
            delete root;
        }

        printf("Building%s oct tree with %d primitives... ", treePurpose,
               int(allPrimitives.size()));
        Timer timer;
        nodes = leafSize = depth = 0;
        root = new Node();
        root->primitives.swap(allPrimitives);

        // std::vector<BBox> primBounds; /////////////////////////////
        // primBounds.reserve(root->primitives.size()); /////////////////
        for (int c = 0; c < root->primitives.size(); c++) {
            root->primitives[c]->expandBox(root->box);
            // primBounds.push_back(root->box); ////////////////////////
        }
        // int sizeP = root->primitives.size();
        // auto box0 = primBounds[0];//////////////////////////////////
        // auto boxLast = primBounds[root->primitives.size() - 1]; //////
        // auto randBox = primBounds[1456]; ////////////////////////////

        build(root);
        printf(" done in %lldms, nodes %d, depth %d, %d leaf size\n", timer.toMs(timer.elapsedNs()),
               nodes, depth, leafSize);
    }

    bool intersect(Node* n, const Ray& ray, float tMin, float& tMax, Intersection& intersection) {
        bool hasHit = false;

        if (n->isLeaf()) {
            for (int c = 0; c < n->primitives.size(); c++) {
                if (n->primitives[c]->intersect(ray, tMin, tMax, intersection)) {
                    tMax = intersection.t;
                    hasHit = true;
                }
            }
        } else {
            for (int c = 0; c < 8; c++) {
                if (n->children[c]->box.testIntersect(ray)) {
                    if (intersect(n->children[c], ray, tMin, tMax, intersection)) {
                        tMax = intersection.t;
                        hasHit = true;
                    }
                }
            }
        }

        return hasHit;
    }

    bool intersect(const Ray& ray, float tMin, float tMax, Intersection& intersection) override {
        return intersect(root, ray, tMin, tMax, intersection);
    }

    bool isBuilt() const override { return root != nullptr; }

    ~OctTree() override { clear(); }
};

////////////////////////////////////////////////////////////////////////////////////
enum class EdgeType { Start, End };

struct BoundEdge {
    // BoundEdge Public Methods
    BoundEdge() {}
    BoundEdge(float t, int primNum, bool starting) : t(t), primNum(primNum) {
        type = starting ? EdgeType::Start : EdgeType::End;
    }

    float t;
    int primNum;
    EdgeType type;
};

struct KDNode {
    void InitLeaf(int* primNums, int np, std::vector<int>* primitiveIndices) {
        flags = 3;
        numPrimitives |= (np << 2);
        // Store primitive ids for leaf node
        if (np == 0)
            onePrimitive = 0;
        else if (np == 1)
            onePrimitive = primNums[0];
        else {
            offsetIdx = primitiveIndices->size();
            for (int i = 0; i < np; ++i) primitiveIndices->push_back(primNums[i]);
        }
    }

    void InitInterior(int axis, int ac, float s) {
        splitAxis = s;
        flags = axis;
        aboveChild |= (ac << 2);
    }

    float SplitPos() const { return splitAxis; }
    int nPrimitives() const { return numPrimitives >> 2; }
    int SplitAxis() const { return flags & 3; }
    bool IsLeaf() const { return (flags & 3) == 3; }
    int AboveChild() const { return aboveChild >> 2; }

    union {
        float splitAxis;
        int onePrimitive;
        int offsetIdx;
    };

private:
    union {
        int flags;  // interior or leaf node
        int numPrimitives;
        int aboveChild;
    };
};

struct KdToDo {
    const KDNode* node;
    float tMin, tMax;
};

struct BuildTreeLoader {};

struct KDTree : IntersectionAccelerator {
    KDTree() = default;

    KDTree(const int _intersectCost, const int _traversalCost)
        : intersectCost(_intersectCost), traversalCost(_traversalCost) {}

    void build(Purpose purpose) override {
        const char* treePurpose = "";
        switch (purpose) {
            case Purpose::Instances:
                maxAllowedPrimitives = 5;
                treePurpose = " instances";
                break;
            case Purpose::Mesh:
                maxAllowedPrimitives = 16;
                treePurpose = " mesh";
                break;
            default:
                printf("Received unsupported purpose");
                break;
        }

        nAllocatedNodes = nextFreeNode = 0;

        if (maxDepth <= 0)
            maxDepth = std::round(8 + 1.3f * log2(int64_t(allPrimitives.size())));

        printf("Building%s KNTree tree with %d primitives... ", treePurpose,
               int(allPrimitives.size()));

        Timer timer;

        // Compute bounds for kd-tree construction
        std::vector<BBox> primBounds;
        primBounds.reserve(allPrimitives.size());
        for (const auto prim : allPrimitives) {
            BBox box = prim->worldBounds();
            prim->expandBox(bounds);
            primBounds.push_back(box);
        }

        // Initialize _primNums_ for kd-tree construction
        std::unique_ptr<int[]> primNums(new int[allPrimitives.size()]);
        for (size_t i = 0; i < allPrimitives.size(); ++i) primNums[i] = i;

        // Allocate working memory for kd-tree construction
        std::unique_ptr<BoundEdge[]> edges[3];
        for (int i = 0; i < 3; ++i) edges[i].reset(new BoundEdge[2 * allPrimitives.size()]);

        std::unique_ptr<int[]> lowerBoundPrims(new int[allPrimitives.size()]);
        std::unique_ptr<int[]> upperBoundPrims(new int[(maxDepth + 1) * allPrimitives.size()]);

        // Start recursive construction of kd-tree
        buildTree(0, bounds, primBounds, primNums.get(), allPrimitives.size(), maxDepth, edges,
                  lowerBoundPrims.get(), upperBoundPrims.get());

        printf("Num created nodes %d, done in %lldms, nodes allocated %d, depth %d\n",
               nextFreeNode - 1, timer.toMs(timer.elapsedNs()), nAllocatedNodes, maxDepth);
    }

private:
    void buildTree(int nodeNum, const BBox& nodeBounds, const std::vector<BBox>& allPrimBounds,
                   int* primNums, int nPrimitives, int depth,
                   const std::unique_ptr<BoundEdge[]> edges[3], int* prims0, int* prims1,
                   int badRefines = 0) {
        // Get next free node from _nodes_ array
        if (nextFreeNode == nAllocatedNodes) {
            int nNewAllocNodes = std::max(2 * nAllocatedNodes, 512);
            KDNode* n = new KDNode[nNewAllocNodes];
            if (nAllocatedNodes > 0) {
                memcpy(n, nodes, nAllocatedNodes * sizeof(KDNode));
                delete[] nodes;
            }
            nodes = n;
            nAllocatedNodes = nNewAllocNodes;
        }
        ++nextFreeNode;

        // Initialize leaf node if termination criteria met
        if (nPrimitives <= maxAllowedPrimitives || depth == 0) {
            nodes[nodeNum].InitLeaf(primNums, nPrimitives, &primitiveIndices);
            return;
        }

        // Choose split axis position for interior node
        int bestAxis = -1, bestOffset = -1;
        float bestCost = INFINITY;
        float oldCost = intersectCost * float(nPrimitives);
        vec3 d = nodeBounds.max - nodeBounds.min;
        float totalSA = 2.f * (d.x * d.y + d.x * d.z + d.y * d.z);
        float invTotalSA = 1 / totalSA;

        // Choose which axis to split along
        int axis;
        if (d.x > d.y && d.x > d.z)
            axis = 0;
        else
            axis = (d.y > d.z) ? 1 : 2;

        int retries = 0;

    retrySplit:
        // Initialize edges for _axis_
        for (int i = 0; i < nPrimitives; ++i) {
            int pn = primNums[i];
            const BBox& bounds = allPrimBounds[pn];
            edges[axis][2 * i] = BoundEdge(bounds.min[axis], pn, true);
            edges[axis][2 * i + 1] = BoundEdge(bounds.max[axis], pn, false);
        }

        // Sort _edges_ for _axis_
        std::sort(&edges[axis][0], &edges[axis][2 * nPrimitives],
                  [](const BoundEdge& e0, const BoundEdge& e1) -> bool {
                      if (e0.t == e1.t)
                          return (int)e0.type < (int)e1.type;
                      else
                          return e0.t < e1.t;
                  });

        // Compute cost of all splits for _axis_ to find best
        int nBelow = 0, nAbove = nPrimitives;
        for (int i = 0; i < 2 * nPrimitives; ++i) {
            if (edges[axis][i].type == EdgeType::End)
                --nAbove;
            float edgeT = edges[axis][i].t;
            if (edgeT > nodeBounds.min[axis] && edgeT < nodeBounds.max[axis]) {
                // Compute cost for split at _i_th edge

                // Compute child surface areas for split at _edgeT_
                int otherAxis0 = (axis + 1) % 3, otherAxis1 = (axis + 2) % 3;
                float belowSA =
                    2 * (d[otherAxis0] * d[otherAxis1] +
                         (edgeT - nodeBounds.min[axis]) * (d[otherAxis0] + d[otherAxis1]));
                float aboveSA =
                    2 * (d[otherAxis0] * d[otherAxis1] +
                         (nodeBounds.max[axis] - edgeT) * (d[otherAxis0] + d[otherAxis1]));

                float pBelow = belowSA * invTotalSA;
                float pAbove = aboveSA * invTotalSA;
                float eb = (nAbove == 0 || nBelow == 0) ? 0.5 : 0;
                float cost =
                    traversalCost + intersectCost * (1 - eb) * (pBelow * nBelow + pAbove * nAbove);

                // Update best split if this is lowest cost so far
                if (cost < bestCost) {
                    bestCost = cost;
                    bestAxis = axis;
                    bestOffset = i;
                }
            }
            if (edges[axis][i].type == EdgeType::Start)
                ++nBelow;
        }

        assert(nBelow == nPrimitives && nAbove == 0);

        // Create leaf if no good splits were found
        if (bestAxis == -1 && retries < 2) {
            ++retries;
            axis = (axis + 1) % 3;
            goto retrySplit;
        }

        if (bestCost > oldCost)
            ++badRefines;

        if ((bestCost > 4 * oldCost && nPrimitives < 32) || bestAxis == -1 || badRefines == 3) {
            nodes[nodeNum].InitLeaf(primNums, nPrimitives, &primitiveIndices);
            return;
        }

        // Classify primitives with respect to split
        int n0 = 0, n1 = 0;
        for (int i = 0; i < bestOffset; ++i)
            if (edges[bestAxis][i].type == EdgeType::Start)
                prims0[n0++] = edges[bestAxis][i].primNum;
        for (int i = bestOffset + 1; i < 2 * nPrimitives; ++i)
            if (edges[bestAxis][i].type == EdgeType::End)
                prims1[n1++] = edges[bestAxis][i].primNum;

        // Recursively initialize children nodes
        float tSplit = edges[bestAxis][bestOffset].t;
        BBox bounds0 = nodeBounds, bounds1 = nodeBounds;
        bounds0.max[bestAxis] = bounds1.min[bestAxis] = tSplit;
        buildTree(nodeNum + 1, bounds0, allPrimBounds, prims0, n0, depth - 1, edges, prims0,
                  prims1 + nPrimitives, badRefines);

        int aboveChild = nextFreeNode;
        nodes[nodeNum].InitInterior(bestAxis, aboveChild, tSplit);
        buildTree(aboveChild, bounds1, allPrimBounds, prims1, n1, depth - 1, edges, prims0,
                  prims1 + nPrimitives, badRefines);
    }

public:
    void addPrimitive(Intersectable* prim) override { allPrimitives.push_back(prim); };

    void clear() override {
        delete[] nodes;
        nodes = nullptr;
        allPrimitives.clear();
    }

    bool isBuilt() const override { return nodes != nullptr; }

    bool intersect(const Ray& ray, float tMin, float tMax, Intersection& intersection) override {
        if (!bounds.testIntersect(ray)) {
            return false;
        }

        vec3 invDir(1 / ray.dir.x, 1 / ray.dir.y, 1 / ray.dir.z);
        constexpr int maxTodo = 64;
        KdToDo todo[maxTodo];
        int todoPos = 0;

        bool hit = false;
        const KDNode* node = &nodes[0];
        while (node != nullptr) {
            if (!node->IsLeaf()) {
                // Compute parametric distance along ray to split plane
                int axis = node->SplitAxis();
                float tPlane = (node->SplitPos() - ray.origin[axis]) * invDir[axis];

                // Get node children pointers for ray
                const KDNode *firstChild, *secondChild;
                int belowFirst = (ray.origin[axis] < node->SplitPos()) ||
                                 (ray.origin[axis] == node->SplitPos() && ray.dir[axis] <= 0);

                if (belowFirst) {
                    firstChild = node + 1;
                    secondChild = &nodes[node->AboveChild()];
                } else {
                    firstChild = &nodes[node->AboveChild()];
                    secondChild = node + 1;
                }

                // Advance to next child node, possibly enqueue other child
                if (tPlane > tMax || tPlane <= 0)
                    node = firstChild;
                else if (tPlane < tMin)
                    node = secondChild;
                else {
                    // Enqueue secondChild in todo list
                    todo[todoPos].node = secondChild;
                    todo[todoPos].tMin = tPlane;
                    todo[todoPos].tMax = tMax;
                    ++todoPos;
                    node = firstChild;
                    tMax = tPlane;
                }
            } else {
                int nPrimitives = node->nPrimitives();
                if (nPrimitives == 1) {
                    const auto p = allPrimitives[node->onePrimitive];
                    if (p->intersect(ray, tMin, tMax, intersection)) {
                        tMax = intersection.t;
                        hit = true;
                    }
                } else {
                    for (int i = 0; i < nPrimitives; ++i) {
                        int index = primitiveIndices[node->offsetIdx + i];
                        const auto p = allPrimitives[index];
                        if (p->intersect(ray, tMin, tMax, intersection)) {
                            tMax = intersection.t;
                            hit = true;
                        }
                    }
                }

                // Grab next node to process from todo list
                if (todoPos > 0) {
                    --todoPos;
                    node = todo[todoPos].node;
                    tMin = todo[todoPos].tMin;
                    tMax = todo[todoPos].tMax;
                } else
                    break;
            }
        }

        return hit;
    }

    ~KDTree() override { clear(); }

private:
    std::vector<Intersectable*> allPrimitives;
    std::vector<int> primitiveIndices;
    BBox bounds;
    KDNode* nodes = nullptr;
    int intersectCost = 1;
    int traversalCost = 1;
    int maxAllowedPrimitives = 1;
    int maxDepth = -1;
    int nAllocatedNodes = 0;
    int nextFreeNode = 0;
};

////////////////////////////////////////////////////////////////////////////////////////

struct BVHTree : IntersectionAccelerator {
    void addPrimitive(Intersectable* prim) override {}
    void clear() override {}
    void build(Purpose purpose) override {}
    bool isBuilt() const override { return false; }
    bool intersect(const Ray& ray, float tMin, float tMax, Intersection& intersection) override {
        return false;
    }
};

AcceleratorPtr makeDefaultAccelerator() {
    // TODO: uncomment or add the acceleration structure you have implemented
    int intersectCost = 80;
    int traversalCost = 1;

    return AcceleratorPtr(new KDTree(intersectCost, traversalCost));

    // return AcceleratorPtr(new KDTree());
    // return AcceleratorPtr(new BVHTree());
    // return AcceleratorPtr(new OctTree());
}
