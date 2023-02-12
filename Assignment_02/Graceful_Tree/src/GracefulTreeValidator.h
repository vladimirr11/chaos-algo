#ifndef GRACEFULTREEVALIDATOR_H
#define GRACEFULTREEVALIDATOR_H

#include <deque>
#include <unordered_set>
#include "GenericTree.h"

struct TreeEdge {
    using GenericTreeNodePtr = std::shared_ptr<GenericTree::TreeNode>;

    GenericTreeNodePtr parent = nullptr;
    GenericTreeNodePtr child = nullptr;

    TreeEdge() = default;

    TreeEdge(const GenericTreeNodePtr _parent, const GenericTreeNodePtr _child)
        : parent(_parent), child(_child) {}
};

class GracefulTreeValidator {
public:
    using SharedNodePtr = TreeEdge::GenericTreeNodePtr;
    using TreeEdgesDeque = std::deque<TreeEdge>;

    GracefulTreeValidator() = delete;

    GracefulTreeValidator(const GenericTree& _tree) : tree(_tree) {
        const size_t doubledTreeSize = (tree.getSize() + 1) << 1;
        for (size_t i = 1; i <= doubledTreeSize; i++) {
            if (i & 1) {
                allowedNodeLabels.insert(i);
            }
        }
    }

    bool isValidTree() {
        TreeEdgesDeque treeEdgesDeq;
        SharedNodePtr root = tree.getRoot();
        root->label = 1;
        for (SharedNodePtr& child : root->descendents) {
            treeEdgesDeq.emplace_back(TreeEdge(root, child));
        }

        const size_t startEdgeLabel = tree.getSize() << 1;
        if (_validate(startEdgeLabel, treeEdgesDeq)) {
            return true;
        }
        return false;
    }

    void recordTree(std::ostream& outputStream) const { tree.print(outputStream); }

private:
    // Deterministic backtracking search
    // References:
    // [1] W. Fang. A computational approach to the graceful tree conjecture
    // [2] M. Horton. Graceful Trees: Statistics and Algorithms
    bool _validate(int edgeLabel, TreeEdgesDeque& treeEdgesDeq) {
        if (edgeLabel == 0 || treeEdgesDeq.empty()) {
            return true;
        }

        const auto [parent, child] = treeEdgesDeq.front();
        if (child->label == -1) {
            const int lowLabel = parent->label - edgeLabel;
            const int highLabel = parent->label + edgeLabel;
            int currNodeLabel{};
            if (allowedNodeLabels.find(lowLabel) != allowedNodeLabels.end() &&
                usedNodeLabels.find(currNodeLabel) == usedNodeLabels.end()) {
                currNodeLabel = lowLabel;
            } else if (allowedNodeLabels.find(highLabel) != allowedNodeLabels.end() &&
                       usedNodeLabels.find(currNodeLabel) == usedNodeLabels.end()) {
                currNodeLabel = highLabel;
            }

            assert(currNodeLabel != 0 && "Failed extracting valid label");

            child->label = currNodeLabel;
            if (!child->descendents.empty()) {
                for (SharedNodePtr& childChild : child->descendents) {
                    treeEdgesDeq.emplace_back(TreeEdge(child, childChild));
                }
            }
            treeEdgesDeq.pop_front();
            usedNodeLabels.insert(currNodeLabel);

            if (_validate(edgeLabel - 2, treeEdgesDeq)) {
                return true;
            }
            
            treeEdgesDeq.emplace_front(TreeEdge(parent, child));
            usedNodeLabels.erase(currNodeLabel);
        }

        return false;
    }

private:
    GenericTree tree;
    std::unordered_set<int> allowedNodeLabels;
    std::unordered_set<int> usedNodeLabels;
};

#endif  // !GRACEFULTREEVALIDATOR_H
