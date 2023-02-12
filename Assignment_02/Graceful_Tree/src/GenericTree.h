#ifndef GENERICTREE_H
#define GENERICTREE_H

#include <cassert>
#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

class GenericTree {
public:
    struct TreeNode {
        int label = -1;
        int id = -1;
        std::vector<std::shared_ptr<TreeNode>> descendents;

        TreeNode() = default;

        TreeNode(const int _label, const int _id) : label(_label), id(_id) {}
    };

public:
    GenericTree() = default;

    void buildTree(std::ifstream& inputFile) {
        assert(inputFile.good() && "Precondition for buildTree() not met - input file not good");

        std::string line;
        while (std::getline(inputFile >> std::ws, line) && line != "---") {
            if (!root) {
                root = std::make_shared<TreeNode>(-1, 0);
                continue;
            }

            treeSize++;
            const int parentId = int(line[0] - '0');
            _insertNode(root, parentId, treeSize);
        }
    }

    std::shared_ptr<TreeNode> getNodeWithId(const int targetNodeId) const {
        assert((size_t)targetNodeId <= treeSize && "Requested node id greater than tree size");
        return _getNode(root, targetNodeId);
    }

    std::shared_ptr<TreeNode> getRoot() const { return root; }

    size_t getSize() const { return treeSize; }

    void print(std::ostream& outputStream) const { _printPreorder(outputStream, root); }

private:
    void _insertNode(std::shared_ptr<TreeNode> node, const int parentId, const int currNodeId) {
        if (node->id == parentId) {
            const auto newChild = std::make_shared<TreeNode>(-1, currNodeId);
            node->descendents.push_back(newChild);
            return;
        }

        for (const std::shared_ptr<TreeNode>& child : node->descendents) {
            _insertNode(child, parentId, currNodeId);
        }
    }

    std::shared_ptr<TreeNode> _getNode(std::shared_ptr<TreeNode> node, const int targetNodeId) const {
        if (node->id == targetNodeId) {
            return node;
        }

        std::shared_ptr<TreeNode> targetNode;
        for (const std::shared_ptr<TreeNode>& child : node->descendents) {
            targetNode = _getNode(child, targetNodeId);
            if (targetNode && targetNode->id == targetNodeId) {
                return targetNode;
            }
        }

        return targetNode;
    }

    void _printPreorder(std::ostream& outputStream, std::shared_ptr<TreeNode> node) const {
        if ((size_t)node->id == treeSize) {
            return;
        }

        for (const std::shared_ptr<TreeNode>& child : node->descendents) {
            const int absDiff = std::abs(node->label - child->label);
            outputStream << "|" << node->label << " - " << child->label << "| = " << absDiff
                      << std::endl;
            _printPreorder(outputStream, child);
        }
    }

private:
    std::shared_ptr<TreeNode> root;
    size_t treeSize{};
};

#endif // !GENERICTREE_H
