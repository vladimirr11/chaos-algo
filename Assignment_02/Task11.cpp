#include <cassert>
#include <deque>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_set>
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

    std::shared_ptr<TreeNode> getNodeWithId(int targetNodeId) const {
        assert((size_t)targetNodeId <= treeSize && "Target node greater than tree size");
        return _getNode(root, targetNodeId);
    }

    std::shared_ptr<TreeNode> getRoot() { return root; }

    size_t getSize() const { return treeSize; }

    void print() const { _printPreorder(root); }

private:
    void _insertNode(std::shared_ptr<TreeNode> node, int parentId, int currNodeId) {
        if (node->id == parentId) {
            auto newChild = std::make_shared<TreeNode>(-1, currNodeId);
            node->descendents.push_back(newChild);
            return;
        }

        for (const std::shared_ptr<TreeNode>& child : node->descendents) {
            _insertNode(child, parentId, currNodeId);
        }
    }

    std::shared_ptr<TreeNode> _getNode(std::shared_ptr<TreeNode> node, int targetNodeId) const {
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

    void _printPreorder(std::shared_ptr<TreeNode> node) const {
        if ((size_t)node->id == treeSize) {
            return;
        }

        for (const std::shared_ptr<TreeNode>& child : node->descendents) {
            const int absDiff = std::abs(node->label - child->label);
            std::cout << "|" << node->label << " - " << child->label << "| = " << absDiff
                      << std::endl;
            _printPreorder(child);
        }
    }

private:
    std::shared_ptr<TreeNode> root;
    size_t treeSize{};
};

class GracefulTreeValidator {
public:
    using sharedNodePtr = std::shared_ptr<GenericTree::TreeNode>;
    using parentChildDeque = std::deque<std::pair<sharedNodePtr, sharedNodePtr>>;

    GracefulTreeValidator(GenericTree& _tree) : tree(_tree) {
        const size_t doubledTreeSize = (tree.getSize() + 1) * 2;
        for (size_t i = 1; i <= doubledTreeSize; ++i) {
            if (i & 1) {
                allowedNodeLabels.insert(i);
            }
        }
    }

    bool isValidTree() {
        parentChildDeque treeNodesDeq;
        sharedNodePtr root = tree.getRoot();
        root->label = 1;
        for (sharedNodePtr& child : root->descendents) {
            treeNodesDeq.push_back(std::make_pair(root, child));
        }

        const size_t startEdgeLabel = tree.getSize() * 2;
        if (_validate(startEdgeLabel, treeNodesDeq)) {
            return true;
        }
        return false;
     }

private:
    bool _validate(int edgeLabel, parentChildDeque treeNodesDeq) {
        if (edgeLabel == 0 || treeNodesDeq.empty()) {
            return true;
        }
        
        const auto [parent, child] = treeNodesDeq.front();
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
                for (sharedNodePtr& childChild : child->descendents) {
                    treeNodesDeq.push_back(std::make_pair(child, childChild));
                }
            }
            treeNodesDeq.pop_front();
            usedNodeLabels.insert(currNodeLabel);
            if (_validate(edgeLabel - 2, treeNodesDeq)) {
                return true;
            }
            usedNodeLabels.erase(currNodeLabel);
            treeNodesDeq.push_front(std::make_pair(parent, child));
        }

        return false;
    }

private:
    GenericTree tree;
    std::unordered_set<int> allowedNodeLabels;
    std::unordered_set<int> usedNodeLabels;
};

int main([[maybe_unused]] int argc, [[maybe_unused]] char* argv[]) {
    if (argc == 1) {
        puts("No input file passed");
        // exit(EXIT_FAILURE);
    }

    // const char* fileName = argv[1];
    const char* fileName = "Data.txt";
    std::ifstream inputFile(fileName, std::ios::in);
    if (!inputFile.is_open()) {
        printf("Failed to open %s\n", fileName);
        exit(EXIT_FAILURE);
    }

    GenericTree tree;
    tree.buildTree(inputFile);

    GracefulTreeValidator validator(tree);
    if (validator.isValidTree()) {
        tree.print();
    } else {
        puts("No graceful labeling found");
    }

    return 0;
}
