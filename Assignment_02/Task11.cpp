#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <memory>
#include <unordered_set>
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

    std::shared_ptr<TreeNode> getNodeWithId(int targetNodeId) const {
        assert((size_t)targetNodeId <= treeSize && "Target node greater than tree size");
        return _getNode(root, targetNodeId);
    }

    std::shared_ptr<TreeNode> getRoot() { return root; }

    size_t getSize() const { return treeSize; }

    void print() const { _print(root); }

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

    void _print(std::shared_ptr<TreeNode> node) const {
        if ((size_t)node->id == treeSize) {
            return;
        }

        for (const std::shared_ptr<TreeNode>& child : node->descendents) {
            const int absDiff = std::abs(node->label - child->label);
            std::cout << "|" << node->label << " - " << child->label << "| = " << absDiff
                      << std::endl;
            _print(child);
        }
    }

private:
    std::shared_ptr<TreeNode> root;
    size_t treeSize{};
};

class GracefulTreeValidator {
public:
    using sharedNodePtr = std::shared_ptr<GenericTree::TreeNode>;

    GracefulTreeValidator(GenericTree& _tree) : tree(_tree) {
        for (size_t i = 1; i <= (tree.getSize() + 1) * 2; i += 2) {
            allowedNodeLabels.insert(i);
        }
    }

    void validate() {
        std::vector<bool> possibleNodes(tree.getSize() + 1, false);
        auto root = tree.getRoot();
        root->label = 1;
        for (size_t c = 0; c < root->descendents.size(); c++) {
            possibleNodes[root->descendents[c]->id] = true;
        }
        const int startEdgeLabel = tree.getSize() * 2;
        if (_validate(startEdgeLabel, root, possibleNodes)) {
            tree.print();
        } else {
            std::cout << "No valid labeling" << std::endl;
        }
    }

private:
    bool _validate(int edgeLabel, sharedNodePtr parent, std::vector<bool> possibleNodes) {
        if (edgeLabel == 0) {
            return true;
        }

        for (size_t i = 0; i < possibleNodes.size(); i++) {
            if (possibleNodes[i] == true) {
                auto possibleNode = tree.getNodeWithId(i);
                int lowLabel = parent->label - edgeLabel;
                int highLabel = parent->label + edgeLabel;
                int potLabel{};
                if (allowedNodeLabels.find(lowLabel) != allowedNodeLabels.end()) {
                    potLabel = lowLabel;
                } else if (allowedNodeLabels.find(highLabel) != allowedNodeLabels.end()) {
                    potLabel = highLabel;
                } else {
                    std::cerr << "NO LABEL" << std::endl;
                    exit(1);
                }

                if (usedNodeLabels.find(potLabel) == usedNodeLabels.end()) {
                    possibleNode->label = potLabel;
                    usedNodeLabels.insert(potLabel);
                    possibleNodes[i] = false;
                    if (_validate(edgeLabel - 2, parent, possibleNodes)) {
                        return true;
                    }
                    for (size_t j = 0; j < possibleNode->descendents.size(); j++) {
                        possibleNodes[possibleNode->descendents[j]->id] = true;
                    }
                    possibleNodes[i] = true;
                    usedNodeLabels.erase(potLabel);
                }
            }
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

    validator.validate();

    return 0;
}
