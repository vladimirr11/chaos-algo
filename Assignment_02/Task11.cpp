#include <algorithm>
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

    GracefulTreeValidator(GenericTree& _tree) : tree(_tree) {
        // const size_t doubledTreeSize = (tree.getSize() + 1) * 2;
        // for (size_t i = 1; i <= doubledTreeSize; ++i) {
        //     if (i & 1) {
        //         allowedNodeLabels.insert(i);
        //     }
        // }
        auto root = tree.getRoot();
        root->label = 1;
        for (size_t i = 1; i <= tree.getSize() + 1; i++) {
            allowedNodeLabels.insert(i);
        }
    }

    bool validate3() {
        std::deque<std::pair<sharedNodePtr, sharedNodePtr>> deq;
        sharedNodePtr root = tree.getRoot();
        root->label = 1;
        for (size_t c = 0; c < root->descendents.size(); c++) {
            deq.push_back(std::make_pair(root, tree.getNodeWithId(root->descendents[c]->id)));
        }

        if (_validate3(tree.getSize(), deq)) {
            return true;
        }
        return false;
     }

    bool _validate3(int edgeLabel, std::deque<std::pair<sharedNodePtr, sharedNodePtr>> deq) {
        if (edgeLabel == 0 || deq.empty()) {
            return true;
        }
        
        auto [parent, currNode] = deq.front();
        if (currNode->label == -1) {
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

            assert(currNodeLabel != 0 && "Failed labeling");
            currNode->label = currNodeLabel;
            if (!currNode->descendents.empty()) {
                for (size_t i = 0; i < currNode->descendents.size(); i++) {
                    deq.push_back(std::make_pair(currNode, tree.getNodeWithId(currNode->descendents[i]->id)));
                }
            }
            deq.pop_front();
            usedNodeLabels.insert(currNodeLabel);
            if (_validate3(edgeLabel -1, deq)) {
                return true;
            }
            usedNodeLabels.erase(currNodeLabel);
            deq.push_front(std::make_pair(parent, currNode));
        }

        return false;
    }

    bool validate2(int edgeLabel, int startIdx1, int startIdx2) {
        if (edgeLabel == 0) {
            return true;
        }

        for (size_t c = startIdx1; c < tree.getSize(); c++) {
            auto parent = tree.getNodeWithId(c);
            if (parent->label != -1) {
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

                assert(currNodeLabel != 0 && "Failed labeling");

                // int highNodeLabel = node->label + edgeLabel;
                for (size_t i = startIdx2; i < parent->descendents.size(); i++) {
                    if (allowedNodeLabels.find(currNodeLabel) != allowedNodeLabels.end() &&
                        usedNodeLabels.find(currNodeLabel) == usedNodeLabels.end()) {
                        parent->descendents[i]->label = currNodeLabel;
                        usedNodeLabels.insert(currNodeLabel);
                        if (validate2(edgeLabel - 1, startIdx1, startIdx2 + 1)) {
                            return true;
                        }
                        parent->descendents[i]->label = -1;
                        usedNodeLabels.erase(currNodeLabel);
                    }
                }
            }
        }

        return false;
    }

    void validate() {
        std::vector<bool> treeNodes(tree.getSize() + 1, false);
        sharedNodePtr root = tree.getRoot();
        root->label = 1;
        for (size_t c = 0; c < root->descendents.size(); c++) {
            treeNodes[root->descendents[c]->id] = true;
        }

        // const int startEdgeLabel = tree.getSize() * 2;
        const int startEdgeLabel = tree.getSize();
        if (_validate(startEdgeLabel, root, treeNodes)) {
            tree.print();
        } else {
            std::cout << "No valid labeling" << std::endl;
        }
    }

private:
    bool _validate(int edgeLabel, sharedNodePtr parent, std::vector<bool>& potentialNodes) {
        if (edgeLabel == 0) {
            return true;
        }

        for (size_t i = 0; i < potentialNodes.size(); i++) {
            if (potentialNodes[i] == true) {
                sharedNodePtr currNode = tree.getNodeWithId(i);
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

                assert(currNodeLabel != 0 && "Failed labeling");

                currNode->label = currNodeLabel;
                usedNodeLabels.insert(currNodeLabel);
                potentialNodes[i] = false;
                for (size_t j = 0; j < currNode->descendents.size(); j++) {
                    potentialNodes[currNode->descendents[j]->id] = true;
                }
                potentialNodes[i] = true;
                if (_validate(edgeLabel - 1, parent, potentialNodes)) {
                    return true;
                }
                for (size_t j = 0; j < currNode->descendents.size(); j++) {
                    potentialNodes[currNode->descendents[j]->id] = true;
                }
                potentialNodes[i] = true;
                usedNodeLabels.erase(currNodeLabel);
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

    // validator.validate();
    // if (validator.validate2(tree.getSize(), 0, 0)) {
    //     tree.print();
    // } else {
    //     tree.print();
    // }

    if (validator.validate3()) {
        tree.print();
    } else {
        tree.print();
    }

    return 0;
}
