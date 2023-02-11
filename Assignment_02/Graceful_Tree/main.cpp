#include "src/GracefulTreeValidator.h"
#include "src/TreeGenerator.h"

int main([[maybe_unused]] int argc, [[maybe_unused]] char* argv[]) {
    std::ifstream inputFile;

    if (argc == 1) {
        puts("No input file provided. Start generating random tree structure...");

        const char* fileName = "TreeStructure.txt";
        std::ofstream treeStructureFile(fileName, std::ios::out | std::ios::trunc);
        if (!treeStructureFile.is_open()) {
            printf("Failed to open %s\n", fileName);
            exit(EXIT_FAILURE);
        }

        TreeGenerator generator(15); // 15 vertices tree structure
        generator.generateTreeStructure(treeStructureFile);
        assert(!treeStructureFile.bad() && "Bad bit set while generating tree structure");

        treeStructureFile.close();

        inputFile.open(fileName, std::ios::in);
        if (!inputFile.is_open()) {
            printf("Failed to open %s\n", fileName);
            exit(EXIT_FAILURE);
        }

    } else {
        const char* fileName = argv[1];
        inputFile.open(fileName, std::ios::in);
        if (!inputFile.is_open()) {
            printf("Failed to open %s\n", fileName);
            exit(EXIT_FAILURE);
        }
    }

    const char* outputFileName = "TreeLabeling.txt";
    std::ofstream outputFile(outputFileName, std::ios::out | std::ios::trunc);
    if (!outputFile.good()) {
        printf("Cannot create output file %s", outputFileName);
        exit(EXIT_FAILURE);
    }

    GenericTree tree;
    tree.buildTree(inputFile);

    GracefulTreeValidator validator(tree);
    if (validator.isValidTree()) {
        validator.recordTree(outputFile);
    } else {
        puts("No graceful labeling found");
    }

    outputFile.close();

    return 0;
}
