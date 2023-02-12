#include "src/GracefulTreeValidator.h"
#include "src/TreeGenerator.h"

int main(int argc, char* argv[]) {
    std::string inputFileName;
    std::ifstream inputFile;

    if (argc == 1) {
        puts("No input file provided. Start generating random tree structure...");

        inputFileName = "TreeStructure.txt";
        std::ofstream treeStructureFile(inputFileName, std::ios::out | std::ios::trunc);
        if (!treeStructureFile.good()) {
            printf("Failed to open %s\n", inputFileName.c_str());
            exit(EXIT_FAILURE);
        }

        const int numVertices = 25;
        TreeGenerator generator(numVertices);
        generator.generateTreeStructure(treeStructureFile);
        assert(!treeStructureFile.bad() && "Bad bit set while generating tree structure");

        treeStructureFile.close();

    } else {
        inputFileName = argv[1];
    }

    inputFile.open(inputFileName, std::ios::in);
    if (!inputFile.is_open()) {
        printf("Failed to open %s\n", inputFileName.c_str());
        exit(EXIT_FAILURE);
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
