#ifndef TREEGENERATOR_H
#define TREEGENERATOR_H

#include <iostream>
#include <random>
#include <string>

struct TreeGenerator {
    int numVertices = 0;

    TreeGenerator() = delete;

    TreeGenerator(const int _numVertices) : numVertices(_numVertices) {}

    void generateTreeStructure(std::ostream& outputStream) {
        std::mt19937 gen(42);
        outputStream << "-1" << std::endl;
        for (int c = 0; c < numVertices; c++) {
            std::uniform_int_distribution<uint8_t> dist(0, c);
            outputStream <<  std::to_string(dist(gen)) << std::endl;
        }
        outputStream << "---" << std::endl;
    }
};

#endif  // !TREEGENERATOR_H
