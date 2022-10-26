#include "config.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>
#include <vector>

bool validatePosition(std::vector<std::vector<glm::vec3>>& img, int a, int b)
{
    return a >= 0 && a < img.size() && b >= 0 && b < img[0].size();
}

std::vector<std::pair<int, int>> moves {
    { -1, -1 }, { -1, 0 }, {-1, 1},
    { 0, -1 }, { 0, 0 }, { 0, 1 },
    { 1, -1 }, { 1, 0 }, {1, 1}
};

//Only takes colors above given threshold
std::vector<std::vector<glm::vec3>> getThresholdedImage(std::vector<std::vector<glm::vec3>>& img, float threshold)
{
    std::vector<std::vector<glm::vec3>> retObj;
    for (int i = 0; i < img.size(); i++) {
        retObj.push_back(std::vector<glm::vec3> {});
        for (int j = 0; j < img[i].size(); j++) {
            if (img[i][j].x < threshold || img[i][j].y < threshold || img[i][j].z < threshold)
                retObj[i].push_back(glm::vec3 { 0.0 });
            else
                retObj[i].push_back(img[i][j]);
        }
    }

    return retObj;
}

//Performs boxFilter operation on given image
std::vector<std::vector<glm::vec3>> boxFilter(std::vector<std::vector<glm::vec3>>& img)
{
    std::vector<std::vector<glm::vec3>> retObj;
    for (int i = 0; i < img.size(); i++) {
        retObj.push_back(std::vector<glm::vec3> {});
        for (int j = 0; j < img[i].size(); j++) {
            double sumX = 0.0, sumY = 0.0, sumZ = 0.0, count = 0.0;

            for (auto move : moves) {
                if (!validatePosition(img, i + move.first, j + move.second))
                    continue;

                count++;
                sumX += img[i + move.first][j + move.second].x;
                sumY += img[i + move.first][j + move.second].y;
                sumZ += img[i + move.first][j + move.second].z;
            }
            retObj[i].push_back(glm::vec3(sumX / count, sumY / count, sumZ / count));
        }
    }

    return retObj;
}