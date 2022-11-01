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

//Performs boxFilter operation of dimensions boxSize x boxSize on given image
std::vector<std::vector<glm::vec3>> boxFilter(std::vector<std::vector<glm::vec3>>& img, const int& boxSize)
{
    std::vector<std::vector<glm::vec3>> retObj;
    for (int i = 0; i < img.size(); i++) {
        retObj.push_back(std::vector<glm::vec3> {});
        for (int j = 0; j < img[i].size(); j++) {
            double sumX = 0.0, sumY = 0.0, sumZ = 0.0, count = 0.0;

            for (int k = i - boxSize; k <= i + boxSize; k++) {
                for (int l = j - boxSize; l <= j + boxSize; l++) {
                    if (!validatePosition(img, k, l))
                        continue;


                    count++;
                    sumX += img[k][l].x;
                    sumY += img[k][l].y;
                    sumZ += img[k][l].z;
                }
            }

            retObj[i].push_back(glm::vec3(sumX / count, sumY / count, sumZ / count));
        }
    }

    return retObj;
}