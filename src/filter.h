#pragma once
#include "config.h"

std::vector<std::vector<glm::vec3>> getThresholdedImage(std::vector<std::vector<glm::vec3>>& img, float threshold);

std::vector<std::vector<glm::vec3>> boxFilter(std::vector<std::vector<glm::vec3>>& img, const int& boxSize = 3);