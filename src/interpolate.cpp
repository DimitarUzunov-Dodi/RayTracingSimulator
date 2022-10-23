#include "interpolate.h"
#include <glm/geometric.hpp>

glm::vec3 computeBarycentricCoord (const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& p)
{
    // TODO: implement this function.
    float ABC = glm::length(glm::cross(v1 - v0, v2 - v0))/ 2;
    float APC = glm::length(glm::cross(p - v0, v2 - v0)) / 2; 
    float BPC = glm::length(glm::cross(v1 - p, v2 - p)) / 2; 


    float u = BPC/ABC;
    float v = APC/ABC;
    float w = (1 - BPC - APC)/ABC;

    glm::vec3 barycentricCoord = glm::vec3(u, v, w);

    return barycentricCoord;
}

glm::vec3 interpolateNormal (const glm::vec3& n0, const glm::vec3& n1, const glm::vec3& n2, const glm::vec3 barycentricCoord)
{
    glm::vec3 vNormal = n0 * barycentricCoord.x + n1 * barycentricCoord.y + n2 * barycentricCoord.z;

    return vNormal;
}

glm::vec2 interpolateTexCoord (const glm::vec2& t0, const glm::vec2& t1, const glm::vec2& t2, const glm::vec3 barycentricCoord)
{
// TODO: implement this function.
    return glm::vec2(0.0);
}
