#include "interpolate.h"
#include <glm/geometric.hpp>

glm::vec3 computeBarycentricCoord (const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& p)
{
    // TODO: implement this function.
    float ABC = glm::length(glm::cross(v1 - v0, v2 - v0))/ 2;
    float APC = glm::length(glm::cross(p - v0, v2 - v0))/ 2; 
    float BPC = glm::length(glm::cross(v1 - p, v2 - p))/ 2; 


    float u = BPC/ABC;
    float v = APC/ABC;
    float w = 1 - u - v;

    glm::vec3 barycentricCoord = glm::vec3(u, v, w);

    return barycentricCoord;
}

glm::vec3 interpolateNormal (const glm::vec3& n0, const glm::vec3& n1, const glm::vec3& n2, const glm::vec3 barycentricCoord)
{
    glm::vec3 a = barycentricCoord.x * n0;
    glm::vec3 b = barycentricCoord.y * n1;
    glm::vec3 c = barycentricCoord.z * n2; 
    
    glm::vec3 intrNormal = a + b + c;

    return intrNormal;
}

glm::vec2 interpolateTexCoord (const glm::vec2& t0, const glm::vec2& t1, const glm::vec2& t2, const glm::vec3 barycentricCoord)
{
    glm::vec2 a = barycentricCoord.x * t0;
    glm::vec2 b = barycentricCoord.y * t1;
    glm::vec2 c = barycentricCoord.z * t2;

    glm::vec2 intrCoord = a + b + c;

    return intrCoord;
}
