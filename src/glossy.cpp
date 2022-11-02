#include "texture.h"
#include <cmath>
#include <glm/geometric.hpp>
#include <shading.h>
#include <glm/glm.hpp>
#include <queue>

#define GLOSSY_RAY_AVG 50
#define SHININES_FACTOR 30.0f


// Computes a random ray in the square around our reflect Ray.
const Ray computePerturbedRay(Ray ray, HitInfo hitInfo, const Features& features ) {
    float a = hitInfo.material.shininess / SHININES_FACTOR ;
    float v = -a / 2.0f + ((float)rand() / RAND_MAX) * a;
    float u = -a / 2.0f + ((float)rand() / RAND_MAX) * a;

    glm::vec3 w = glm::normalize(ray.direction);
    glm::vec3 t = w;
    //we take a random non-colinear and then cross it with our direction 
    

    if (std::abs(w.x) > std::abs(w.y)) {
        if (std::abs(w.y) > std::abs(w.z)) {
            t.z = 1.0f;
        } else {
            t.y = 1.0f;
        }

    } else if (std::abs(w.x) > std::abs(w.z)) {
        t.z = 1.0f;
    } else {
        t.z = 1.0f;
    }
    glm::vec3 x = glm::normalize(glm::cross(w, glm::normalize(t)));    
    glm::vec3 y = glm::cross(x,w);



    Ray preturbedRay
    {
            ray.origin,
            ray.direction + (u * x) + (v * y)

    };
    if (glm::dot(glm::normalize(preturbedRay.direction), hitInfo.normal) > 0) {
        return preturbedRay;
    }
    return Ray {ray.origin,glm::vec3(0,0,0)};

}

const Ray computeGlossyRay(Ray ray, HitInfo hitInfo, const Features& features)
{
    glm::vec3 sumDirections = ray.direction;
    for (int i = 0; i < GLOSSY_RAY_AVG; i++) {
        sumDirections += computePerturbedRay(ray, hitInfo, features).direction;
    }

    glm::vec3 avgDir = sumDirections / float(GLOSSY_RAY_AVG+1);
    
    return Ray { ray.origin, avgDir };
}