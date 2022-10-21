#include "texture.h"
#include <cmath>
#include <glm/geometric.hpp>
#include <shading.h>

const glm::vec3 computeShading(const glm::vec3& lightPosition, const glm::vec3& lightColor, const Features& features, Ray ray, HitInfo hitInfo)
{

        glm::vec3 diffuse = glm::vec3(0, 0, 0);
        glm::vec3 specular = glm::vec3(0, 0, 0);

        glm::vec3 vertex = ray.origin + ray.direction * ray.t;
        glm::vec3 lightDirection = glm::normalize(lightPosition - vertex);
        float dot = glm::dot(hitInfo.normal, lightDirection);

        if (dot > 0) {
            diffuse = lightColor * hitInfo.material.kd * dot;
        }


        glm::vec3 reflectVector = glm::normalize(-lightDirection + 2 * glm::dot(lightDirection, hitInfo.normal) * glm::normalize(hitInfo.normal));
        glm::vec3 viewDirection = glm::normalize(-ray.direction);
        

        if (glm::dot(lightDirection, hitInfo.normal) > 0 && glm::dot(reflectVector, viewDirection) > 0)
        {
            specular = lightColor * hitInfo.material.ks * pow(glm::dot(reflectVector, viewDirection), hitInfo.material.shininess);
        }
        return  diffuse + specular;
    
}


const Ray computeReflectionRay (Ray ray, HitInfo hitInfo)
{
    // Do NOT use glm::reflect!! write your own code.
    Ray reflectionRay {};
    // TODO: implement the reflection ray computation.
    return reflectionRay;
}