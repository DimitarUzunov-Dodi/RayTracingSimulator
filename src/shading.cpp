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
        float dotKd = glm::dot(hitInfo.normal, lightDirection);

        if (dotKd > 0) {
            diffuse = lightColor * hitInfo.material.kd * dotKd;
        }

        
        glm::vec3 reflectVector = glm::normalize(-lightDirection + 2 * dotKd * glm::normalize(hitInfo.normal));
        glm::vec3 viewDirection = glm::normalize(-ray.direction);
        float dotKs = glm::dot(reflectVector, viewDirection);

        if (dotKd > 0 && dotKs > 0)
        {
            specular = lightColor * hitInfo.material.ks * pow(dotKs, hitInfo.material.shininess);
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