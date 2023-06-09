#include "light.h"
#include "config.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>
#include <iostream>
#include <random>

std::default_random_engine defEngine;
//Provides random double number in range (start,end) with uniform distribution
double getRandomNumInRange(const double &start, const double &end)
{
    std::uniform_real_distribution<double> dDistro(start, end);
    return dDistro(defEngine);
}

//Computes area of triangle plane
float area(glm::vec3 v0, glm::vec3 v1, glm::vec3 v2)
{
    return glm::length(glm::cross(v1 - v0, v2 - v0)) / 2.0;
}

//Computes average of two doubles
double avg(double x, double y) {
    return (x + y) / 2;
}

//Provides random point on line between "start" and "end"
glm::vec3 getRandomPoint(const glm::vec3& start, const glm::vec3& end) {
    return glm::vec3(getRandomNumInRange(start.x, end.x), getRandomNumInRange(start.y, end.y), getRandomNumInRange(start.z, end.z));
}

//Provides t for which normalized * t = finalVec
float getMaxTVal(glm::vec3 finalVec, glm::vec3 normalized) {
    float maxT = 0.0;

    if (normalized.x) {
        maxT = glm::max(maxT, (finalVec.x / normalized.x));
    }

    if (normalized.y) {
        maxT = glm::max(maxT, (finalVec.y / normalized.y)); 
    }

    if (normalized.z) {
        maxT = glm::max(maxT, (finalVec.z / normalized.z));
    }

    return maxT;
}

// samples a segment light source
// you should fill in the vectors position and color with the sampled position and color
void sampleSegmentLight(const SegmentLight& segmentLight, glm::vec3& position, glm::vec3& color)
{
    glm::vec3 direction = glm::normalize(segmentLight.endpoint1 - segmentLight.endpoint0);
    float maxT = getMaxTVal(segmentLight.endpoint1 - segmentLight.endpoint0, direction);
    float t = getRandomNumInRange(0.0f, maxT);

    position = segmentLight.endpoint0 + t * direction;
    glm::vec3 avgV = ((position - segmentLight.endpoint0) / (segmentLight.endpoint1 - segmentLight.endpoint0));
    double percentage = avg(avg(avgV.x, avgV.y), avgV.z);
    color = (float)(1.0 - percentage) * segmentLight.color0 + (float)percentage * segmentLight.color1;
}

// samples a parallelogram light source
// you should fill in the vectors position and color with the sampled position and color
void sampleParallelogramLight(const ParallelogramLight& parallelogramLight, glm::vec3& position, glm::vec3& color)
{   
    glm::vec3 edge01 = glm::normalize(parallelogramLight.edge01),
              edge02 = glm::normalize(parallelogramLight.edge02); //Normalized edge01/02

    float maxT1 = getMaxTVal(parallelogramLight.edge01, edge01),
          maxT2 = getMaxTVal(parallelogramLight.edge02, edge02); // t values such maxT1 * edge01 = parallelogramLight.edge01

    float t1 = getRandomNumInRange(0, maxT1), t2 = getRandomNumInRange(0, maxT2); //Random t_val
    glm::vec3 pointA = parallelogramLight.v0 + t1 * edge01,
              pointB = parallelogramLight.v0 + t2 * edge02; //Sampled points on the edges

    position = parallelogramLight.v0 + (t1 * edge01) + (t2 * edge02); //Sampled point  
  
    float parallelogramArea = area(parallelogramLight.v0, parallelogramLight.v0 + parallelogramLight.edge01, parallelogramLight.v0 + parallelogramLight.edge02) + 
        area(parallelogramLight.v0 + parallelogramLight.edge01, parallelogramLight.v0 + parallelogramLight.edge02, parallelogramLight.v0 + parallelogramLight.edge01 + parallelogramLight.edge02);
    //Area of whole parallelogram
    
    color = 
          (area(parallelogramLight.v0, pointA, pointB) 
            + area(pointA, pointB, position)
          ) / parallelogramArea * parallelogramLight.color0 

        + (area(pointA, parallelogramLight.v0 + parallelogramLight.edge01, position) 
            + area(position, parallelogramLight.v0 + parallelogramLight.edge01, pointB + parallelogramLight.edge01)
          ) / parallelogramArea * parallelogramLight.color1

        + (area(pointB, position, parallelogramLight.v0 + parallelogramLight.edge02) 
            + area(position, parallelogramLight.v0 + parallelogramLight.edge02, pointA + parallelogramLight.edge02)
          ) / parallelogramArea * parallelogramLight.color2

        + (area(position, pointB + parallelogramLight.edge01, pointA + parallelogramLight.edge02)
              + area(parallelogramLight.v0 + parallelogramLight.edge01 + parallelogramLight.edge02, pointB + parallelogramLight.edge01, pointA + parallelogramLight.edge02)
          ) / parallelogramArea * parallelogramLight.color3;

}

// test the visibility at a given light sample
// returns 1.0 if sample is visible, 0.0 otherwise
float testVisibilityLightSample(const glm::vec3& samplePos, const glm::vec3& debugColor, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo)
{
    Ray shadowRay;
    glm::vec3 point = (ray.origin + ray.direction * ray.t);
    shadowRay.direction = glm::normalize(point - samplePos);
    shadowRay.origin = samplePos;
    shadowRay.t = glm::distance(samplePos, point) - 0.00001;
    Ray debugRay = shadowRay;

    // TODO: implement this function.
    if (features.enableHardShadow) {
        if (bvh.intersect(shadowRay, hitInfo, features) || glm::dot(shadowRay.direction, hitInfo.normal) * glm::dot(ray.direction, hitInfo.normal) < 0)
        {
            drawRay(debugRay, glm::vec3(1,  0, 0));
            return 0.0;
        }
        drawRay(debugRay, debugColor);
        return 1.0;

    }

    return 1.0;
}

// given an intersection, computes the contribution from all light sources at the intersection point
// in this method you should cycle the light sources and for each one compute their contribution
// don't forget to check for visibility (shadows!)

// Lights are stored in a single array (scene.lights) where each item can be either a PointLight, SegmentLight or ParallelogramLight.
// You can check whether a light at index i is a PointLight using std::holds_alternative:
// std::holds_alternative<PointLight>(scene.lights[i])
//
// If it is indeed a point light, you can "convert" it to the correct type using std::get:
// PointLight pointLight = std::get<PointLight>(scene.lights[i]);
//
//
// The code to iterate over the lights thus looks like this:
// for (const auto& light : scene.lights) {
//     if (std::holds_alternative<PointLight>(light)) {
//         const PointLight pointLight = std::get<PointLight>(light);
//         // Perform your calculations for a point light.
//     } else if (std::holds_alternative<SegmentLight>(light)) {
//         const SegmentLight segmentLight = std::get<SegmentLight>(light);
//         // Perform your calculations for a segment light.
//     } else if (std::holds_alternative<ParallelogramLight>(light)) {
//         const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
//         // Perform your calculations for a parallelogram light.
//     }
// }
//
// Regarding the soft shadows for **other** light sources **extra** feature:
// To add a new light source, define your new light struct in scene.h and modify the Scene struct (also in scene.h)
// by adding your new custom light type to the lights std::variant. For example:
// std::vector<std::variant<PointLight, SegmentLight, ParallelogramLight, MyCustomLightType>> lights;
//
// You can add the light sources programmatically by creating a custom scene (modify the Custom case in the
// loadScene function in scene.cpp). Custom lights will not be visible in rasterization view.

#define NUM_OF_SAMPLES 100
glm::vec3 computeLightContribution(const Scene& scene, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo)
{
    if (features.enableShading) {
        glm::vec3 result = glm::vec3(0, 0, 0);
        glm::vec3 lightContribution;
         //If shading is enabled, compute the contribution from all lights.
         for (const auto& light : scene.lights) {
             if (std::holds_alternative<PointLight>(light)) {

                const PointLight pointLight = std::get<PointLight>(light);    
                lightContribution = computeShading(pointLight.position, pointLight.color, features, ray, hitInfo);
                lightContribution *= testVisibilityLightSample(pointLight.position, lightContribution, bvh, features, ray, hitInfo);

             } else if (std::holds_alternative<SegmentLight>(light)) {

                 const SegmentLight segmentLight = std::get<SegmentLight>(light);
                 if (!features.enableAreaLightSampling) {
                     lightContribution = 0.5f * computeShading(segmentLight.endpoint0, segmentLight.color0, features, ray, hitInfo) 
                         + 0.5f * computeShading(segmentLight.endpoint1, segmentLight.color1, features, ray, hitInfo);
                 } else {
                     lightContribution = glm::vec3(0.0f);
                     for (int i = 0; i < NUM_OF_SAMPLES; i++) {
                         glm::vec3 sampledPosition, sampledColor;
                         sampleSegmentLight(segmentLight, sampledPosition, sampledColor);
                         
                         lightContribution += computeShading(sampledPosition, sampledColor, features, ray, hitInfo);
                     }

                     lightContribution /= (float)NUM_OF_SAMPLES;
                 }
             } else if (std::holds_alternative<ParallelogramLight>(light)) {
                 const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
                 if (!features.enableAreaLightSampling) {
                     lightContribution = 0.25f * computeShading(parallelogramLight.v0, parallelogramLight.color0, features, ray, hitInfo) 
                         + 0.25f * computeShading(parallelogramLight.v0 + parallelogramLight.edge01, parallelogramLight.color1, features, ray, hitInfo) 
                         + 0.25f * computeShading(parallelogramLight.v0 + parallelogramLight.edge02, parallelogramLight.color2, features, ray, hitInfo) 
                         + 0.25f * computeShading(parallelogramLight.v0 + parallelogramLight.edge01 + parallelogramLight.edge02, parallelogramLight.color0, features, ray, hitInfo);
                 } else {
                     lightContribution = glm::vec3(0.0f);
                     for (int i = 0; i < NUM_OF_SAMPLES; i++) {
                         glm::vec3 sampledPosition, sampledColor;
                         sampleParallelogramLight(parallelogramLight, sampledPosition, sampledColor);

                         lightContribution += computeShading(sampledPosition, sampledColor, features, ray, hitInfo);
                     }

                     lightContribution /= (float)NUM_OF_SAMPLES;
                 }
             }
             result += lightContribution;
             
         }
        // TODO: replace this by your own implementation of shading
         
        return result;

    } else {
        // If shading is disabled, return the albedo of the material.
        return hitInfo.material.kd;
    }
}
