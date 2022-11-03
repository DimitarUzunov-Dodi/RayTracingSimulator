#include "render.h"
#include "intersect.h"
#include "light.h"
#include "screen.h"
#include "texture.h"
#include <framework/trackball.h>
#ifdef NDEBUG
#include <omp.h>
#endif
#include <iostream>

#define MAX_RENDER_DEPTH 50 

glm::vec3 getFinalColor(const Scene& scene, const BvhInterface& bvh, Ray& ray, const Features& features, int rayDepth)
{
    HitInfo hitInfo;
    if (rayDepth <= MAX_RENDER_DEPTH && bvh.intersect(ray, hitInfo, features)) {

        glm::vec3 Lo = computeLightContribution(scene, bvh, features, ray, hitInfo);

        if (features.enableRecursive) {
            Ray reflection = computeReflectionRay(ray, hitInfo);
            reflection.origin += hitInfo.normal * std::numeric_limits<float>::epsilon();
            Lo += getFinalColor(scene, bvh, reflection, features, rayDepth + 1);
        }

        // Draw a white debug ray if the ray hits.
        drawRay(ray, Lo);

        //apply the texture to the textured objects
        if (features.enableTextureMapping && hitInfo.material.kdTexture) {

            return acquireTexel(*hitInfo.material.kdTexture, hitInfo.texCoord, features);   
        }
        // Set the color of the pixel to white if the ray hits.
        return Lo;
    } else {
        // Draw a red debug ray if the ray missed.
        drawRay(ray, glm::vec3(1.0f, 0.0f, 0.0f));
        // Set the color of the pixel to black if the ray misses.
        return glm::vec3(0.0f);
    }
}

void renderRayTracing(const Scene& scene, const Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features)
{
    glm::ivec2 windowResolution = screen.resolution();

    // Enable multi threading in Release mode
#ifdef NDEBUG
#pragma omp parallel for schedule(guided)
#endif
    for (int y = 0; y < windowResolution.y; y++) {
        for (int x = 0; x != windowResolution.x; x++) {
            // NOTE: (-1, -1) at the bottom left of the screen, (+1, +1) at the top right of the screen.
            const glm::vec2 normalizedPixelPos {
                float(x) / float(windowResolution.x) * 2.0f - 1.0f,
                float(y) / float(windowResolution.y) * 2.0f - 1.0f
            };


            Ray cameraRay = camera.generateRay(normalizedPixelPos);
            screen.setPixel(x, y, getFinalColor(scene, bvh, cameraRay, features));

            if (features.extra.enableMotionBlur && !screen.firstRender) {
                glm::ivec2 pixelPosition(x, y);
                glm::vec3 worldPosition = cameraRay.origin + cameraRay.direction * cameraRay.t;
                if (!screen.firstRender) {
                    screen.setVelocityBuffer(pixelPosition, windowResolution, worldPosition, cameraRay.t != std::numeric_limits<float>::max());
                }
            } 
        }
    }

    if (features.extra.enableMotionBlur) {
        if (!screen.firstRender) {
            screen.motionBlur(6);
        } else {
            screen.firstRender = false;
            screen.initVelocityBuffer(windowResolution.x * windowResolution.y);
        }
        screen.setPreviousCameraMatrix(camera);
    }
}