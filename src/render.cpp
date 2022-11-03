#include "render.h"
#include "intersect.h"
#include "light.h"
#include "screen.h"
#include "texture.h"
#include <framework/trackball.h>
#ifdef NDEBUG
#include <omp.h>
#endif

#define MAX_RENDER_DEPTH 5 

glm::vec3 getFinalColor(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth)
{
    HitInfo hitInfo, hitInfo2;
    if (rayDepth <= MAX_RENDER_DEPTH && bvh.intersect(ray, hitInfo, features)) {
        glm::vec3 Lo = computeLightContribution(scene, bvh, features, ray, hitInfo);

        if (features.enableRecursive) {
            Ray reflection = computeReflectionRay(ray, hitInfo);
            reflection.origin += hitInfo.normal * std::numeric_limits<float>::epsilon();
            glm::vec3 recursiveCallResult = getFinalColor(scene, bvh, reflection, features, rayDepth + 1);


            Lo += recursiveCallResult;  
        }

        if (features.extra.enableTransparency) {
            Ray transparentRay { ray.origin + ray.direction * ray.t, glm::normalize(ray.direction) };
            transparentRay.origin += glm::normalize(ray.direction) * std::numeric_limits<float>::epsilon();
            glm::vec3 transparentColor = glm::vec3(0.0f);

            if (features.enableRecursive) {
                transparentColor = getFinalColor(scene, bvh, transparentRay, features, rayDepth + 1);
                drawRay(transparentRay, transparentColor); // VISUAL DEBUG
            } else if (bvh.intersect(transparentRay, hitInfo2, features)) {
                transparentColor = computeLightContribution(scene, bvh, features, transparentRay, hitInfo2);
                drawRay(transparentRay, transparentColor); // VISUAL DEBUG
            }

            Lo = hitInfo.material.transparency * Lo + (1 - hitInfo.material.transparency) * transparentColor;
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
            const Ray cameraRay = camera.generateRay(normalizedPixelPos);
            screen.setPixel(x, y, getFinalColor(scene, bvh, cameraRay, features));
        }
    }
}