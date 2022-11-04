#include "intersect.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/vector_relational.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>
#include <limits>

float area(glm::vec3& v0, glm::vec& v1, glm::vec& v2)
{
    return glm::length(glm::cross(v1 - v0, v2 - v0)) / 2.0;
}

bool pointInTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& n, const glm::vec3& p) {
    return glm::dot(glm::cross(v1 - v0, p - v0), n) >= 0
        && glm::dot(glm::cross(v2 - v1, p - v1), n) >= 0
        && glm::dot(glm::cross(v0 - v2, p - v2), n) >= 0;
}

bool intersectRayWithPlane(const Plane& plane, Ray& ray)
{
    float num = plane.D - glm::dot(plane.normal, ray.origin),
          det = glm::dot(plane.normal, ray.direction);
    if (!det)
        return false; // orthogonal

    if (num / det >= ray.t || num / det < 0)
        return false;

    ray.t = num / det;
    return true;
}

Plane trianglePlane(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2)
{
    Plane plane;
    plane.normal = glm::normalize(glm::cross(v1 - v0, v2 - v0));
    plane.D = glm::dot(plane.normal, v0);
    return plane;
}

/// Input: the three vertices of the triangle
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, Ray& ray, HitInfo& hitInfo)
{
    Plane plane = trianglePlane(v0, v1, v2);
    float prevT = ray.t;

    if (intersectRayWithPlane(plane, ray)) {
        if (pointInTriangle(v0, v1, v2, plane.normal, ray.origin + ray.direction * ray.t)) {
            float a = area(ray.origin + ray.direction * ray.t, v1, v2) / area(v0, v1, v2);
            float b = area(ray.origin + ray.direction * ray.t, v0, v2) / area(v0, v1, v2);
            float c = area(ray.origin + ray.direction * ray.t, v0, v1) / area(v0, v1, v2);

            hitInfo.normal = glm::normalize(a * n1 + b * n2 + c * n3); //Interpolates normal
            if (glm::dot(plane.normal, -raydirection) <= 0)
                hitInfo.normal *= -1;

            return true;
        }
        ray.t = prevT;
    }

    return false;
}

/// Input: a sphere with the following attributes: sphere.radius, sphere.center
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const Sphere& sphere, Ray& ray, HitInfo& hitInfo)
{
    float a = glm::length(ray.direction) * glm::length(ray.direction),
          b = 2 * glm::dot(ray.direction, ray.origin - sphere.center),
          c = glm::length(ray.origin - sphere.center) * glm::length(ray.origin - sphere.center) - sphere.radius * sphere.radius;

    float delta = b * b - 4 * a * c;
    if (delta < 0)
        return false;

    float s1 = (-b - (float)std::sqrt(delta)) / (2 * a),
          s2 = (-b + (float)std::sqrt(delta)) / (2 * a);

    if (s1 >= 0 && s1 < ray.t) {
        ray.t = s1;
        hitInfo.normal = glm::normalize(ray.origin + ray.direction * ray.t - sphere.center);
        return true;
    } else if (s2 >= 0 && s2 < ray.t) {
        ray.t = s2;
        hitInfo.normal = glm::normalize(ray.origin + ray.direction * ray.t - sphere.center);
        return true;
    }

    return false;
}

/// Input: an axis-aligned bounding box with the following parameters: minimum coordinates box.lower and maximum coordinates box.upper
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const AxisAlignedBox& box, Ray& ray)
{
    glm::vec3 tLower = (box.lower - ray.origin) / ray.direction;
    glm::vec3 tUpper = (box.upper - ray.origin) / ray.direction;

    float tLowerCoeff = glm::max(
        glm::max(
            glm::min(tLower.x, tUpper.x),
            glm::min(tLower.y, tUpper.y)),
        glm::min(tLower.z, tUpper.z));

    float tUpperCoeff = glm::min(
        glm::min(
            glm::max(tLower.x, tUpper.x),
            glm::max(tLower.y, tUpper.y)),
        glm::max(tLower.z, tUpper.z));

    if (tLowerCoeff > tUpperCoeff || tUpperCoeff < 0)
        return false;
    else if (tLowerCoeff < 0 && tUpperCoeff < ray.t) {
        ray.t = tUpperCoeff;
        return true;
    } else if (tLowerCoeff < ray.t) {
        ray.t = tLowerCoeff;
        return true;
    }

    return false;
}
