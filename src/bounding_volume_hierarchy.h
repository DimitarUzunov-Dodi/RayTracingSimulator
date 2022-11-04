#pragma once
#include "common.h"
#include <array>
#include <framework/ray.h>
#include <vector>

#define MAX_LEVELS 20
extern bool enableDebugNonChecked;

// Forward declaration.
struct Scene;

class TriangleOrChild {
public:
    int meshOrChild;
    int triangle;
    TriangleOrChild(int meshOrChild, int triangle);
};

class Node {
public:
    Node();
    bool internal = false;
    int level = 0;
    std::vector<TriangleOrChild> trianglesOrChildren;
    AxisAlignedBox bounds;
};

class BoundingVolumeHierarchy {
public:
    // Constructor. Receives the scene and builds the bounding volume hierarchy.
    BoundingVolumeHierarchy(Scene* pScene);

    void subdivide(int limit, int node);

    static void quickSortStepByAxis(std::vector<TriangleOrChild>& triangles, int axis, Scene* scene);

    // Return how many levels there are in the tree that you have constructed.
    [[nodiscard]] int numLevels() const;

    // Return how many leaf nodes there are in the tree that you have constructed.
    [[nodiscard]] int numLeaves() const;

    // Visual Debug 1: Draw the bounding boxes of the nodes at the selected level.
    void debugDrawLevel(int level);

    // Visual Debug 2: Draw the triangles of the i-th leaf
    void debugDrawLeaf(int leafIdx);

    // Return true if something is hit, returns false otherwise.
    // Only find hits if they are closer than t stored in the ray and the intersection
    // is on the correct side of the origin (the new t >= 0).
    bool intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const;


private:
    int m_numLevels;
    int m_numLeaves;
    Scene* m_pScene;
    std::vector<Node> nodes;
};
