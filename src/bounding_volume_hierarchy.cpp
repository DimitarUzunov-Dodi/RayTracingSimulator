#include "bounding_volume_hierarchy.h"
#include "draw.h"
#include "intersect.h"
#include "scene.h"
#include "texture.h"
#include "interpolate.h"
#include "bounding_volume_hierarchy.h"
#include <glm/glm.hpp>
#include <queue>
#include <iostream>

float max(float f1, float f2) {
    return f1 > f2 ? f1 : f2;
}


float min(float f1, float f2)
{
    return f1 < f2 ? f1 : f2;
}

AxisAlignedBox getBounds(std::vector<Triangle> &triangles, Scene* scene) {
    float xmax = -FLT_MAX;
    float xmin = FLT_MAX;
    float ymax = -FLT_MAX;
    float ymin = FLT_MAX;
    float zmax = -FLT_MAX;
    float zmin = FLT_MAX;

    for (const auto& tri : triangles) {
        const auto v0 = scene->meshes.at(tri.mesh).vertices.at(scene->meshes.at(tri.mesh).triangles.at(tri.triangle).x);
        const auto v1 = scene->meshes.at(tri.mesh).vertices.at(scene->meshes.at(tri.mesh).triangles.at(tri.triangle).y);
        const auto v2 = scene->meshes.at(tri.mesh).vertices.at(scene->meshes.at(tri.mesh).triangles.at(tri.triangle).z);

        xmax = max(xmax, v0.position.x);
        xmax = max(xmax, v1.position.x);
        xmax = max(xmax, v2.position.x);

        ymax = max(ymax, v0.position.y);
        ymax = max(ymax, v1.position.y);
        ymax = max(ymax, v2.position.y);

        zmax = max(zmax, v0.position.z);
        zmax = max(zmax, v1.position.z);
        zmax = max(zmax, v2.position.z);

        xmin = min(xmin, v0.position.x);
        xmin = min(xmin, v1.position.x);
        xmin = min(xmin, v2.position.x);

        ymin = min(ymin, v0.position.y);
        ymin = min(ymin, v1.position.y);
        ymin = min(ymin, v2.position.y);

        zmin = min(zmin, v0.position.z);
        zmin = min(zmin, v1.position.z);
        zmin = min(zmin, v2.position.z);
    }
    AxisAlignedBox ret;
    ret.lower = glm::vec3(xmin, ymin, zmin);
    ret.upper = glm::vec3(xmax, ymax, zmax);
    return ret;
}

BoundingVolumeHierarchy::BoundingVolumeHierarchy(Scene* pScene)
    : m_pScene(pScene)
{

    std::vector<Triangle> triangles;
    int m = 0;
    for (const auto& mesh : m_pScene->meshes) {
        for (int t = 0; t < mesh.triangles.size(); t++) { 
            triangles.push_back(Triangle(m, t));
        }
        m++;
    }
    AxisAlignedBox box = getBounds(triangles, m_pScene);
    Node root;
    root.bounds = box;
    root.internal = false;
    root.triangles = triangles;
    root.level = 0;
    m_numLevels = 1;
    m_numLeaves = 1;
    nodes.push_back(root);
    subdivide(20, 0);
}


void BoundingVolumeHierarchy::subdivide(int limit, int node)
{
    int level = nodes.at(node).level;
    bool divided = false;
    if (level >= limit)
        return;


    std::vector<Triangle> otherHalf;
    
    int size = nodes.at(node).triangles.size();
    if (size > 1) {

        sortByAxis(nodes.at(node).triangles, level % 3, m_pScene);
        divided = true;
        m_numLeaves ++;

        std::vector<Triangle> left;
        std::vector<Triangle> right;
        std::vector<Triangle>::iterator middleItr(nodes.at(node).triangles.begin() + nodes.at(node).triangles.size() / 2);

        for (auto it = nodes.at(node).triangles.begin(); it != nodes.at(node).triangles.end(); ++it) {
            if (std::distance(it, middleItr) > 0) {
                left.push_back(*it);
            } else {
                right.push_back(*it);
            }
        }

        Node n1;
        n1.internal = false;
        n1.triangles = left;
        n1.level = level + 1;
        n1.bounds = getBounds(n1.triangles, m_pScene);

        Node n2;
        n2.internal = false;
        n2.triangles = right;
        n2.level = level + 1;
        n2.bounds = getBounds(n2.triangles, m_pScene);

        nodes.at(node).internal = true;       
        nodes.at(node).triangles.clear();
        nodes.at(node).children.push_back(nodes.size());
        int n1pos = nodes.size();
        nodes.push_back(n1);
        nodes.at(node).children.push_back(nodes.size());

        int n2pos = nodes.size();
        nodes.push_back(n2);

        if(n1.triangles.size() > 1) subdivide(limit, n1pos);
        if(n2.triangles.size() > 1) subdivide(limit, n2pos);
       
    }

    if (divided) {
        m_numLevels = max(m_numLevels, level + 1);
    }
}

// Return the depth of the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 1.
int BoundingVolumeHierarchy::numLevels() const
{
    return m_numLevels;
}

// Return the number of leaf nodes in the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 2.
int BoundingVolumeHierarchy::numLeaves() const
{
    return m_numLeaves;
}

// Use this function to visualize your BVH. This is useful for debugging. Use the functions in
// draw.h to draw the various shapes. We have extended the AABB draw functions to support wireframe
// mode, arbitrary colors and transparency.
void BoundingVolumeHierarchy::debugDrawLevel(int level)
{
    // Draw the AABB as a transparent green box.
    //AxisAlignedBox aabb{ glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    //drawShape(aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);
    for (auto& n : nodes) {
        if (n.level == level) {
            drawAABB(n.bounds, DrawMode::Wireframe, glm::vec3(1.00f, 0.05f, 0.05f), 0.9f);
        }
    }
}


// Use this function to visualize your leaf nodes. This is useful for debugging. The function
// receives the leaf node to be draw (think of the ith leaf node). Draw the AABB of the leaf node and all contained triangles.
// You can draw the triangles with different colors. NoteL leafIdx is not the index in the node vector, it is the
// i-th leaf node in the vector.
void BoundingVolumeHierarchy::debugDrawLeaf(int leafIdx)
{
    // Draw the AABB as a transparent green box.
    //AxisAlignedBox aabb{ glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    //drawShape(aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);

    // Draw the AABB as a (white) wireframe box.
    //drawAABB(aabb, DrawMode::Wireframe);
    int i = 1;
    for (auto& n : nodes) {
        if (n.internal == false) {
            if (i == leafIdx) {
                drawAABB(n.bounds, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.3f);
                for (auto& t : n.triangles) {
                    drawTriangle(m_pScene->meshes.at(t.mesh).vertices.at(m_pScene->meshes.at(t.mesh).triangles.at(t.triangle).x), 
                                 m_pScene->meshes.at(t.mesh).vertices.at(m_pScene->meshes.at(t.mesh).triangles.at(t.triangle).y), 
                                 m_pScene->meshes.at(t.mesh).vertices.at(m_pScene->meshes.at(t.mesh).triangles.at(t.triangle).z));
                }
            }
            i++;
        }
    }
    //drawAABB(aabb, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);

    // once you find the leaf node, you can use the function drawTriangle (from draw.h) to draw the contained primitives
}

Scene* tempScene;
int tempAxis;

struct CompareTriangles {
    CompareTriangles(Scene* tempScene, int tempAxis)
    {
        this->tempScene = tempScene;
        this->tempAxis = tempAxis;
    }
    bool operator()(Triangle& t1, Triangle& t2)
    {
        auto mesh1 = &tempScene->meshes.at(t1.mesh);

        glm::vec3 center1 = (mesh1->vertices.at(mesh1->triangles.at(t1.triangle).x).position + 
                             mesh1->vertices.at(mesh1->triangles.at(t1.triangle).y).position + 
                             mesh1->vertices.at(mesh1->triangles.at(t1.triangle).z).position) / 3.0f;

        auto mesh2 = &tempScene->meshes.at(t2.mesh);


        glm::vec3 center2 = (mesh2->vertices.at(mesh2->triangles.at(t2.triangle).x).position + 
                             mesh2->vertices.at(mesh2->triangles.at(t2.triangle).y).position + 
                             mesh2->vertices.at(mesh2->triangles.at(t2.triangle).z).position) / 3.0f;

        switch (tempAxis) {
        case 0:
            return center1.x < center2.x;
        case 1:
            return center1.y < center2.y;
        case 2:
            return center1.z < center2.z;
        default:
            return false;
        }
    }

    Scene *tempScene;
    int tempAxis;
};


void BoundingVolumeHierarchy::sortByAxis(std::vector<Triangle>& triangles, int axis, Scene* scene) {
    tempScene = scene;
    tempAxis = axis;
    sort(triangles.begin(), triangles.end(), CompareTriangles(scene, axis));
}

// Return true if something is hit, returns false otherwise. Only find hits if they are closer than t stored
// in the ray and if the intersection is on the correct side of the origin (the new t >= 0). Replace the code
// by a bounding volume hierarchy acceleration structure as described in the assignment. You can change any
// file you like, including bounding_volume_hierarchy.h.
bool BoundingVolumeHierarchy::intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const
{
    // If BVH is not enabled, use the naive implementation.
    if (!features.enableAccelStructure) {
        bool hit = false;
        // Intersect with all triangles of all meshes.
        for (const auto& mesh : m_pScene->meshes) {
            for (const auto& tri : mesh.triangles) {
                const auto v0 = mesh.vertices[tri[0]];
                const auto v1 = mesh.vertices[tri[1]];
                const auto v2 = mesh.vertices[tri[2]];
                if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                    hitInfo.material = mesh.material;
                    hitInfo.normal = glm::normalize(glm::cross(v1.position - v0.position, v2.position - v0.position));
                    hit = true;
                }
            }
        }
        // Intersect with spheres.
        for (const auto& sphere : m_pScene->spheres)
            hit |= intersectRayWithShape(sphere, ray, hitInfo);
        return hit;
    } else {
        // TODO: implement here the bounding volume hierarchy traversal.
        // Please note that you should use `features.enableNormalInterp` and `features.enableTextureMapping`
        // to isolate the code that is only needed for the normal interpolation and texture mapping features.

        //nodes 
        bool hit = false;
        std::queue<Node> correctNodes;
        std::queue<Node> tree;
        std::queue<Node> coloredNodes;
        Ray ourRay = ray;

        
        float originalT = ray.t;

        // Traverses all the nodes in the tree and ads the final nodes to a queue
        if (intersectRayWithShape(nodes.at(0).bounds, ray)) {
            ray.t = originalT;
            tree.push(nodes.at(0));
            
            while (tree.size() != 0) {

                Node curr = tree.front();
                tree.pop();

                if (intersectRayWithShape(curr.bounds, ray)) {
                    if (correctNodes.empty())  
                         correctNodes.push(nodes.at(0));

                    drawAABB(curr.bounds, DrawMode::Wireframe, glm::vec3(0.00f, 1.0f, 0.5f), 1.0f);
                    for (auto n : curr.children) {
                        tree.push(nodes.at(n));
                        if (curr.internal) {
                            correctNodes.push(nodes.at(n));
                        }
                    }
                    ray.t = originalT;
                } else {
                    drawAABB(curr.bounds, DrawMode::Wireframe, glm::vec3(1.00f, 0.5f, 0.0f), 0.3f);
                    ray.t = originalT;
                }
            }
        }
        

        float hitT = originalT;
        Node hitNode;
       

        Vertex v1;
        Vertex v2;
        Vertex v3;
        Node correctNode;

        //finds the node in which the hit happens
        while (correctNodes.size() > 0) {
            Node curr = correctNodes.front();
            correctNodes.pop();
            if (!curr.triangles.empty()) {
                for (auto t : curr.triangles) {
                    v1 = m_pScene->meshes.at(t.mesh).vertices.at(m_pScene->meshes.at(t.mesh).triangles.at(t.triangle).x);
                    v2 = m_pScene->meshes.at(t.mesh).vertices.at(m_pScene->meshes.at(t.mesh).triangles.at(t.triangle).y);
                    v3 = m_pScene->meshes.at(t.mesh).vertices.at(m_pScene->meshes.at(t.mesh).triangles.at(t.triangle).z);

                    if (intersectRayWithTriangle(v1.position, v2.position, v3.position, ray, hitInfo)) {
                        if (hitT > ray.t) {
                            hitT = ray.t;
                            correctNode = curr;
                            hit = true;
                        }
                    } else {
                        intersectRayWithShape(curr.bounds, ray);
                        if (hitT >= ray.t) {
                            coloredNodes.push(curr);
                            ray.t = originalT;
                        }
                    }
                }
            }
        }
        if (hit) {
            //intersectRayWithTriangle(tX.position, tY.position, tZ.position, ray, hitInfo)
            auto t = correctNode.triangles.at(0);
            v1 = m_pScene->meshes.at(t.mesh).vertices.at(m_pScene->meshes.at(t.mesh).triangles.at(t.triangle).x);
            v2 = m_pScene->meshes.at(t.mesh).vertices.at(m_pScene->meshes.at(t.mesh).triangles.at(t.triangle).y);
            v3 = m_pScene->meshes.at(t.mesh).vertices.at(m_pScene->meshes.at(t.mesh).triangles.at(t.triangle).z);
            drawTriangle(v1, v2, v3);
            ray.t = originalT;
            intersectRayWithTriangle(v1.position, v2.position, v3.position, ray, hitInfo);
        }
        
        return hit;
    }
}




Triangle::Triangle(int mesh, int triangle)
{
    this->mesh = mesh;
    this->triangle = triangle;
}

Node::Node()
{
}