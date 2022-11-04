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
#include <stack>
#include <functional>

float max(float f1, float f2) {
    return f1 > f2 ? f1 : f2;
}


float min(float f1, float f2)
{
    return f1 < f2 ? f1 : f2;
}

AxisAlignedBox getBounds(std::vector<TriangleOrChild> &triangles, Scene* scene) {
    float xmax = -FLT_MAX;
    float xmin = FLT_MAX;
    float ymax = -FLT_MAX;
    float ymin = FLT_MAX;
    float zmax = -FLT_MAX;
    float zmin = FLT_MAX;

    for (const auto& tri : triangles) {
        const auto v0 = scene->meshes.at(tri.meshOrChild).vertices.at(scene->meshes.at(tri.meshOrChild).triangles.at(tri.triangle).x);
        const auto v1 = scene->meshes.at(tri.meshOrChild).vertices.at(scene->meshes.at(tri.meshOrChild).triangles.at(tri.triangle).y);
        const auto v2 = scene->meshes.at(tri.meshOrChild).vertices.at(scene->meshes.at(tri.meshOrChild).triangles.at(tri.triangle).z);

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

    std::vector<TriangleOrChild> triangles;
    int m = 0;
    int bytes = 0;
    for (const auto& mesh : m_pScene->meshes) {
        bytes += mesh.triangles.size();
    }
    triangles.reserve(bytes);
    for (const auto& mesh : m_pScene->meshes) {
        for (int t = 0; t < mesh.triangles.size(); t++) { 
            triangles.push_back(TriangleOrChild(m, t));
        }
        m++;
    }
    AxisAlignedBox box = getBounds(triangles, m_pScene);
    Node root;
    root.bounds = box;
    root.internal = false;
    root.trianglesOrChildren = triangles;
    root.level = 0;
    m_numLevels = 1;
    m_numLeaves = 1;
    nodes.push_back(root);
    subdivide(MAX_LEVELS, 0);
}

void BoundingVolumeHierarchy::subdivide(int limit, int node)
{
    int level = nodes.at(node).level;
    bool divided = false;
    if (level >= limit)
        return;


    std::vector<TriangleOrChild> otherHalf;
    
    int size = nodes.at(node).trianglesOrChildren.size();
    if (size > 1) {

        quickSortStepByAxis(nodes.at(node).trianglesOrChildren, level % 3, m_pScene);
        divided = true;
        m_numLeaves ++;

        std::vector<TriangleOrChild> left;
        std::vector<TriangleOrChild> right;
        std::vector<TriangleOrChild>::iterator middleItr = nodes.at(node).trianglesOrChildren.begin() + nodes.at(node).trianglesOrChildren.size() / 2;

        for (auto it = nodes.at(node).trianglesOrChildren.begin(); it != nodes.at(node).trianglesOrChildren.end(); ++it) {
            if (std::distance(it, middleItr) > 0) {
                left.push_back(*it);
            } else {
                right.push_back(*it);
            }
        }

        Node n1;
        n1.internal = false;
        n1.trianglesOrChildren = left;
        n1.level = level + 1;
        n1.bounds = getBounds(n1.trianglesOrChildren, m_pScene);

        Node n2;
        n2.internal = false;
        n2.trianglesOrChildren = right;
        n2.level = level + 1;
        n2.bounds = getBounds(n2.trianglesOrChildren, m_pScene);

        nodes.at(node).internal = true;       
        nodes.at(node).trianglesOrChildren.clear();
        nodes.at(node).trianglesOrChildren.push_back(TriangleOrChild(nodes.size(), 0));
        int n1pos = nodes.size();
        nodes.push_back(n1);
        nodes.at(node).trianglesOrChildren.push_back(TriangleOrChild(nodes.size(), 0));

        int n2pos = nodes.size();
        nodes.push_back(n2);

        if (n1.trianglesOrChildren.size() > 1)
            subdivide(limit, n1pos);
        if (n2.trianglesOrChildren.size() > 1)
            subdivide(limit, n2pos);
       
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
                for (auto& t : n.trianglesOrChildren) {
                    drawTriangle(m_pScene->meshes.at(t.meshOrChild).vertices.at(m_pScene->meshes.at(t.meshOrChild).triangles.at(t.triangle).x), 
                                 m_pScene->meshes.at(t.meshOrChild).vertices.at(m_pScene->meshes.at(t.meshOrChild).triangles.at(t.triangle).y), 
                                 m_pScene->meshes.at(t.meshOrChild).vertices.at(m_pScene->meshes.at(t.meshOrChild).triangles.at(t.triangle).z));
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
bool compareTriangles(TriangleOrChild& t1, TriangleOrChild&t2) {

    auto mesh1 = &tempScene->meshes.at(t1.meshOrChild);
    auto mesh2 = &tempScene->meshes.at(t2.meshOrChild);
    float c1, c2;
    switch (tempAxis) {
    case 0:
        c1 = (mesh1->vertices.at(mesh1->triangles.at(t1.triangle).x).position.x + mesh1->vertices.at(mesh1->triangles.at(t1.triangle).y).position.x + mesh1->vertices.at(mesh1->triangles.at(t1.triangle).z).position.x); 
       
        c2 = (mesh2->vertices.at(mesh2->triangles.at(t2.triangle).x).position.x + mesh2->vertices.at(mesh2->triangles.at(t2.triangle).y).position.x + mesh2->vertices.at(mesh2->triangles.at(t2.triangle).z).position.x); 
        return c1 < c2;
    case 1:
        c1 = (mesh1->vertices.at(mesh1->triangles.at(t1.triangle).x).position.y + mesh1->vertices.at(mesh1->triangles.at(t1.triangle).y).position.y + mesh1->vertices.at(mesh1->triangles.at(t1.triangle).z).position.y);

        c2 = (mesh2->vertices.at(mesh2->triangles.at(t2.triangle).x).position.y + mesh2->vertices.at(mesh2->triangles.at(t2.triangle).y).position.y + mesh2->vertices.at(mesh2->triangles.at(t2.triangle).z).position.y);
        return c1 < c2;
    case 2:
        c1 = (mesh1->vertices.at(mesh1->triangles.at(t1.triangle).x).position.z + mesh1->vertices.at(mesh1->triangles.at(t1.triangle).y).position.z + mesh1->vertices.at(mesh1->triangles.at(t1.triangle).z).position.z);

        c2 = (mesh2->vertices.at(mesh2->triangles.at(t2.triangle).x).position.z + mesh2->vertices.at(mesh2->triangles.at(t2.triangle).y).position.z + mesh2->vertices.at(mesh2->triangles.at(t2.triangle).z).position.z);
        return c1 < c2;
    default:
        return false;
    }
}


void BoundingVolumeHierarchy::quickSortStepByAxis(std::vector<TriangleOrChild>& triangles, int axis, Scene* scene) {
    tempScene = scene;
    tempAxis = axis;
    std::vector<TriangleOrChild> left;
    std::vector<TriangleOrChild> right;
    std::vector<TriangleOrChild>::iterator pivot(triangles.begin() + triangles.size() / 2);
    std::nth_element(triangles.begin(), triangles.begin() + triangles.size() / 2, triangles.end(), compareTriangles);
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
        Ray norm;
        
        Ray normal;
        Ray normalA;
        Ray normalB;
        Ray normalC;
       
       
        // Intersect with all triangles of all meshes.
        for (const auto& mesh : m_pScene->meshes) {
            for (const auto& tri : mesh.triangles) {
                const auto v0 = mesh.vertices[tri[0]];
                const auto v1 = mesh.vertices[tri[1]];
                const auto v2 = mesh.vertices[tri[2]];
                if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {

                    hitInfo.material = mesh.material;
                    hitInfo.normal = glm::normalize(glm::cross(v1.position - v0.position, v2.position - v0.position));
                    hitInfo.barycentricCoord = computeBarycentricCoord(v0.position, v1.position, v2.position, ray.origin + ray.direction * ray.t);
                    hitInfo.texCoord = interpolateTexCoord(v0.texCoord, v1.texCoord, v2.texCoord, hitInfo.barycentricCoord);
                    hit = true;

                    norm = Ray(ray.origin + ray.direction * ray.t, hitInfo.normal, 0.5f);
                    
                        normalA = Ray(v0.position, v0.normal, 0.5f);
                        normalB = Ray(v1.position, v1.normal, 0.5f);
                        normalC = Ray(v2.position, v2.normal, 0.5f);
                        normal = Ray(ray.origin + ray.direction * ray.t, interpolateNormal(v0.normal, v1.normal, v2.normal, hitInfo.barycentricCoord), 0.5f);
                      
                }
            }
            if (features.enableNormalInterp && hit == true) {
                hitInfo.normal = normal.direction;
            }
        }
        
        if (features.enableNormalInterp && hit == true) {

            drawRay(normalA, glm::vec3(0.0f, 1.0f, 1.0f));
            drawRay(normalB, glm::vec3(0.0f, 1.0f, 0.5f));
            drawRay(normalC, glm::vec3(0.0f, 0.0f, 1.0f));
            drawRay(normal, glm::vec3(0.0f, 0.5f, 1.0f));

        } 
        // Intersect with spheres.
        for (const auto& sphere : m_pScene->spheres)
            hit |= intersectRayWithShape(sphere, ray, hitInfo);
        return hit;
    } else {
        // TODO: implement here the bounding volume hierarchy traversal.
        // Please note that you should use `features.enableNormalInterp` and `features.enableTextureMapping`
        // to isolate the code that is only needed for the normal interpolation and texture mapping features.

        // nodes
        bool hit = false;
        std::queue<Node> colloredNodes;
        std::queue<Node> tree;
        
        Ray ourRay = ray;
        std::stack < std::pair<float, int> > traverse ;

        Vertex v1;
        Vertex v2;
        Vertex v3;
        Node correctNode;

        float originalT = ray.t;

       

        if (intersectRayWithShape(nodes.at(0).bounds, ray)) {

            traverse.push(std::pair(ray.t,0));
            
            float hitT = originalT;

            while (!traverse.empty() ) {
                
                ray.t = originalT;
                const Node& curr = nodes.at(traverse.top().second);
                float dist = traverse.top().first;
                if (enableDebugDraw) {
                     drawAABB(curr.bounds, DrawMode::Wireframe, glm::vec3(0.00f, 1.f, 0.50f), 0.3f);
                }
                traverse.pop();
                if (hitT < dist) {
                    continue;
                }

                if (curr.internal) {
                    float hit1;
                    float hit2;
                    bool intr1 = false;
                    bool intr2 = true;
                    const Node& child1 = nodes.at(curr.children[0]);
                    const Node& child2 = nodes.at(curr.children[1]);
                    
                    if (intersectRayWithShape(child1.bounds, ray)) { 
                        hit1 = ray.t;
                        ray.t = originalT;
                        intr1 = true;   
                    }
                    if (intersectRayWithShape(child2.bounds, ray)) {
                        hit2 = ray.t;
                        ray.t = originalT;
                        intr2 = true;
                    }

                    if (intr1 && intr2) {
                        if (hit1 < hit2) {
                            traverse.push(std::pair(hit2, curr.children[1]));
                            traverse.push(std::pair(hit1, curr.children[0]));
                        } else {
                            traverse.push(std::pair(hit1, curr.children[0]));
                            traverse.push(std::pair(hit2, curr.children[1]));
                        }
                    } else if (intr1) {
                        traverse.push(std::pair(hit1, curr.children[0]));
                    } else if (intr2) {
                        traverse.push(std::pair(hit2, curr.children[1]));
                    }
                    
                    
                                                                
                } else {
                    if (!curr.triangles.empty()) { 

                    
                        for (const Triangle &t : curr.triangles) {

                        const Mesh& mesh = m_pScene->meshes[t.mesh];
                        const glm::vec3& tr = mesh.triangles.at(t.triangle);

                        v1 = m_pScene->meshes.at(t.mesh).vertices.at(m_pScene->meshes.at(t.mesh).triangles.at(t.triangle).x);
                        v2 = m_pScene->meshes.at(t.mesh).vertices.at(m_pScene->meshes.at(t.mesh).triangles.at(t.triangle).y);
                        v3 = m_pScene->meshes.at(t.mesh).vertices.at(m_pScene->meshes.at(t.mesh).triangles.at(t.triangle).z);
                        ray.t = originalT;
                        if (intersectRayWithTriangle(v1.position, v2.position, v3.position, ray, hitInfo)) {
                            if (hitT > ray.t) {
                                hitT = ray.t;
                                correctNode = curr;
                                hitInfo.material = m_pScene->meshes.at(t.mesh).material;
                                hitInfo.normal = glm::normalize(glm::cross(v2.position - v1.position, v3.position - v1.position));
                                hitInfo.barycentricCoord = computeBarycentricCoord(v1.position, v2.position, v3.position, ray.origin + ray.direction * ray.t);
                                hitInfo.texCoord = interpolateTexCoord(v1.texCoord, v2.texCoord, v3.texCoord, hitInfo.barycentricCoord);
                                hit = true;

                           } // TODO add the normalInterpolation case
                        }

                    }
                    
                }    
            }


        }
            ray.t = originalT;
        
            if (hit) {
              
                auto t = correctNode.triangles.at(0);
                v1 = m_pScene->meshes.at(t.mesh).vertices.at(m_pScene->meshes.at(t.mesh).triangles.at(t.triangle).x);
                v2 = m_pScene->meshes.at(t.mesh).vertices.at(m_pScene->meshes.at(t.mesh).triangles.at(t.triangle).y);
                v3 = m_pScene->meshes.at(t.mesh).vertices.at(m_pScene->meshes.at(t.mesh).triangles.at(t.triangle).z);
                if (enableDebugDraw) {
                    drawTriangle(v1, v2, v3);
                }
                ray.t = originalT;

                intersectRayWithTriangle(v1.position, v2.position, v3.position, ray, hitInfo);
            }

           
        }
        return hit;
    }
    
}

TriangleOrChild::TriangleOrChild(int mesh, int triangle)
{
    this->meshOrChild = mesh;
    this->triangle = triangle;
}

Node::Node()
{
}

