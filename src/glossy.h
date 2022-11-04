#pragma once
#include "common.h"
#include <framework/ray.h>

//Computes a random ray to use for the glossy feature
const Ray computePerturbedRay(Ray ray, HitInfo hitInfo, const Features& features);

const Ray computeGlossyRay(Ray ray, HitInfo hitInfo, const Features& features);
 