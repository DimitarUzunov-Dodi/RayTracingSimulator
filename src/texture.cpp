#include "texture.h"
#include <framework/image.h>

glm::vec3 acquireTexel(const Image& image, const glm::vec2& texCoord, const Features& features)
{
    // TODO: implement this function.
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)

    const auto xCoord = (int)(texCoord.x * image.width);
    const auto yCoord = (int)(texCoord.y * image.height);
    const auto coord = image.pixels[yCoord * image.width + xCoord];
    
    return coord;
}