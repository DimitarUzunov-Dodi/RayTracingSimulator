#include "texture.h"
#include <framework/image.h>

glm::vec3 acquireTexel(const Image& image, const glm::vec2& texCoord, const Features& features)
{
    // TODO: implement this function.
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)

     auto xCoord = (int)(texCoord.x * image.width);
     auto yCoord = (int)((1-texCoord.y) * image.height);

    if (xCoord < 0) {
        xCoord = 0;
    }
    if (xCoord >= image.width) {
        xCoord = image.width-1;
    }
    if (yCoord < 0) {
        yCoord = 0;
    }
    if (yCoord >= image.height) {
        yCoord = image.height-1;
    }
    //if (xCoord > 0 && xCoord < image.width && yCoord > 0 && yCoord < image.height) {
        const auto coord = image.pixels[yCoord * image.width + xCoord];
    //} 
    
    return coord;
}