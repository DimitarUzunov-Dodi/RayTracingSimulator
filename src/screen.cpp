#include "screen.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/common.hpp>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb/stb_image_write.h>
DISABLE_WARNINGS_POP()
#include <algorithm>
#include <framework/opengl_includes.h>
#include <string>
#include "glm/gtx/string_cast.hpp"
#include <iostream>

Screen::Screen(const glm::ivec2& resolution, bool presentable)
    : m_presentable(presentable)
    , m_resolution(resolution)
    , m_textureData(size_t(resolution.x * resolution.y), glm::vec3(0.0f))
{
    // Create OpenGL texture if we want to present the screen.
    if (m_presentable) {
        // Generate texture
        glGenTextures(1, &m_texture);
        glBindTexture(GL_TEXTURE_2D, m_texture);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glBindTexture(GL_TEXTURE_2D, 0);
    }
    firstRender = true;
}

void Screen::clear(const glm::vec3& color)
{
    std::fill(std::begin(m_textureData), std::end(m_textureData), color);
}

void Screen::setPixel(int x, int y, const glm::vec3& color)
{
    // In the window/camera class we use (0, 0) at the bottom left corner of the screen (as used by GLFW).
    // OpenGL / stbi like the origin / (-1,-1) to be at the TOP left corner so transform the y coordinate.
    const int i = (m_resolution.y - 1 - y) * m_resolution.x + x;
    m_textureData[i] = glm::vec4(color, 1.0f);
}

void Screen::setPreviousCameraMatrix(const Trackball& camera)
{
    m_previousViewMatrix = camera.viewMatrix();
    m_previousProjectionMatrix = camera.projectionMatrix();
}

glm::vec2 getScreenPositionFromWorldPosition(glm::mat4 viewMatrix, glm::mat4 projectionMatrix, const glm::vec2& resolution, const glm::vec3& position)
{
    glm::vec4 normalized = projectionMatrix * (viewMatrix * glm::vec4(position, 1));
    normalized = normalized / normalized.w;
    glm::vec2 pixelPosition {
        (normalized.x + 1) / 2 * float(resolution.x),
        (normalized.y + 1) / 2 * float(resolution.y)
    };
    return pixelPosition;
}

void Screen::setVelocityBuffer(const glm::ivec2& pixelPosition, const glm::ivec2& resolution, const glm::vec3& worldPosition, int sampleCount, bool hit)
{

    const int i = (m_resolution.y - 1 - pixelPosition.y) * m_resolution.x + pixelPosition.x;
    if (!hit) {
        m_velocityBuffer.at(i) = glm::vec2(0, 0);
        return;
    }
    glm::vec2 oldPixelPosition = getScreenPositionFromWorldPosition(m_previousViewMatrix, m_previousProjectionMatrix, resolution, worldPosition);
    m_velocityBuffer.at(i) = (oldPixelPosition - glm::vec2(pixelPosition)) / (float)sampleCount;
}

void Screen::initVelocityBuffer(int size) 
{
    for (int i = 0; i < size; i++)
        m_velocityBuffer.push_back(glm::vec2(0, 0));
}

void Screen::motionBlur(int sampleCount)
{
    for (int y = 0; y < m_resolution.y; y++) {
        for (int x = 0; x != m_resolution.x; x++) {
            int i = (m_resolution.y - 1 - y) * m_resolution.x + x;
            glm::vec2 texCoords = glm::vec2(x, y) + m_velocityBuffer.at(i);    
            int j;
            glm::vec3 before = m_textureData.at(i);
            for (j = 1; j < sampleCount; j++) {
                int i2 = (m_resolution.y - 1 - (int)texCoords.y) * m_resolution.x + (int)texCoords.x;
                if (i2 < m_textureData.size())
                    m_textureData.at(i) += m_textureData.at(i2);
                else
                    break; // potential bug with dividing wrong number of samples
                texCoords += m_velocityBuffer.at(i);
            }
            m_textureData.at(i) /= j;
        }
    }
}

void Screen::writeBitmapToFile(const std::filesystem::path& filePath)
{
    std::vector<glm::u8vec4> textureData8Bits(m_textureData.size());
    std::transform(std::begin(m_textureData), std::end(m_textureData), std::begin(textureData8Bits),
        [](const glm::vec3& color) {
            const glm::vec3 clampedColor = glm::clamp(color, 0.0f, 1.0f);
            return glm::u8vec4(glm::vec4(clampedColor, 1.0f) * 255.0f);
        });

    std::string filePathString = filePath.string();
    stbi_write_bmp(filePathString.c_str(), m_resolution.x, m_resolution.y, 4, textureData8Bits.data());
}

void Screen::draw()
{
    if (m_presentable) {
        glPushAttrib(GL_ALL_ATTRIB_BITS);

        glBindTexture(GL_TEXTURE_2D, m_texture);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, m_resolution.x, m_resolution.y, 0, GL_RGB, GL_FLOAT, m_textureData.data());

        glDisable(GL_LIGHTING);
        glDisable(GL_LIGHT0);
        glDisable(GL_COLOR_MATERIAL);
        glDisable(GL_NORMALIZE);
        glColor3f(1.0f, 1.0f, 1.0f);

        glEnable(GL_TEXTURE_2D);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, m_texture);

        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        glLoadIdentity();
        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
        glLoadIdentity();

        glBegin(GL_QUADS);
        glTexCoord2f(0.0f, 1.0f);
        glVertex3f(-1.0f, -1.0f, 0.0f);
        glTexCoord2f(1.0f, 1.0f);
        glVertex3f(+1.0f, -1.0f, 0.0f);
        glTexCoord2f(1.0f, 0.0f);
        glVertex3f(+1.0f, +1.0f, 0.0f);
        glTexCoord2f(0.0f, 0.0f);
        glVertex3f(-1.0f, +1.0f, 0.0f);
        glEnd();

        glMatrixMode(GL_PROJECTION);
        glPopMatrix();
        glMatrixMode(GL_MODELVIEW);
        glPopMatrix();

        glPopAttrib();
    } else {
        std::cerr << "Screen::draw() called on non-presentable screen" << std::endl;
    }
}

glm::ivec2 Screen::resolution() const
{
    return m_resolution;
}

const std::vector<glm::vec3>& Screen::pixels() const
{
    return m_textureData;
}

std::vector<glm::vec3>& Screen::pixels()
{
    return m_textureData;
}

int Screen::indexAt(int x, int y) const
{
    return (m_resolution.y - 1 - y) * m_resolution.x + x;
}
