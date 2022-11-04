#pragma once
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
DISABLE_WARNINGS_POP()
#include <filesystem>
#include <vector>
#include <framework/trackball.h>

class Screen {
public:
    Screen(const glm::ivec2& resolution, bool presentable = true);

    void clear(const glm::vec3& color);
    void setPixel(int x, int y, const glm::vec3& color);
    void setPreviousCameraMatrix(const Trackball& camera);
    void setVelocityBuffer(const glm::ivec2& pixelPosition, const glm::ivec2& resolution, const glm::vec3& worldPosition, int sampleCount, bool hit);
    void initVelocityBuffer(int size);
    glm::vec3 getPixelVelocity(int x, int y);
    void motionBlur(int sampleCount, float strength);
    void debugMotionBlur(int sampleCount, float strength, float density);
    void writeBitmapToFile(const std::filesystem::path& filePath);
    void draw();

    [[nodiscard]] glm::ivec2 resolution() const;

    /// Calculates the index of a pixel in the `m_textureData` vector.
    /// Pixels are stored from bottom to top, left to right (this is to facilitate writing as a bmp).
    [[nodiscard]] int indexAt(int x, int y) const;
    [[nodiscard]] const std::vector<glm::vec3>& pixels() const;
    [[nodiscard]] std::vector<glm::vec3>& pixels();

    
    bool firstRender = true;

private:
    bool m_presentable;
    glm::ivec2 m_resolution;
    std::vector<glm::vec3> m_textureData;
    uint32_t m_texture;

    // Added by Filip
    glm::mat4 m_previousViewMatrix;
    glm::mat4 m_previousProjectionMatrix;
    std::vector<glm::vec2> m_velocityBuffer;
};
