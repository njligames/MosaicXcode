//
//  Color.cpp
//  SDLSkeleton
//
//  Created by James Folk on 3/31/20.
//

#include "Color.h"

#include <algorithm>

Color::Color() : mRGBColor(glm::vec4(1, 0, 0, 1)) {}
Color::Color(const Color &rhs) : mRGBColor(rhs.mRGBColor) {}
const Color &Color::operator=(const Color &rhs) {
    if (this != &rhs) {
        mRGBColor = rhs.mRGBColor;
    }
    return *this;
}

Color::~Color() {}

void Color::RGBtoHSV(const float fR, const float fG, const float fB, float &fH,
                     float &fS, float &fV) {

    float fCMax = std::max(std::max(fR, fG), fB);
    float fCMin = std::min(std::min(fR, fG), fB);
    float fDelta = fCMax - fCMin;

    if (fDelta > 0) {
        if (fCMax == fR) {
            fH = 60 * (fmod(((fG - fB) / fDelta), 6));
        } else if (fCMax == fG) {
            fH = 60 * (((fB - fR) / fDelta) + 2);
        } else if (fCMax == fB) {
            fH = 60 * (((fR - fG) / fDelta) + 4);
        }

        if (fCMax > 0) {
            fS = fDelta / fCMax;
        } else {
            fS = 0;
        }

        fV = fCMax;
    } else {
        fH = 0;
        fS = 0;
        fV = fCMax;
    }

    if (fH < 0) {
        fH = 360 + fH;
    }
}

void Color::HSVtoRGB(const float fH, const float fS, const float fV, float &fR,
                     float &fG, float &fB) {

    float fC = fV * fS; // Chroma
    float fHPrime = fmod(fH / 60.0, 6);
    float fX = fC * (1 - fabs(fmod(fHPrime, 2) - 1));
    float fM = fV - fC;

    if (0 <= fHPrime && fHPrime < 1) {
        fR = fC;
        fG = fX;
        fB = 0;
    } else if (1 <= fHPrime && fHPrime < 2) {
        fR = fX;
        fG = fC;
        fB = 0;
    } else if (2 <= fHPrime && fHPrime < 3) {
        fR = 0;
        fG = fC;
        fB = fX;
    } else if (3 <= fHPrime && fHPrime < 4) {
        fR = 0;
        fG = fX;
        fB = fC;
    } else if (4 <= fHPrime && fHPrime < 5) {
        fR = fX;
        fG = 0;
        fB = fC;
    } else if (5 <= fHPrime && fHPrime < 6) {
        fR = fC;
        fG = 0;
        fB = fX;
    } else {
        fR = 0;
        fG = 0;
        fB = 0;
    }

    fR += fM;
    fG += fM;
    fB += fM;
}

glm::vec4 Color::rgb2hsv(const glm::vec4 &rgb)
// void rgb2hsv(const unsigned char &src_r, const unsigned char &src_g, const
// unsigned char &src_b, unsigned char &dst_h, unsigned char &dst_s, unsigned
// char &dst_v)
{

    float h, s, v;

    RGBtoHSV(rgb.x, rgb.y, rgb.z, h, s, v);

    return glm::vec4(h, s, v, rgb.w);
}

glm::vec4 Color::hsv2rgb(const glm::vec4 &hsv) {
    float dst_r;
    float dst_g;
    float dst_b;

    HSVtoRGB(hsv.x, hsv.y, hsv.z, dst_r, dst_g, dst_b);

    return glm::vec4(dst_r, dst_g, dst_b, hsv.w);
}

Color Color::modifyHue(const Color &color, const float h) {
    Color c(color);
    glm::vec4 hsv(c.getHSV());
    c.setHSV(hsv.x + h, hsv.y, hsv.z);
    return c;
}

Color Color::modifySaturation(const Color &color, const float s) {
    Color c(color);
    glm::vec4 hsv(c.getHSV());
    float _s = hsv.y + s;
    hsv.y = _s;
    if (_s < 0)
        hsv.y = 0;
    if (_s > 1)
        hsv.y = 1;
    return c;
}

Color Color::modifyValue(const Color &color, const float v) {
    Color c(color);
    glm::vec4 hsv(c.getHSV());
    float _v = hsv.z + v;
    hsv.z = _v;
    if (_v < 0)
        hsv.z = 0;
    if (_v > 1)
        hsv.z = 1;
    return c;
}

void Color::setHSV(const float fH, const float fS, const float fV) {
    float fR, fG, fB;
    HSVtoRGB(fH, fS, fV, fR, fG, fB);
    mRGBColor = glm::vec4(fR, fG, fB, mRGBColor.a);
}

void Color::setHSV(const glm::vec4 &rgb) {}

void Color::setRGB(const float fR, const float fG, const float fB) {
    mRGBColor = glm::vec4(fR, fG, fB, mRGBColor.a);
}

void Color::setRGB(const std::string &hex) {
    std::string _hex(hex);
    while (_hex.size() < 6) {
        _hex.append("0");
    }

    int r, g, b;
    sscanf(_hex.c_str(), "%02x%02x%02x", &r, &g, &b);
    mRGBColor.r = r / 255.f;
    mRGBColor.g = g / 255.f;
    mRGBColor.b = b / 255.f;
}
void Color::setRGB(const int hex) {
    mRGBColor.r = ((hex >> 16) & 0xFF) / 255.0; // Extract the RR byte
    mRGBColor.g = ((hex >> 8) & 0xFF) / 255.0;  // Extract the GG byte
    mRGBColor.b = ((hex)&0xFF) / 255.0;         // Extract the BB byte
}
void Color::setRGB(const glm::vec4 &rgb) {
    mRGBColor.r = (rgb.r > 1.f) ? 1.f : (rgb.r < 0.f) ? 0.f : rgb.r;
    mRGBColor.g = (rgb.g > 1.f) ? 1.f : (rgb.g < 0.f) ? 0.f : rgb.g;
    mRGBColor.b = (rgb.b > 1.f) ? 1.f : (rgb.b < 0.f) ? 0.f : rgb.b;
    mRGBColor.a = (rgb.a > 1.f) ? 1.f : (rgb.a < 0.f) ? 0.f : rgb.a;
}

glm::vec4 Color::getHSV() const {
    float h, s, v;
    RGBtoHSV(mRGBColor.r, mRGBColor.g, mRGBColor.b, h, s, v);
    return glm::vec4(h, s, v, mRGBColor.a);
}
glm::vec4 Color::getRGB() const { return mRGBColor; }

float Color::red() const { return getRGB().x; }
float Color::green() const { return getRGB().y; }
float Color::blue() const { return getRGB().z; }

float Color::hue() const { return getHSV().x; }
float Color::saturation() const { return getHSV().y; }
float Color::value() const { return getHSV().z; }
