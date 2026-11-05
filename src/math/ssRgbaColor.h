//desc.h 
#pragma once


#include <cstdint>


namespace SolveSpace
{

 
#define RGBi(r, g, b) RgbaColor::From((r), (g), (b))
#define RGBf(r, g, b) RgbaColor::FromFloat((float)(r), (float)(g), (float)(b))

// Note: sizeof(class RgbaColor) should be exactly 4
//
class RgbaColor {
public:
    uint8_t red, green, blue, alpha;

    float redF()   const { return (float)red   / 255.0f; }
    float greenF() const { return (float)green / 255.0f; }
    float blueF()  const { return (float)blue  / 255.0f; }
    float alphaF() const { return (float)alpha / 255.0f; }

    bool IsEmpty() const { return alpha == 0; }

    bool Equals(RgbaColor c) const {
        return
            c.red   == red   &&
            c.green == green &&
            c.blue  == blue  &&
            c.alpha == alpha;
    }

    RgbaColor WithAlpha(uint8_t newAlpha) const {
        RgbaColor color = *this;
        color.alpha = newAlpha;
        return color;
    }

    uint32_t ToPackedIntBGRA() const {
        return
            blue |
            (uint32_t)(green << 8) |
            (uint32_t)(red << 16) |
            (uint32_t)((255 - alpha) << 24);
    }

    uint32_t ToPackedInt() const {
        return
            red |
            (uint32_t)(green << 8) |
            (uint32_t)(blue << 16) |
            (uint32_t)((255 - alpha) << 24);
    }

    uint32_t ToARGB32() const {
        return
            blue |
            (uint32_t)(green << 8) |
            (uint32_t)(red << 16) |
            (uint32_t)(alpha << 24);
    }

    static RgbaColor From(int r, int g, int b, int a = 255) {
        RgbaColor c;
        c.red   = (uint8_t)r;
        c.green = (uint8_t)g;
        c.blue  = (uint8_t)b;
        c.alpha = (uint8_t)a;
        return c;
    }

    static RgbaColor FromFloat(float r, float g, float b, float a = 1.0) {
        return From(
            (int)(255.1f * r),
            (int)(255.1f * g),
            (int)(255.1f * b),
            (int)(255.1f * a));
    }

    static RgbaColor FromPackedInt(uint32_t rgba) {
        return From(
            (int)((rgba)       & 0xff),
            (int)((rgba >> 8)  & 0xff),
            (int)((rgba >> 16) & 0xff),
            (int)(255 - ((rgba >> 24) & 0xff)));
    }

    static RgbaColor FromPackedIntBGRA(uint32_t bgra) {
        return From(
            (int)((bgra >> 16) & 0xff),
            (int)((bgra >> 8)  & 0xff),
            (int)((bgra)       & 0xff),
            (int)(255 - ((bgra >> 24) & 0xff)));
    }
};

struct RgbaColorCompare {
    bool operator()(RgbaColor a, RgbaColor b) const {
        return a.ToARGB32() < b.ToARGB32();
    }
};

}//ss