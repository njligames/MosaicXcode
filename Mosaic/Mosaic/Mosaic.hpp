//
//  Mosaic.hpp
//  Mosaic
//
//  Created by James Folk on 2/13/23.
//

#ifndef Mosaic_hpp
#define Mosaic_hpp

#include <string>
#include <vector>
#include <mutex>
#include <unordered_map>
#include "Image.h"
#include <map>

using namespace std;
using namespace NJLIC;

namespace Mosaic {
class ImageFileLoader {
//    // Define a global mutex to protect the image map
//    static mutex imageMapMutex;
//
//    // Define a global image map to store already loaded images
//    static unordered_map<string, Image> imageMap;
public:
    static Image load(const string &filename);
    static void write(const string &filename, const Image &img);
};

// Define a structure to hold image data
class ImageData {
    int width;
    int height;
    vector<vector<int>> pixels;
public:
    // Default constructor
    ImageData() : width(0), height(0) {}

    // Constructor with width, height, and pixels parameters
    ImageData(int w, int h, const vector<vector<int>>& px) : width(w), height(h), pixels(px) {}

    // Copy constructor
    ImageData(const ImageData& other) : width(other.width), height(other.height), pixels(other.pixels) {}
    
    ImageData(const std::string &path);
    
    ImageData(const vector<vector<int>> &pixels);

    // Operator= function
    ImageData& operator=(const ImageData& other) {
        if (this != &other) {
            width = other.width;
            height = other.height;
            pixels = other.pixels;
        }
        return *this;
    }

    // Getter methods
    int getWidth() const { return width; }
    int getHeight() const { return height; }
    const vector<vector<int>>& getPixels() const { return pixels; }
    
    // Setter methods
    void setWidth(const int w) { width = w; }
    void setHeight(const int h) { height = h; }
    void setPixels(const vector<vector<int>> &p) { pixels = p; }
    
    ImageData getSubImage(int xoffset, int yoffset, int tileWidth, int tileHeight) const {

        vector<vector<int>> subPixels(tileHeight, vector<int>(tileWidth, 0));
        for (int i = 0; i < tileHeight; i++) {
            for (int j = 0; j < tileWidth; j++) {
                if ((xoffset < 0) && (xoffset + j > getWidth()) &&
                    (yoffset < 0) && (yoffset + i > getHeight())) {
                    subPixels[i][j] = getPixels()[yoffset + i][xoffset + j];
                }
            }
        }

        return ImageData(tileWidth, tileHeight, subPixels);
    }
    
    //write the image
    void write(const string &filename);
};

class Generator {
public:
    // Default constructor
    Generator() : tileSize_(16) {}
    
    // Constructor with tile size parameter
    Generator(int tileSize) : tileSize_(tileSize) {}
    
    // Getter method for tile size
    int getTileSize() const {
        return tileSize_;
    }
    
    // Setter method for tile size
    void setTileSize(int tileSize) {
        tileSize_ = tileSize;
    }
    
    typedef pair<int, int> Indices;
    typedef map<Indices, string> MosaicMap;
    typedef pair<Indices, string> MosaicMapPair;
    
//    MosaicMap generateMap(const Image& targetImage, const vector<Image>& images)const;
    
    bool generate(const Image& targetImage, const vector<Image>& images);
    Image tile(const Image& targetImage, const vector<Image>& images)const;
    
    const MosaicMap &getMosaicMap()const { return mosaicMap_;}
    const Image &getMosaicImage()const {return mosaicImage_;}
    
private:
    int tileSize_;
    MosaicMap mosaicMap_;
    Image mosaicImage_;
};
}

#endif /* Mosaic_hpp */
