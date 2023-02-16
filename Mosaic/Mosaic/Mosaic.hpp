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

namespace Mosaic {

struct ImageData {
    int width;
    int height;
    int components;
    string id;
    unsigned char *data;
};

typedef pair<int, int> Indices;
typedef map<Indices, string> MosaicMap;
typedef pair<Indices, string> MosaicMapPair;

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
    
    bool generate(const ImageData& targetImage);
    const ImageData &getMosaicImage()const { return mosaicImage_;}
    string getMosaicMap()const;
    
    void addTileImage(const ImageData &image);
    void clear();
    
    
private:
    int tileSize_;
    MosaicMap mosaicMap_;
    ImageData mosaicImage_;
    vector<ImageData> tileImages_;
};
}

#endif /* Mosaic_hpp */
