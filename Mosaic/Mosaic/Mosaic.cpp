//
//  Mosaic.cpp
//  Mosaic
//
//  Created by James Folk on 2/13/23.
//

#include "Mosaic.hpp"

#include <iostream>
#include <vector>
#include <string>
#include <thread>

#include "Color.h"

#include "nlohmann/json.hpp"
using json = nlohmann::json;
using namespace NJLIC;




using namespace std;

namespace Mosaic {
class Generator;
}

static Image tileImage(const Image& targetImage, const vector<Image>& images, int tileSize, int numThreads);
static Image generateMosaic(const Image& targetImage, const vector<Image>& images, int tileSize, int numThreads, Mosaic::MosaicMap&);

static int getMaxThreads() {
    int maxThreads = std::thread::hardware_concurrency();
    if (maxThreads == 0) {
        std::cout << "Unable to determine the number of hardware threads. Using default value of 1." << std::endl;
        maxThreads = 1;
    }
    return maxThreads;
}

namespace Mosaic {

static Image sImage;

bool Generator::generate(const ImageData& img) {
    
    int numThreads = getMaxThreads();
    
    Image targetImage;
    targetImage.copyData(img.data, img.width, img.height, img.components, img.id);
    
    vector<Image> resizedImages;
    for(auto iter = tileImages_.begin(); iter != tileImages_.end(); iter++) {
        ImageData tileImgData = *iter;
        Image tileImg;
        
        tileImg.copyData(tileImgData.data, tileImgData.width, tileImgData.height, tileImgData.components, tileImgData.id);
        
        tileImg.resize(getTileSize(), getTileSize());
        resizedImages.push_back(tileImg);
    }
    
    mosaicMap_.clear();
    sImage = generateMosaic(targetImage, resizedImages, tileSize_, numThreads, mosaicMap_);
    
    mosaicImage_.data = (unsigned char*)sImage.getDataPtr();
    mosaicImage_.width = sImage.getWidth();
    mosaicImage_.height = sImage.getHeight();
    mosaicImage_.components = sImage.getNumberOfComponents();
    mosaicImage_.id = sImage.getFilename();
    
    return true;
}

// Convert a MosaicMap to a JSON object
string MosaicMapToJson(const MosaicMap& mosaicMap) {
    json j;

    // iterate through the MosaicMap and add each pair to the JSON object
    for (const auto& pair : mosaicMap) {
        Indices indices = pair.first;
        string filename = pair.second;

        // create a JSON object for the pair and add it to the JSON array
        j.push_back({
            {"x", indices.first},
            {"y", indices.second},
            {"filename", filename}
        });
    }

    return j.dump();
}

string Generator::getMosaicMap()const {
    return MosaicMapToJson(mosaicMap_);
}

void Generator::addTileImage(const ImageData &image) {
    tileImages_.push_back(image);
}

void Generator::clear() {
    tileImages_.clear();
}

//Image Generator::tile(const Image& targetImage, const vector<Image>& images)const {
//    int numThreads = getMaxThreads();
//    auto t = tileImage(targetImage, images, tileSize_, numThreads);
//    Image imgData(t);
//    return imgData;
//}

}

// Define a global mutex to protect the image map
static mutex imageMapMutex;

// Define a function to calculate the similarity between two images
static double calculateSimilarity(const Image& target, int image1OffsetX, int image1OffsetY, int image1SizeX, int image1SizeY,
                               const Image& image) {
    int sum = 0;
    for (int y = 0; y < image.getHeight(); y++) {
        for (int x = 0; x < image.getWidth(); x++) {
            
            glm::vec4 pixel1;
            image.getPixel(glm::vec2(x, y), pixel1);
            Color color1;
            color1.setRGB(pixel1);

            glm::vec4 pixel2;
            target.getPixel(glm::vec2(x + image1OffsetX, y + image1OffsetY), pixel2);
            Color color2;
            color2.setRGB(pixel2);

            int diff = color1.distance(color2);
            
            sum += diff * diff;
        }
    }

    double mse = (double) sum / (image.getWidth() * image.getHeight() * image.getNumberOfComponents());
    double rmse = sqrt(mse);
    return (int) rmse;
}

static Image tileImage(const Image& targetImage, const vector<Image>& images, int tileSize, int numThreads) {
    Image img;
    img.generate(targetImage.getWidth(), targetImage.getHeight(), targetImage.getNumberOfComponents());
    int targetWidth = targetImage.getWidth();
    int targetHeight = targetImage.getHeight();
    int tileCols = targetWidth / tileSize;
    int tileRows = targetHeight / tileSize;
    
    
    
    // Define a function to find the best match for each tile in parallel
    auto findBestMatch = [&](int threadIndex) {
        for (int j = threadIndex; j < tileCols * tileRows; j += numThreads) {

            // Get the tile coordinates
            int tileRow = j / tileCols;
            int tileCol = j % tileCols;
            int tileX = tileCol * tileSize;
            int tileY = tileRow * tileSize;
            
            int bestMatchIndex = j % images.size();
            img.setPixels(glm::vec2(tileX, tileY), images[bestMatchIndex]);
        }
    };

    // Find the best match for each tile in parallel
    vector<thread> threads;
    for (int i = 0; i < numThreads; i++) {
        threads.push_back(thread(findBestMatch, i));
    }
    for (auto& t : threads) {
        t.join();
    }
    
    return img;
}

// Define a function to generate a mosaic image
static Image generateMosaic(const Image& targetImage, const vector<Image>& images, int tileSize, int numThreads, Mosaic::MosaicMap &mmap) {
    int targetWidth = targetImage.getWidth();
    int targetHeight = targetImage.getHeight();
    int tileCols = targetWidth / tileSize;
    int tileRows = targetHeight / tileSize;

    Image mosaicPixels;
    mosaicPixels.copyData(targetImage.getDataPtr(), targetImage.getWidth(), targetImage.getHeight(), targetImage.getNumberOfComponents(), targetImage.getFilename());
    
    // Create a vector to hold the similarity scores
    vector<vector<double>> similarityScores(images.size(), vector<double>(tileCols * tileRows));

    // Define a function to calculate the similarity scores in parallel
    auto calculateSimilarityScores = [&](int threadIndex) {
        for (int i = threadIndex; i < images.size(); i += numThreads) {
            for (int j = 0; j < tileCols * tileRows; j++) {
                
                int tileRow = j / tileCols;
                int tileCol = j % tileCols;
                int tileX = tileCol * tileSize;
                int tileY = tileRow * tileSize;
                
                double similarity = calculateSimilarity(targetImage, tileX, tileY, tileSize, tileSize, images[i]);
                similarityScores[i][j] = similarity;
            }
        }
    };

    // Calculate the similarity scores in parallel
    vector<thread> threads;
    for (int i = 0; i < numThreads; i++) {
        threads.push_back(thread(calculateSimilarityScores, i));
    }
    for (auto& t : threads) {
        t.join();
    }

    // Define a function to find the best match for each tile in parallel
    auto findBestMatch = [&](int threadIndex) {
        for (int j = threadIndex; j < tileCols * tileRows; j += numThreads) {

            // Get the tile coordinates
            int tileRow = j / tileCols;
            int tileCol = j % tileCols;
            int tileX = tileCol * tileSize;
            int tileY = tileRow * tileSize;

            // Find the best matching image for the tile
            int bestMatchIndex = 0;
            int bestMatchScore = similarityScores[0][j];
            for (int i = 1; i < images.size(); i++) {
                int similarity = similarityScores[i][j];
                if (similarity < bestMatchScore) {
                    bestMatchIndex = i;
                    bestMatchScore = similarity;
                }
            }
            
            mosaicPixels.setPixels(glm::vec2(tileX, tileY), images[bestMatchIndex]);
            
            imageMapMutex.lock();
            mmap.insert(Mosaic::MosaicMapPair(Mosaic::Indices(tileX / tileSize, tileY / tileSize), images[bestMatchIndex].getFilename()));
            imageMapMutex.unlock();
        }
    };

    // Find the best match for each tile in parallel
    threads.clear();
    for (int i = 0; i < numThreads; i++) {
        threads.push_back(thread(findBestMatch, i));
    }
    for (auto& t : threads) {
        t.join();
    }
    
    return mosaicPixels;
}
