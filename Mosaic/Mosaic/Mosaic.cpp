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


#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

using namespace std;

namespace Mosaic {
class ImageData;
class Generator;
}

static Mosaic::ImageData loadImage(const string& filename);
static Image tileImage(const Image& targetImage, const vector<Image>& images, int tileSize, int numThreads);
static Image generateMosaic(const Image& targetImage, const vector<Image>& images, int tileSize, int numThreads, Mosaic::Generator::MosaicMap&);
static vector<vector<int>> generateMosaic(const Mosaic::ImageData& targetImage, const vector<Mosaic::ImageData>& images, int tileSize, int numThreads);

static int getMaxThreads() {
    int maxThreads = std::thread::hardware_concurrency();
    if (maxThreads == 0) {
        std::cout << "Unable to determine the number of hardware threads. Using default value of 1." << std::endl;
        maxThreads = 1;
    }
    return maxThreads;
}

namespace Mosaic {















//static Image load(const string &filename);
//static void write(const string &filename, const Image &img);
// Define a function to load an image
Image ImageFileLoader::load(const string& filename) {
    // Check if the image is already loaded in the map
//    imageMapMutex.lock();
//    auto it = imageMap.find(filename);
//    if (it != imageMap.end()) {
//        Image image = it->second;
//        imageMapMutex.unlock();
//        return image;
//    }
//    imageMapMutex.unlock();
    
    int width, height, channels;
    unsigned char *data = stbi_load(filename.c_str(), &width, &height, &channels, 0);
    if (!data) {
        throw runtime_error("Failed to load image: " + filename);
    }

    Image image;
    image.copyData(data, width, height, channels, filename);

    stbi_image_free(data);
    
    // Add the loaded image to the image map
//    imageMapMutex.lock();
//    imageMap[filename] = image;
//    imageMapMutex.unlock();

    return image;
}

void ImageFileLoader::write(const string &filename, const Image &img) {
    stbi_write_png(filename.c_str(),
                   img.getWidth(),
                   img.getHeight(),
                   img.getNumberOfComponents(),
                   img.getDataPtr(),
                   img.getNumberOfComponents() * img.getWidth());
}

ImageData::ImageData(const std::string &path) { *this = loadImage(path);}
ImageData::ImageData(const vector<vector<int>> &pixels) {
    this->width = pixels[0].size();
    this->height = pixels.size();
    this->pixels = pixels;
}
void ImageData::write(const string &filename) {
    int channels = 3; // assuming RGB color channels
    stbi_write_png(filename.c_str(), width, height, channels, &pixels[0][0], channels*width);
}

//Generator::MosaicMap Generator::generateMap(const Image& targetImage, const vector<Image>& images)const {
//    int numThreads = getMaxThreads();
//    MosaicMap mmap;
//    
//    generateMosaic(targetImage, images, tileSize_, numThreads, mmap);
//    
//    return mmap;
//}

bool Generator::generate(const Image& targetImage, const vector<Image>& images) {
    
    int numThreads = getMaxThreads();
//    MosaicMap mmap;
    
//    typedef map<string, Image> ImageMap;
//    typedef pair<string, Image> ImagePairMap;
//    ImageMap imagemap;
    
//    int maxWidth = std::numeric_limits<int>::min();
//    int maxHeight = std::numeric_limits<int>::min();
    vector<Image> resizedImages;
    for(auto iter = images.begin(); iter != images.end(); iter++) {
//        if(iter->getWidth() > maxWidth)maxWidth = iter->getWidth();
//        if(iter->getHeight() > maxHeight)maxHeight = iter->getHeight();
        resizedImages.push_back(iter->resize(getTileSize(), getTileSize()));
//        imagemap.insert(ImagePairMap(iter->getFilename(), *iter));
    }
    
//    int targetWidth = targetImage.getWidth();
//    int targetHeight = targetImage.getHeight();
//    int tileCols = targetWidth / getTileSize();
//    int tileRows = targetHeight / getTileSize();
    
    mosaicMap_.clear();
    mosaicImage_ = generateMosaic(targetImage, resizedImages, tileSize_, numThreads, mosaicMap_);
//    Image imgData(t);
    
//    Image imgData(targetImage);
//    imgData.resize(maxWidth*tileCols, maxWidth*tileRows);
//    for(auto iter = mmap.begin(); iter != mmap.end(); iter++) {
//        int tileX = iter->first.first;
//        int tileY = iter->first.second;
//        string filename = iter->second;
//        const Image &img = imagemap[filename];
//
//        imgData.setPixels(glm::vec2(tileX, tileY), img);
//    }
    
    return true;
    
}

Image Generator::tile(const Image& targetImage, const vector<Image>& images)const {
    int numThreads = getMaxThreads();
    auto t = tileImage(targetImage, images, tileSize_, numThreads);
    Image imgData(t);
    return imgData;
}

}

/// This function takes in the filename of an image file and returns a vector<vector<int>> that contains the pixel values of the image. The function also sets the width, height, and channels variables to the corresponding values of the loaded image.

/// This function first calls stbi_load to load the image file and store the pixel values in a linear array of unsigned chars. It then creates a 2D vector of int values with the same dimensions as the loaded image and loops through the pixels to store the pixel values in the vector. Finally, the function frees the memory allocated by stbi_load and returns the vector of pixel values.
/// - Parameters:
///   - filename: filename description
///   - width: <#width description#>
///   - height: <#height description#>
///   - channels: <#channels description#>
static std::vector<std::vector<int>> loadPixels(const char* filename, int& width, int& height, int& channels) {
    unsigned char* image = stbi_load(filename, &width, &height, &channels, STBI_rgb);

    if (!image) {
        throw std::runtime_error("Failed to load image file.");
    }

    std::vector<std::vector<int>> pixels(height, std::vector<int>(width * channels));

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            for (int c = 0; c < channels; c++) {
                pixels[y][x * channels + c] = image[y * width * channels + x * channels + c];
            }
        }
    }

    stbi_image_free(image);

    return pixels;
}






// Define a global mutex to protect the image map
static mutex imageMapMutex;

// Define a global image map to store already loaded images
static unordered_map<string, Mosaic::ImageData> imageMap;

// Define a function to load an image
static Mosaic::ImageData loadImage(const string& filename) {
    // Check if the image is already loaded in the map
    imageMapMutex.lock();
    auto it = imageMap.find(filename);
    if (it != imageMap.end()) {
        Mosaic::ImageData image = it->second;
        imageMapMutex.unlock();
        return image;
    }
    imageMapMutex.unlock();
    
    
    
    
    
    int width, height, channels;
    unsigned char *data = stbi_load(filename.c_str(), &width, &height, &channels, 0);
    if (!data) {
        throw runtime_error("Failed to load image: " + filename);
    }

    vector<vector<int>> pixels(height, vector<int>(width));
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            pixels[i][j] = data[i * width + j];
        }
    }

    stbi_image_free(data);
    Mosaic::ImageData image(pixels);



//    // Load image from disk
//    vector<vector<int>> pixels;
//    int x,y,n;
//
//    unsigned char *data = stbi_load(filename.c_str(), &x, &y, &n, 0);
//
//    for(int i = 0; i < y; i++) {
//        vector<int> row;
//        for(int j = 0; j < x; j++) {
//            row.push_back(data[x + (i*y)]);
//        }
//        pixels.push_back(row);
//    }
//    stbi_image_free(data);
//
//    // Convert the loaded image data to the ImageData structure
//    Mosaic::ImageData image;
////    image.width = x;
////    image.height = y;
////    image.pixels = pixels;
//
//    image.setWidth(x);
//    image.setHeight(y);
//    image.setPixels(pixels);
    // Add the loaded image to the image map
    imageMapMutex.lock();
    imageMap[filename] = image;
    imageMapMutex.unlock();

    return image;
}

// Define a function to calculate the similarity between two images
static double calculateSimilarity(const Image& target, int image1OffsetX, int image1OffsetY, int image1SizeX, int image1SizeY,
                               const Image& image) {
    int sum = 0;
//    double totaldiff = 0.0;
    for (int y = 0; y < image.getHeight(); y++) {
        for (int x = 0; x < image.getWidth(); x++) {
//            unsigned int c1, c2;
//            image.getPixel(glm::vec2(x, y), c1);
//            target.getPixel(glm::vec2(x + image1OffsetX, y + image1OffsetY), c2);
            
            
            
//            int diff = image1.getPixels()[y][x] - image2.getPixels()[y][x];
            
//            {
//                glm::vec4 pixel1;
//                image.getPixel(glm::vec2(x, y), pixel1);
//
//                glm::vec4 pixel2;
//                target.getPixel(glm::vec2(x + image1OffsetX, y + image1OffsetY), pixel2);
//
//                int rFirst = pixel1.r * 255;
//                int gFirst = pixel1.g * 255;
//                int bFirst = pixel1.b * 255;
//                int rSecond = pixel2.r * 255;
//                int gSecond = pixel2.g * 255;
//                int bSecond = pixel2.b * 255;
//                totaldiff += std::abs( rFirst - rSecond ) / 255.0 ;
//                totaldiff += std::abs( gFirst - gSecond ) / 255.0 ;
//                totaldiff += std::abs( bFirst -bSecond ) / 255.0 ;
//
//            }
            {
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
//            {
//                glm::vec4 pixel1;
//                image.getPixel(glm::vec2(x, y), pixel1);
//                Color color1;
//                color1.setRGB(pixel1);
//                int hue1 = color1.getHSV().x ;
//
//                glm::vec4 pixel2;
//                target.getPixel(glm::vec2(x + image1OffsetX, y + image1OffsetY), pixel2);
//                Color color2;
//                color2.setRGB(pixel2);
//                int hue2 = color2.getHSV().x ;
//
//                float avgHue = (hue1 + hue2) / 2.f;
//
//                int diff = abs(hue1 - avgHue);
//
//                sum += diff * diff;
//            }
            
        }
    }

    double mse = (double) sum / (image.getWidth() * image.getHeight() * image.getNumberOfComponents());
    double rmse = sqrt(mse);
    return (int) rmse;
//    return (totaldiff * 100)  / (image.getWidth() * image.getHeight() * image.getNumberOfComponents());
}


// Define a function to calculate the similarity between two images
static int calculateSimilarity(const Mosaic::ImageData& image1, const Mosaic::ImageData& image2) {
    if (image1.getWidth() != image2.getWidth() || image1.getHeight() != image2.getHeight()) {
        return INT_MAX;
    }

    int sum = 0;
    for (int y = 0; y < image1.getHeight(); y++) {
        for (int x = 0; x < image1.getWidth(); x++) {
            int diff = image1.getPixels()[y][x] - image2.getPixels()[y][x];
            sum += diff * diff;
        }
    }

    double mse = (double) sum / (image1.getWidth() * image1.getHeight());
    double rmse = sqrt(mse);
    return (int) rmse;
}

// Define a function to calculate the similarity between two images
static int calculateSimilarity2(const Mosaic::ImageData& image1, const Mosaic::ImageData& image2) {
    if (image1.getWidth() != image2.getWidth() || image1.getHeight() != image2.getHeight()) {
        return INT_MAX;
    }

    int sum = 0;
    for (int y = 0; y < image1.getHeight(); y++) {
        for (int x = 0; x < image1.getWidth(); x++) {
            int diff = abs(image1.getPixels()[y][x] - image2.getPixels()[y][x]);
            sum += diff;
        }
    }

    return sum;
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
static Image generateMosaic(const Image& targetImage, const vector<Image>& images, int tileSize, int numThreads, Mosaic::Generator::MosaicMap &mmap) {
    int targetWidth = targetImage.getWidth();
    int targetHeight = targetImage.getHeight();
    int tileCols = targetWidth / tileSize;
    int tileRows = targetHeight / tileSize;

    // Create a vector to hold the mosaic image data
//    vector<vector<int>> mosaicPixels(targetHeight, vector<int>(targetWidth));
    Image mosaicPixels;
//    mosaicPixels.generate(targetImage.getWidth(), targetImage.getHeight(), targetImage.getNumberOfComponents());
    mosaicPixels.copyData(targetImage.getDataPtr(), targetImage.getWidth(), targetImage.getHeight(), targetImage.getNumberOfComponents(), targetImage.getFilename());
    

//    // Create a vector to hold the loaded images
//    vector<ImageData> images;
//    for (const string& filename : imageFilenames) {
//        images.push_back(loadImage(filename));
//    }

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

//                const Mosaic::ImageData& tileImage = images[i].getSubImage(tileX, tileY, tileSize, tileSize);
//                int similarity = calculateSimilarity(targetImage, tileImage);
//                similarityScores[i][j] = similarity;


            }
        }
    };
    
    
    
//    for(int i = 0; i < images.size(); i++) {
//        for(int j = 0; j < tileCols * tileRows; j++) {
//            int tileRow = j / tileCols;
//            int tileCol = j % tileCols;
//            int tileX = tileCol * tileSize;
//            int tileY = tileRow * tileSize;
//
//            int similarity = calculateSimilarity(targetImage, tileX, tileY, tileSize, tileSize, images[i]);
//            cout << similarity << endl;
//            similarityScores[i][j] = similarity;
//        }
//    }

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

            // Copy the best matching image data into the mosaic image data
//            const Mosaic::ImageData& bestMatchImage = images[bestMatchIndex];
//            for (int y = 0; y < tileSize; y++) {
//                for (int x = 0; x < tileSize; x++) {
//                    int sourceX = x;
//                    int sourceY = y;
//                    int targetX = tileX + x;
//                    int targetY = tileY + y;
//                    mosaicPixels[targetY][targetX] = bestMatchImage.getPixels()[sourceY][sourceX];
//                }
//            }
            mosaicPixels.setPixels(glm::vec2(tileX, tileY), images[bestMatchIndex]);
            
            
            
            
            imageMapMutex.lock();
            mmap.insert(Mosaic::Generator::MosaicMapPair(Mosaic::Generator::Indices(tileX / tileSize, tileY / tileSize), images[bestMatchIndex].getFilename()));
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

    // Output the mosaic image data
    // ...
    return mosaicPixels;
}

















// Define a function to generate a mosaic image
static vector<vector<int>> generateMosaic(const Mosaic::ImageData& targetImage, const vector<Mosaic::ImageData>& images, int tileSize, int numThreads) {
    int targetWidth = targetImage.getWidth();
    int targetHeight = targetImage.getHeight();
    int tileCols = targetWidth / tileSize;
    int tileRows = targetHeight / tileSize;

    // Create a vector to hold the mosaic image data
    vector<vector<int>> mosaicPixels(targetHeight, vector<int>(targetWidth));

//    // Create a vector to hold the loaded images
//    vector<ImageData> images;
//    for (const string& filename : imageFilenames) {
//        images.push_back(loadImage(filename));
//    }

    // Create a vector to hold the similarity scores
    vector<vector<int>> similarityScores(images.size(), vector<int>(tileCols * tileRows));

    // Define a function to calculate the similarity scores in parallel
    auto calculateSimilarityScores = [&](int threadIndex) {
        for (int i = threadIndex; i < images.size(); i += numThreads) {
            for (int j = 0; j < tileCols * tileRows; j++) {
                int tileRow = j / tileCols;
                int tileCol = j % tileCols;
                int tileX = tileCol * tileSize;
                int tileY = tileRow * tileSize;
                
                const Mosaic::ImageData& tileImage = images[i].getSubImage(tileX, tileY, tileSize, tileSize);
                int similarity = calculateSimilarity(targetImage, tileImage);
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
                if (similarity > bestMatchScore) {
                    bestMatchIndex = i;
                    bestMatchScore = similarity;
                }
            }

            // Copy the best matching image data into the mosaic image data
            const Mosaic::ImageData& bestMatchImage = images[bestMatchIndex];
            for (int y = 0; y < tileSize; y++) {
                for (int x = 0; x < tileSize; x++) {
                    int sourceX = x;
                    int sourceY = y;
                    int targetX = tileX + x;
                    int targetY = tileY + y;
                    mosaicPixels[targetY][targetX] = bestMatchImage.getPixels()[sourceY][sourceX];
                }
            }
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

    // Output the mosaic image data
    // ...
    return mosaicPixels;
}


//void Generator::generate(const string& targetImageFilename, const vector<string>& imageFilenames, const string &outputFilename) {
////        // Load the target image data
////        ImageData targetImage = loadImage("target.jpg");
////
////        // Load the image filenames
////        vector<string> imageFilenames = {"image1.jpg", "image2.jpg", "image3.jpg"};
////
////        // Generate the mosaic image
////        int tileSize = 50;
////        int numThreads = getMaxThreads();
////        generateMosaic(targetImage, imageFilenames, tileSize, numThreads);
//}
