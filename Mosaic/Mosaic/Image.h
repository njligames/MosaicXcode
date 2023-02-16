//
//  Image.h
//  JLIGameEngineTest
//
//  Created by James Folk on 1/8/15.
//  Copyright (c) 2015 James Folk. All rights reserved.
//

#ifndef __JLIGameEngineTest__Image__
#define __JLIGameEngineTest__Image__

#include "glm.hpp"
#include <string>
#include <vector>

namespace NJLIC {

    /**
     *  <#Description#>
     */
    class Image {
      public:
        Image();
        virtual ~Image();

        Image(const Image &);

        Image &operator=(const Image &);

      public:
        /**
         *  <#Description#>
         *
         *  @param imageSource <#imageSource description#>
         *  @param position    <#position description#>
         *  @param dimensions  <#dimensions description#>
         *
         *  @return <#return value description#>
         */
        //        static Image *createSubImage(const Image &imageSource,
        //                                     const glm::vec2 &position);

        // TODO: fill in specific methods for Image
        /**
         *  <#Description#>
         *
         *  @param position <#position description#>
         *  @param color    <#color description#>
         *
         *  @return <#return value description#>
         */
        bool setPixel(const glm::vec2 &position, const glm::vec4 &color);
        /**
         *  <#Description#>
         *
         *  @param position    <#position description#>
         *  @param imageSource <#imageSource description#>
         *
         *  @return <#return value description#>
         */
        bool setPixel(const glm::vec2 &position, const Image &imageSource);
        /**
         *  <#Description#>
         *
         *  @param position <#position description#>
         *  @param pixel    <#pixel description#>
         *
         *  @return <#return value description#>
         */
        bool getPixel(const glm::vec2 &position, glm::vec4 &pixel) const;
        glm::vec4 getPixel(const glm::vec2 &position) const;
        bool getPixel(const glm::vec2 &position, unsigned int &pixel) const;

        /**
         *  <#Description#>
         *
         *  @param position   <#position description#>
         *  @param dimensions <#dimensions description#>
         *  @param color      <#color description#>
         *
         *  @return <#return value description#>
         */
        //        bool setPixels(const glm::vec2 &position, const glm::vec2
        //        &dimensions,
        //                       const glm::vec4 &color);

        /**
         *  <#Description#>
         *
         *  @param position    <#position description#>
         *  @param dimensions  <#dimensions description#>
         *  @param sourceImage <#sourceImage description#>
         *  @param 0           <#0 description#>
         *  @param 0           <#0 description#>
         *
         *  @return <#return value description#>
         */
        //        bool setPixels(const glm::vec2 &position, const glm::vec2
        //        &dimensions,
        //                       const Image &sourceImage,
        //                       const glm::vec2 &sourceImageOffset =
        //                       glm::vec2(0, 0));

        bool setPixels(const glm::vec2 &toPosition, const Image &fromImage,
                       const glm::vec2 &fromOffsetPosition = glm::vec2(0, 0));

        /**
         *  <#Description#>
         *
         *  @param position         <#position description#>
         *  @param dimension        <#dimension description#>
         *  @param destinationImage <#destinationImage description#>
         *  @param 0                <#0 description#>
         *  @param 0                <#0 description#>
         *
         *  @return <#return value description#>
         */
        //        bool getPixels(
        //            const glm::vec2 &position, const glm::vec2 &dimension,
        //            Image &destinationImage,
        //            const glm::vec2 &destinationImageOffset = glm::vec2(0, 0))
        //            const;

        /**
         *  <#Description#>
         *
         *  @return <#return value description#>
         */
        int getNumberOfComponents() const;
        /**
         *  <#Description#>
         *
         *  @return <#return value description#>
         */
        int getWidth() const;
        /**
         *  <#Description#>
         *
         *  @return <#return value description#>
         */
        int getHeight() const;
        /**
         *  <#Description#>
         *
         *  @return <#return value description#>
         */
        const char *getFilename() const;
        /**
         *  <#Description#>
         *
         *  @return <#return value description#>
         */
        int getBytesPerPixel() const;

        int *getDataPtr() const;
        long getDataSize() const;

        /**
         *  <#Description#>
         *
         *  @return <#return value description#>
         */
        //    const int* getDataRaw() const;

        //        bool isPvr() const;
        //        bool isCompressed() const;

        //        njliImageType getImageType() const;

        /**
         *  <#Description#>
         *
         *  @return <#return value description#>
         */
        //    int* getDataRaw();
        /**
         *  <#Description#>
         *
         *  @param width              <#width description#>
         *  @param height             <#height description#>
         *  @param numberOfComponents <#numberOfComponents description#>
         *  @param const              <#const description#>
         *  @param 1.0f               <#1.0f description#>
         *  @param 1.0f               <#1.0f description#>
         *  @param 1.0f               <#1.0f description#>
         *  @param 0.0f               <#0.0f description#>
         */
        //    void setDataRaw(int width, int height, int numberOfComponents,
        //    const int* const, const glm::vec4& fillColor =
        //    glm::vec4(1.0f, 1.0f, 1.0f, 0.0f)); void copy(const char
        //    *fileName);

        /**
         *  <#Description#>
         *
         *  @param width              <#width description#>
         *  @param height             <#height description#>
         *  @param numberOfComponents <#numberOfComponents description#>
         *  @param 1.0f               <#1.0f description#>
         *  @param 1.0f               <#1.0f description#>
         *  @param 1.0f               <#1.0f description#>
         *  @param 1.0f               <#1.0f description#>
         */
        void generate(int width, int height, int numberOfComponents,
                      const glm::vec4 &fillColor = glm::vec4(1.0f, 0.0f, 0.0f,
                                                             1.0f));
        /**
         *  <#Description#>
         *
         *  @return <#return value description#>
         */
        //    u64 getDataRawLength() const;

        /**
         *  <#Description#>
         *
         *  @param from  <#from description#>
         *  @param to    <#to description#>
         *  @param color <#color description#>
         */
        void drawLine(const glm::vec2 &from, const glm::vec2 &to,
                      const glm::vec4 &color);

        /**
         *  <#Description#>
         *
         *  @param from        <#from description#>
         *  @param to          <#to description#>
         *  @param imageSource <#imageSource description#>
         */
        void drawLine(const glm::vec2 &from, const glm::vec2 &to,
                      const Image &imageSource);

        /**
         *  <#Description#>
         */
        void blur();
        /**
         *  <#Description#>
         *
         *  @param image <#image description#>
         *
         *  @return <#return value description#>
         */
        bool setAlpha(const Image &image);
        /**
         *  <#Description#>
         *
         *  @param alpha <#alpha description#>
         *
         *  @return <#return value description#>
         */
        bool setAlpha(float alpha);
        /**
         *  <#Description#>
         *
         *  @return <#return value description#>
         */
        bool hasAlpha() const;
        /**
         *  <#Description#>
         */
        void preMultiplyAlpha();
        /**
         *  <#Description#>
         */
        void flip();

        void rotate();
        
        void resize(int width, int height);
        Image resize(int width, int height)const;

        //    void setPVRData(int *pvrData, unsigned long dataSize, const char
        //    *fileName);
        //    int *getCompressedData()const;

        // get the pixel row
        std::vector<glm::vec4> operator[](int row);

        bool isWidthHeightPowerOf2() const;

        bool copyData(void *dataPtr, int width, int height, int components,
                      const std::string &filename);

        void setId(int _id) { mId = _id; }
        int getId() const { return mId; }

      protected:
        static void copyPixelRow(const Image &src, unsigned int srcY,
                                 Image &dest, unsigned int destY,
                                 unsigned int srcXOffset = 0,
                                 unsigned int destXOffset = 0);
        static bool setPixel(Image &image, const glm::vec2 &pos,
                             const glm::vec4 &color);
        static bool getPixel(const Image &image, const glm::vec2 &pos,
                             glm::vec4 &color);

        //        bool copyData(const WorldResourceLoader::ImageFileData
        //        *fileData);
        int getClosestValidGLDim(const int dim) const;
        //    void setDataRawFromWorldResourceLoader(int*, int x, int y, int
        //    numberOfComponents, const char* filename);
        //    bool isInWorldResourceLoader() const;

        //        void setPixelRow(unsigned char *data, int row, int width, int
        //        xOffset = 0); void getPixelRow(unsigned char *data, int row,
        //        int width);

        unsigned char *createFillRow_createsmemory(
            int xOffset, int fillWidth,
            const glm::vec4 &fillColor = glm::vec4(1.0f, 1.0f, 1.0f, 1.0f));

      private:
        // An output image with N components has the following components
        // interleaved
        // in this order in each pixel:
        //
        //     N=#comp     components
        //       1           grey
        //       2           grey, alpha
        //       3           red, green, blue
        //       4           red, green, blue, alpha

        unsigned char *m_RawData;

        long m_RawDataSize;

        int m_Width;
        int m_Height;
        int m_Componenents;

        std::string m_Filename;

        int mId;
        //    bool m_hasAlpha;

        //    bool m_IsInWorldResourceLoader;
        //    unsigned long m_pvrDataSize;
    };
} // namespace NJLIC

#endif /* defined(__JLIGameEngineTest__Image__) */
