#pragma once
#include "fermat2.h"
#include <stdio.h>
#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <fstream>
#include "file_system.h"


using namespace std;
class Image_Conv
{
public:
	Image_Conv();
	~Image_Conv();

	void addOutlierOnImage(cv::Mat& srcImage);
	void filterImage(cv::Mat& srcImage);
	void KMeansSegImage(cv::Mat& srcImage);
	void doExpandAndErosionImage(cv::Mat& srcBuildingImage);
	void doLaplaceKernel();
	void doCannyKernel();
	void doSharpenKernel();
	void makeBgBlack(int i, int j, int cPointB, int cPointG, int cPointR);
	void makeBgBlackAndForeWhite(int i, int j, int cPointB, int cPointG, int cPointR);
	void makeObjectDepartBgBlack(int i, int j, int c_pixel_b, int c_pixel_g, int c_pixel_r);

	void doSplitBuildAndGrasslake(cv::Mat srcGrassAndLakeImage, cv::Mat srcBuildingImage, int i, int j, int cPointB,
	                              int cPointG, int cPointR);
	void doMergeTwoImages(cv::Mat srcGrassAndLakeImage, cv::Mat srcBuildingImage, float weigthOne = 1.0, float weightTwo = 1.0);
	int getContoursByCplus(cv::Mat& srcImage, double minarea, double whRatio);
	void reserveInterestRegion(int i, int j, int cPointB, int cPointG, int cPointR);
	void init();

	//grid map
	void doGridMap(cv::Mat& srcImage);
	void biInterpolation(cv::Mat& src);
private:

	string baseFileName;
	cv::Mat srcMapImage;

	cv::Mat outExpandImage;
	cv::Mat expandElement;
	cv::Mat erosionElement;
	cv::Mat outErosionImage; //½øÐÐ¸¯Ê´²Ù×÷ 
	cv::Mat dstSharpenImage;
	cv::Mat dstLapImage;
	cv::Mat grayImage , blurImage , canny_detection;
	cv::Mat dstCannyKernel;
};

