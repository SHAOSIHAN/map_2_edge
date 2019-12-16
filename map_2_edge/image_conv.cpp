#include "image_conv.h"
#include <opencv2/core/core.hpp>
using namespace cv;
Image_Conv::Image_Conv()
{

}

Image_Conv::~Image_Conv()
{

}

//Canny边缘检测相关变量
Mat g_cannyDetectedEdges;
int g_cannyLowThreshold = 1;//TrackBar位置参数  
RNG rng;

void Image_Conv::addOutlierOnImage(cv::Mat& srcImage)
{
	int height = srcImage.rows;
	int width = srcImage.cols;

	for (int i=0; i< height;i++)
	{
		srcImage.at<Vec3b>(i, 0)[0] = 0;
		srcImage.at<Vec3b>(i, 0)[1] = 0;
		srcImage.at<Vec3b>(i, 0)[2] = 0;

	}
	for (int i = 0; i < height; i++)
	{
		srcImage.at<Vec3b>(i, width-1)[0] = 0;
		srcImage.at<Vec3b>(i, width-1)[1] = 0;
		srcImage.at<Vec3b>(i, width-1)[2] = 0;

	}

	for (int i = 0; i < width; i++)
	{
		srcImage.at<Vec3b>(0, width - 1)[0] = 0;
		srcImage.at<Vec3b>(0, width - 1)[1] = 0;
		srcImage.at<Vec3b>(0, width - 1)[2] = 0;

	}
	for (int i = 0; i < width; i++)
	{
		srcImage.at<Vec3b>(height-1, i)[0] = 0;
		srcImage.at<Vec3b>(height-1, i)[1] = 0;
		srcImage.at<Vec3b>(height-1, i)[2] = 0;

	}

	imshow("add outlier",srcImage);
}

void Image_Conv::filterImage(cv::Mat& srcImage)
{
	//https://zhuanlan.zhihu.com/p/40325840
	Mat src, src2, dst1, dst2, dst3, dst4, maskResult;
	src = imread("./resources/tmp_saliency_2.bmp");
	if (src.empty()) {
		printf("could not load image...\n");
		return;
	}
	char input_title[] = "input image";
	char output_title[] = "blur image";
	namedWindow(input_title, CV_WINDOW_AUTOSIZE);
	namedWindow(output_title, CV_WINDOW_AUTOSIZE);
	imshow(input_title, src);
	blur(src, dst1, Size(5, 5), Point(-1, -1)); //均值滤波，窗口大小，中心点，borderType默认为4不用管
	imshow(output_title, dst1);

	GaussianBlur(src, dst2, Size(5, 5), 3, 3); //高斯模糊，高斯窗口大小x y坐标，σx σy正态分布情况
	imshow("GaussianBlur", dst2);
	medianBlur(src, dst3, 3); //中值滤波，窗口大小3*3
	imshow("medianBlur", dst3);
	bilateralFilter(src, dst4, 15, 100, 3); //高斯双边滤波，半径，像素差值，空间大小
	imshow("bilateralFilter", dst4);
	Mat kernel = (Mat_<int>(3, 3) << 0, -1, 0, -1, 5, -1, 0, -1, 0); //定义掩膜
	filter2D(dst4, maskResult, -1, kernel, Point(-1, -1), 0); //矩阵掩膜，提高对比度
	imshow("bilateralFilter and maskResult", maskResult);
	imwrite("./resources/bilateralFilterand maskResult.bmp", maskResult);
	waitKey(0);
}

void Image_Conv::KMeansSegImage(cv::Mat& srcImage)
{
	Mat img = imread("./resources/gaode-szu-big.png");
	namedWindow("Source Image", 0);
	imshow("Source Image", img);
	//生成一维采样点,包括所有图像像素点,注意采样点格式为32bit浮点数。 
	Mat samples(img.cols*img.rows, 1, CV_32FC3);
	//标记矩阵，32位整形 
	Mat labels(img.cols*img.rows, 1, CV_32SC1);
	uchar* p;
	int i, j, k = 0;
	for (i = 0; i < img.rows; i++)
	{
		p = img.ptr<uchar>(i);
		for (j = 0; j < img.cols; j++)
		{
			samples.at<Vec3f>(k, 0)[0] = float(p[j * 3]);
			samples.at<Vec3f>(k, 0)[1] = float(p[j * 3 + 1]);
			samples.at<Vec3f>(k, 0)[2] = float(p[j * 3 + 2]);
			k++;
		}
	}

	int clusterCount = 6;
	Mat centers(clusterCount, 1, samples.type());
	kmeans(samples, clusterCount, labels,
		TermCriteria(CV_TERMCRIT_EPS + CV_TERMCRIT_ITER, 10, 1.0),
		3, KMEANS_PP_CENTERS, centers);
	//我们已知有3个聚类，用不同的灰度层表示。 
	Mat img1(img.rows, img.cols, CV_8UC1);
	float step = 255 / (clusterCount - 1);
	k = 0;
	for (i = 0; i < img1.rows; i++)
	{
		p = img1.ptr<uchar>(i);
		for (j = 0; j < img1.cols; j++)
		{
			int tt = labels.at<int>(k, 0);
			k++;
			p[j] = 255 - tt*step;
		}
	}

	namedWindow("K-Means分割效果", 0);
	imshow("K-Means分割效果", img1);
	waitKey();
}

void Image_Conv::reserveInterestRegion(int i, int j, int cPointB, int cPointG, int cPointR)
{
	//将非建筑物、grass、lake的背景涂黑
	if (!(
		   (cPointB == 249 && cPointG == 219 && cPointR == 190) //lake
		|| (cPointB == 244 && cPointG == 247 && cPointR == 249) //building
		|| (cPointB == 204 && cPointG == 237 && cPointR == 208) //grass
		|| (cPointB == 255 && cPointG == 255 && cPointR == 255) //road
		))
	{
		srcMapImage.at<Vec3b>(i, j)[0] = 0;
		srcMapImage.at<Vec3b>(i, j)[1] = 0;
		srcMapImage.at<Vec3b>(i, j)[2] = 0;
	}	
}
void Image_Conv::doExpandAndErosionImage(cv::Mat& srcImage)
{	
	for (int k=0; k<1; k++)
	{
		expandElement = getStructuringElement(MORPH_RECT, Size(3, 3));
		dilate(srcImage, outExpandImage, expandElement);        //显示效果图
		imshow("建筑物【效果图】膨胀操作", outExpandImage);

		erosionElement = getStructuringElement(MORPH_RECT, Size(3, 3));
		erode(outExpandImage, outErosionImage, erosionElement);        //显示效果图 
		imshow("建筑物【效果图】腐蚀操作", outErosionImage);

		//erosionElement = getStructuringElement(MORPH_RECT, Size(3, 3));
		//erode(srcImage, outErosionImage, erosionElement);        //显示效果图 
		//imshow("建筑物【效果图】腐蚀操作", outErosionImage);
		//
		//expandElement = getStructuringElement(MORPH_RECT, Size(3, 3));
		//dilate(outErosionImage, outExpandImage, expandElement);        //显示效果图
		//imshow("建筑物【效果图】膨胀操作", outExpandImage);
	}
	imshow("建筑物【效果图】操作", outErosionImage);
	cvWaitKey(0);

}

void Image_Conv::doLaplaceKernel()
{
	imwrite("./resources/" + baseFileName +"_abstract.bmp", srcMapImage);
	imshow("原图", srcMapImage);
	cvWaitKey(0);

	Mat lap_kernel = (Mat_<double>(3, 3) <<
		1, 1, 1,
		1, -8, 1,
		1, 1, 1);
	filter2D(srcMapImage, dstLapImage, srcMapImage.depth(), lap_kernel);
	imwrite("./resources/" + baseFileName +"_abstract-lap.jpg",dstLapImage);
	imshow("lap卷积图", dstLapImage);

	Mat src = imread("./resources/" + baseFileName +"_abstract-lap.jpg");
	int m, n;
	for (m = 0; m < src.rows; m++)
		for (n = 0; n < src.cols; n++)
		{
			int cPixelB = src.at<Vec3b>(m, n)[0];
			int cPixelG = src.at<Vec3b>(m, n)[1];
			int cPixelR = src.at<Vec3b>(m, n)[2];
			if (cPixelB > 15 && cPixelG > 15 && cPixelR > 15)
			{
				dstLapImage.at<Vec3b>(m, n)[0] = 0;
				dstLapImage.at<Vec3b>(m, n)[1] = 128;
				dstLapImage.at<Vec3b>(m, n)[2] = 0;
			}

		}
	imwrite("./resources/" + baseFileName +"_abstract-lap_green.bmp", dstLapImage);
	imshow("da", dstLapImage);
}

void Image_Conv::doCannyKernel()
{
	cvtColor(srcMapImage,grayImage,COLOR_BGR2GRAY);
	blur(grayImage,blurImage,Size(5,5));
	Canny(blurImage,canny_detection,1,3,3);//边缘检测

	imwrite("./resources/" + baseFileName +"_canny.bmp",canny_detection);
	imshow("canny算子边缘检测",canny_detection);
	dstCannyKernel = imread("./resources/" + baseFileName +"_canny.bmp");
	imshow("a", dstCannyKernel);
	int i, j;
	int cPixelR, cPixelG, cPixelB;//currentPoint;
	for (i = 0; i < dstCannyKernel.rows; i++)
		for (j = 0; j < dstCannyKernel.cols; j++)
		{
			cPixelB = dstCannyKernel.at<Vec3b>(i, j)[0];
			cPixelG = dstCannyKernel.at<Vec3b>(i, j)[1];
			cPixelR = dstCannyKernel.at<Vec3b>(i, j)[2];
			if (cPixelB ==255 && cPixelR ==255 && cPixelG == 255)
			{
				dstCannyKernel.at<Vec3b>(i, j)[0] = 0;
				dstCannyKernel.at<Vec3b>(i, j)[1] = 128;
				dstCannyKernel.at<Vec3b>(i, j)[2] = 0;
			}else
			{
				dstCannyKernel.at<Vec3b>(i, j)[0] = 0;
				dstCannyKernel.at<Vec3b>(i, j)[1] = 0;
				dstCannyKernel.at<Vec3b>(i, j)[2] = 0;
			}


		}
	imwrite("./resources/" + baseFileName +"_canny_green.bmp", dstCannyKernel);
	imshow("da", dstCannyKernel);
}

void Image_Conv::doSharpenKernel()
{
	Mat sharpen_kernel = (Mat_<double>(3, 3) <<
		-1, -1, -1,
		-1, 9, -1,
		-1, -1, -1);
	filter2D(srcMapImage, dstSharpenImage, srcMapImage.depth(), sharpen_kernel);
	imwrite("./resources/" + baseFileName +"_sharpen.jpg", dstSharpenImage);
	imshow("锐化卷积图", dstSharpenImage);
}


void Image_Conv::makeObjectDepartBgBlack(int i, int j, int cPointB, int cPointG, int cPointR)
{
	//背景为黑 前景分块
	if ((cPointB == 167 && cPointG ==232 && cPointR ==207))//grass->green
	{
		srcMapImage.at<Vec3b>(i, j)[0] = 0;
		srcMapImage.at<Vec3b>(i, j)[1] = 255;
		srcMapImage.at<Vec3b>(i, j)[2] = 0;
	}
	else if ((cPointB == 237 && cPointG ==241 && cPointR ==239)||(cPointB == 235 && cPointG ==241 && cPointR ==243))
	{
		srcMapImage.at<Vec3b>(i, j)[0] = 0;
		srcMapImage.at<Vec3b>(i, j)[1] = 0;
		srcMapImage.at<Vec3b>(i, j)[2] = 255;
	}
	else if((cPointB == 255 && cPointG ==210 && cPointR ==173))
	{
		srcMapImage.at<Vec3b>(i, j)[0] = 255;
		srcMapImage.at<Vec3b>(i, j)[1] = 0;
		srcMapImage.at<Vec3b>(i, j)[2] = 0;
	}
	else
	{
		srcMapImage.at<Vec3b>(i, j)[0] = 0;
		srcMapImage.at<Vec3b>(i, j)[1] = 0;
		srcMapImage.at<Vec3b>(i, j)[2] = 0;
	}

}

void Image_Conv::makeBgBlack(int i, int j, int cPointB, int cPointG, int cPointR)
{
	//将非建筑物、grass、lake的背景涂黑
	if (!((cPointB == 167 && cPointG == 232 && cPointR == 207) || (cPointB == 237 && cPointG == 241 && cPointR == 239) || (cPointB == 255 && cPointG == 210 && cPointR == 173) || (cPointB == 235 && cPointG == 241 && cPointR == 243)))
	{
		srcMapImage.at<Vec3b>(i, j)[0] = 0;
		srcMapImage.at<Vec3b>(i, j)[1] = 0;
		srcMapImage.at<Vec3b>(i, j)[2] = 0;
	}
}

void Image_Conv::makeBgBlackAndForeWhite(int i, int j, int cPointB, int cPointG, int cPointR)
{
	//背景为黑 前景为白
	if (!((cPointB == 167 && cPointG ==232 && cPointR ==207)||(cPointB == 237 && cPointG ==241 && cPointR ==239)||(cPointB == 255 && cPointG ==210 && cPointR ==173)||(cPointB == 235 && cPointG ==241 && cPointR ==243)))
	{
		srcMapImage.at<Vec3b>(i, j)[0] = 0;
		srcMapImage.at<Vec3b>(i, j)[1] = 0;
		srcMapImage.at<Vec3b>(i, j)[2] = 0;
	}
	else
	{
		srcMapImage.at<Vec3b>(i, j)[0] = 255;
		srcMapImage.at<Vec3b>(i, j)[1] = 255;
		srcMapImage.at<Vec3b>(i, j)[2] = 255;
	}
}

void Image_Conv::doSplitBuildAndGrasslake(cv::Mat srcGrassAndLakeImage, cv::Mat srcBuildingImage, int i, int j, int cPointB, int cPointG, int cPointR)
{
	//if ((cPointB == 237 && cPointG == 241 && cPointR == 239) || (cPointB == 235 && cPointG == 241 && cPointR == 243))
	//{
	//	srcBuildingImage.at<Vec3b>(i, j)[0] = 255;
	//	srcBuildingImage.at<Vec3b>(i, j)[1] = 255;
	//	srcBuildingImage.at<Vec3b>(i, j)[2] = 255;
	//}
	//else if ((cPointB == 167 && cPointG ==232 && cPointR ==207) || (cPointB == 255 && cPointG == 210 && cPointR == 173))
	//{
	//	srcGrassAndLakeImage.at<Vec3b>(i, j)[0] = 255;
	//	srcGrassAndLakeImage.at<Vec3b>(i, j)[1] = 255;
	//	srcGrassAndLakeImage.at<Vec3b>(i, j)[2] = 255;
	//}


	if ((cPointB == 237 && cPointG == 241 && cPointR == 239) || (cPointB == 235 && cPointG == 241 && cPointR == 243))//building
	{
		srcBuildingImage.at<Vec3b>(i, j)[0] = 0;
		srcBuildingImage.at<Vec3b>(i, j)[1] = 0;
		srcBuildingImage.at<Vec3b>(i, j)[2] = 255;
	}
	else if((cPointB == 167 && cPointG ==232 && cPointR ==207))//grass
	{
		srcGrassAndLakeImage.at<Vec3b>(i, j)[0] = 0;
		srcGrassAndLakeImage.at<Vec3b>(i, j)[1] = 255;
		srcGrassAndLakeImage.at<Vec3b>(i, j)[2] = 0;
	}
	else if ((cPointB == 255 && cPointG == 210 && cPointR == 173))//lake
	{
		srcGrassAndLakeImage.at<Vec3b>(i, j)[0] = 255;
		srcGrassAndLakeImage.at<Vec3b>(i, j)[1] = 0;
		srcGrassAndLakeImage.at<Vec3b>(i, j)[2] = 0;
	}
}

void Image_Conv::doMergeTwoImages(cv::Mat srcGrassAndLakeImage, cv::Mat srcBuildingImage,float weigthOne, float weightTwo)
{
	for (int i = 0; i < srcMapImage.rows; i++)
		for (int j = 0; j < srcMapImage.cols; j++)
		{
			
			srcMapImage.at<Vec3b>(i, j)[0] = srcBuildingImage.at<Vec3b>(i, j)[0] * weigthOne + srcGrassAndLakeImage.at<Vec3b>(i, j)[0] * weightTwo;
			srcMapImage.at<Vec3b>(i, j)[1] = srcBuildingImage.at<Vec3b>(i, j)[1] * weigthOne+ srcGrassAndLakeImage.at<Vec3b>(i, j)[1] * weightTwo;
			srcMapImage.at<Vec3b>(i, j)[2] = srcBuildingImage.at<Vec3b>(i, j)[2] * weigthOne+ srcGrassAndLakeImage.at<Vec3b>(i, j)[2] * weightTwo;
		}
	imshow("overlap srcimage", srcMapImage);

}

static void doExportContoursPosToPolyFile(vector<vector<Point>> contours)
{
	ofstream contoursPosition;
	contoursPosition.open(".\\position130.poly", ios_base::app);
	for (int i = 0; i < contours.size(); i++)
	{
		contoursPosition << contours[i].size() << " " << 2 << " " << 0 << " " << 1 << endl;
		for (int j_i=0; j_i < contours[i].size();j_i++)
		{
			contoursPosition << j_i <<" "<< contours[i][j_i].x <<" "<< contours[i][j_i].y << " "<<1 <<endl;
		}
		contoursPosition << contours[i].size() << " " <<  1 << endl;
		for (int j_i = 0; j_i < contours[i].size(); j_i++)
		{
			if (j_i == (contours[i].size() - 1))
			{
				contoursPosition << j_i << " " << j_i << " " << 0 << " "<< 1<<endl;

			}
			else
			{
				contoursPosition << j_i << " " << j_i << " " << j_i + 1 <<  " "<< 1<<endl;
			}
		}
		contoursPosition <<  0 << endl;

	}
	contoursPosition.close();
}

static void doExportContoursPosToPolyFileSimple(vector<vector<Point>> contours)
{
	ofstream contoursPosition;
	contoursPosition.open(".\\position.poly", ios_base::app);
	for (int i = 0; i < contours.size(); i++)
	{
		contoursPosition << contours[i].size() <<endl;
		for (int j_i = 0; j_i < contours[i].size(); j_i++)
		{
			contoursPosition << contours[i][j_i].x << " " << contours[i][j_i].y << endl;
		}
	}
	contoursPosition.close();
}

/*采用cvFindContours提取轮廓，并过滤掉小面积轮廓，最后将轮廓保存*/
int Image_Conv::getContoursByCplus(cv::Mat& srcImage, double minarea, double whRatio)
{
	cv::Mat src, dst, canny_output,grayImage,img_binary;
	/// Load source image and convert it to gray
	//src = srcImage;

	//if (!src.data)
	//{
	//	std::cout << "read data error!" << std::endl;
	//	return -1;
	//}
	cvtColor(srcImage,grayImage,COLOR_BGR2GRAY);
	
	threshold(grayImage, img_binary, 1, 255, THRESH_BINARY);

	//imshow("grayImage", grayImage);
	
	//blur(grayImage, src, Size(3, 3));
	//imshow("模糊",src);
	/// Detect edges using canny
	//Canny(grayImage, canny_output, 1,3,3);
	//imshow("canny" ,img_binary);
	//cvWaitKey(0);
	//Mat lap_kernel = (Mat_<double>(3, 3) <<
	//	1, 1, 1,
	//	1, -8, 1,
	//	1, 1, 1);
	//filter2D(src, canny_output, src.depth(), lap_kernel);
	//the pram. for findContours,
	vector<vector<Point> > contours;
	vector<Vec4i> hierarchy;


	/// Find contours
	findContours(img_binary, contours, hierarchy, CV_RETR_LIST, CV_CHAIN_APPROX_NONE, Point(0, 0));
	//CV_RETR_CCOMP, CV_CHAIN_APPROX_SIMPLE  CV_RETR_LIST CV_RETR_EXTERNAL

	//std::cout << "size of image canny_output" << canny_output.rows<<" " <<canny_output.cols<< std::endl;

	double maxarea = 0;
	int maxAreaIdx = 0;



	for (int i = 0; i < contours.size(); i++)
	{

		double tmparea = fabs(contourArea(contours[i]));
		if (tmparea > maxarea)
		{
			maxarea = tmparea;
			maxAreaIdx = i;
			continue;
		}

		if (tmparea < minarea)
		{
			//删除面积小于设定值的轮廓
			contours.erase(contours.begin() + i);
			std::wcout << "delete a small area" << std::endl;
			continue;
		}
		//计算轮廓的直径宽高
		Rect aRect = boundingRect(contours[i]);
		if ((aRect.width / aRect.height) < whRatio)
		{
			//删除宽高比例小于设定值的轮廓
			contours.erase(contours.begin() + i);
			std::wcout << "delete a unnomalRatio area" << std::endl;
			continue;
		}
	}

	//Because the lines from dt_distance_field result is thick, and use findContours function will get two contours for each line,so here we use m_i%2 to save the outside of each line
	vector<vector<Point> > tmpContours;
	for (int m_i = contours.size()-1; m_i >= 0; --m_i)
	{
		if ((m_i % 2) != 0)
		{
			tmpContours.push_back(contours[m_i]);
		}
	}

	contours = tmpContours;


	// Draw contours,绿色轮廓
	dst = Mat::zeros(srcImage.size(), CV_8UC3);

	//std::cout << "size of image dst" << dst.rows << " " << dst.cols << std::endl;
	for (int i = 0; i < contours.size(); i++)
	{
		//contours[i].resize(contours[i].size() - 5);
		//随机颜色
		//Scalar color = Scalar( rng.uniform(0, 255), rng.uniform(0, 255), rng.uniform(0, 255));
		Scalar color = Scalar( 0,128,0);
		drawContours(dst, contours, i, color, 2, 8, hierarchy, 0, Point());
		
	}


	//write contours position to poly file
	//poly file format：https://www.cs.cmu.edu/~quake/triangle.poly.html
	
	vector<vector<Point> > contoursTemp(contours.size());//((contours[0].size()) % 10);

	if (contours.size()>0)
	{
		for (int contours_x = 0; contours_x < contours.size(); ++contours_x)
		{
			for (int controus_y = 0; controus_y < contours[contours_x].size(); ++controus_y)
			{
				if ((controus_y % 37) == 0)//按照高度37pxiel的粗略间隔
				{
					cv::Point ttx = contours[contours_x][controus_y];
					contoursTemp[contours_x].push_back(ttx);
				}
			}
		}
	}
	else
	{
		std::cout << "There is no controus" << endl;
	}

	//doExportContoursPosToPolyFile(contours);
	//doExportContoursPosToPolyFileSimple(contours);

	//Fermat::init_fermat(contoursTemp);//contours
	Fermat::get_fermat_from_contour(contoursTemp);//contours
	cv::imwrite("./resources/" + baseFileName +"_green.bmp",dst);
	cv::imshow("countors", dst);
	cv::waitKey(0);

	return 0;
}
void Image_Conv::doGridMap(cv::Mat& srcImage)
{
	int pixelOfGrid = 40;
	std::cout << "开始grid"<< std::endl;
	int gridColNum =  srcImage.cols / pixelOfGrid ;
	int gridRowNum = srcImage.rows / pixelOfGrid;

	cv::Mat gridImage(srcImage.rows, srcImage.cols, CV_8UC3, Scalar(0, 0, 0));
	vector<vector<int> > gridColor(gridRowNum,vector<int>(gridColNum));

	for (int i = 0; i < gridRowNum; i++)
		for (int j = 0; j < gridColNum; j++)
		{
			float cNumOfBuildingPixel = 0;
			float cNumOfGrassPixel = 0;
			float cNumOfLakePixel = 0;
			for (int m = i*pixelOfGrid; m < (i+1)*pixelOfGrid-1; m++)
				for (int n = j*pixelOfGrid; n < (j+1)*pixelOfGrid-1; n++)
				{
					int cPointB = srcImage.at<Vec3b>(m, n)[0];
					int cPointG = srcImage.at<Vec3b>(m, n)[1];
					int cPointR = srcImage.at<Vec3b>(m, n)[2];
					if ((cPointB == 237 && cPointG ==241 && cPointR ==239)||(cPointB == 235 && cPointG ==241 && cPointR ==243))
					{
						cNumOfBuildingPixel++;
					} 
					else if ((cPointB == 167 && cPointG ==232 && cPointR ==207))
					{
						cNumOfGrassPixel++;
					}
					else if ((cPointB == 255 && cPointG ==210 && cPointR ==173))
					{
						cNumOfLakePixel++;
					}
				}
			int cColorBOfGrid = cNumOfLakePixel / pixelOfGrid*pixelOfGrid * 255;
			int cColorROfGrid = cNumOfBuildingPixel / pixelOfGrid*pixelOfGrid * 255;
			int cColorGOfGrid = cNumOfGrassPixel / pixelOfGrid*pixelOfGrid * 255;
			
			for (int gm = i * pixelOfGrid; gm < (i+1) * pixelOfGrid - 1; gm++)
				for (int gn = j * pixelOfGrid; gn < (j+1) * pixelOfGrid - 1; gn++)
				{
					gridImage.at<Vec3b>(gm, gn)[0] = cColorBOfGrid;
					gridImage.at<Vec3b>(gm, gn)[1] = cColorGOfGrid;
					gridImage.at<Vec3b>(gm, gn)[2] = cColorROfGrid;
				}

		}
	std::cout << "结束grid" << std::endl;
	imshow("gridImage", gridImage);
	imwrite("./resources/gridImage.bmp", gridImage);
	//biInterpolation(gridImage);
	cvWaitKey(0);
	//biInterpolation(gridImage);
}

void Image_Conv::biInterpolation(cv::Mat& src)
{
	//调整对比度和亮度
	Mat matSrc, matDst1, matDst2;
	matSrc = src;
	matDst1 = Mat(Size(src.cols, src.rows), matSrc.type());
	matDst2 = Mat(matDst1.size(), matSrc.type());
	double scale_x = (double)matSrc.cols / matDst2.cols;
	double scale_y = (double)matSrc.rows / matDst2.rows;

	uchar* dataDst = matDst2.data;
	int stepDst = matDst2.step;//1024*3=3072
	uchar* dataSrc = src.data;
	int stepSrc = src.step;//1536
	int iWidthSrc = src.cols;//512
	int iHiehgtSrc = src.rows;//512
	for (int j = 0; j < matDst2.rows; ++j)
	{
		float fy = (float)((j + 0.5) * scale_y - 0.5);//-0.25
		int sy = cvFloor(fy);//-1
		fy -= sy;//0.75
		sy = std::min(sy, iHiehgtSrc - 2);//-1
		sy = std::max(0, sy);//0

		short cbufy[2];
		cbufy[0] = cv::saturate_cast<short>((1.f - fy) * 2048);//512
		cbufy[1] = 2048 - cbufy[0];//1536

		for (int i = 0; i < matDst2.cols; ++i)
		{
			float fx = (float)((i + 0.5) * scale_x - 0.5);//-0.25
			int sx = cvFloor(fx);//-1
			fx -= sx;//0.75

			if (sx < 0) {
				fx = 0, sx = 0;
			}
			if (sx >= iWidthSrc - 1) {
				fx = 0, sx = iWidthSrc - 2;
			}
			//这里不大懂
			short cbufx[2];
			cbufx[0] = cv::saturate_cast<short>((1.f - fx) * 2048);//2048
			cbufx[1] = 2048 - cbufx[0];//0

			for (int k = 0; k < src.channels(); ++k)
			{
				*(dataDst + j*stepDst + 3 * i + k) =
					(*(dataSrc + sy*stepSrc + 3 * sx + k) * cbufx[0] * cbufy[0] +
					*(dataSrc + (sy + 1)*stepSrc + 3 * sx + k) * cbufx[0] * cbufy[1] +
					*(dataSrc + sy*stepSrc + 3 * (sx + 1) + k) * cbufx[1] * cbufy[0] +
					*(dataSrc + (sy + 1)*stepSrc + 3 * (sx + 1) + k) * cbufx[1] * cbufy[1]) >> 22;
			}
		}
	}
	imshow("matDst2",matDst2);


}

void Image_Conv::init()
{
	
	string fileSrcImagePath = "./resources/white_contour_spiral.bmp";
	baseFileName = FileSystem::base_name(fileSrcImagePath);
	srcMapImage = imread(fileSrcImagePath);

	//imshow("原图", srcMapImage);
	cv::Mat srcGrassAndLakeImage(srcMapImage.rows, srcMapImage.cols, CV_8UC3, Scalar(0, 0, 0));
	cv::Mat srcBuildingImage(srcMapImage.rows, srcMapImage.cols, CV_8UC3, Scalar(0, 0, 0));

	//doGridMap(srcMapImage);
	//addOutlierOnImage(srcMapImage);
	//cvWaitKey(0);

	//filterImage(srcMapImage);
	//KMeansSegImage(srcMapImage);
	int i, j;
	for (i = 0; i < srcMapImage.rows; i++)
		for (j = 0; j < srcMapImage.cols; j++)
		{
			int cPixelB = srcMapImage.at<Vec3b>(i, j)[0];
			int cPixelG = srcMapImage.at<Vec3b>(i, j)[1];
			int cPixelR = srcMapImage.at<Vec3b>(i, j)[2];
			
			//gaode
			//bgr:
			//building:235,241,243	237,241,239
			//grass:167,232,207
			//lake:255,210,173

			//makeObjectDepartBgBlack(i,j,cPixelB, cPixelG, cPixelR);
			//makeBgBlack(i, j, cPixelB, cPixelG, cPixelR);
			//makeBgBlackAndForeWhite(i, j, cPixelB, cPixelG, cPixelR);
			//建筑物与grassAndLake分离处理
			//doSplitBuildAndGrasslake(srcGrassAndLakeImage, srcBuildingImage, i, j, cPixelB, cPixelG, cPixelR);

			//reserveInterestRegion(i,j,cPixelB, cPixelG, cPixelR);

		}
	//doLaplaceKernel();
	//imshow("region interest", srcMapImage);
	//getContoursByCplus(dstLapImage, 0, 0);

	//imshow("建筑物原图", srcBuildingImage);
	//imshow("绿地和湖水原图", srcGrassAndLakeImage);
	//cvWaitKey(0);


	
	//doExpandAndErosionImage(srcMapImage);
	//srcMapImage = outErosionImage;
	//doLaplaceKernel();
	getContoursByCplus(srcMapImage, 2500, 0);

	//doMergeTwoImages(srcGrassAndLakeImage,outErosionImage);
	
	//imwrite("./resources/" + baseFileName +"_BgBlackAndrgb.bmp",srcMapImage);

	
	//for (i = 0; i < srcMapImage.rows; i++)
	//	for (j = 0; j < srcMapImage.cols; j++)
	//	{
	//		int cPixelB = srcMapImage.at<Vec3b>(i, j)[0];
	//		int cPixelG = srcMapImage.at<Vec3b>(i, j)[1];
	//		int cPixelR = srcMapImage.at<Vec3b>(i, j)[2];

	//		if ((cPixelB == 255 && cPixelG ==255 && cPixelR ==255))
	//		{
	//			srcMapImage.at<Vec3b>(i, j)[0] = 0;
	//			srcMapImage.at<Vec3b>(i, j)[1] = 0;
	//			srcMapImage.at<Vec3b>(i, j)[2] = 0;
	//		}
	//		else
	//		{
	//			srcMapImage.at<Vec3b>(i, j)[0] = 255;
	//			srcMapImage.at<Vec3b>(i, j)[1] = 255;
	//			srcMapImage.at<Vec3b>(i, j)[2] = 255;
	//		}

	//	}

	//imshow("建筑物反转", srcMapImage);
	//imwrite("./resources/" + baseFileName +"_BgBlackAndForeWhite.bmp",srcMapImage);

	//Mat src2 = imread("./resources/gaode_szu-big_abstract_distance-grass.bmp");
	//Mat src1 = imread("./resources/gaode_szu-big_abstract_distance-build.bmp");
	//doMergeTwoImages(src1, src2,0.7,0.3);

	//getContoursByCplus(srcGrassAndLakeImage, 0, 0);
	//getContoursByCplus(srcBuildingImage,0,0);

	//doLaplaceKernel();
	//doSharpenKernel();
	//doCannyKernel();

	
	//waitKey(0);
	return ;
}

