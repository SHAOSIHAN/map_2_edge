## 地图边检测

输入：2D地图

输出：建筑物、绿地和水的外部轮廓（膨胀与腐蚀之后的结果，填补空洞及删除噪音）

算法流程：

1. 对原始图片进行分离，一类是建筑物等关键信息，另一类是绿地与水等非关键信息；
2. 对建筑物图进行膨胀收缩；
3. 将上述两张图融合；
4. 边检测等算法，有效提取外部轮廓

![image](http://wx2.sinaimg.cn/large/6e529308gy1fvp8mgykprj214e0h3424.jpg)

## 代码结构

### image_cov

图片处理相关代码

init:代码入口

```c++
addOutlierOnImage(cv::Mat& srcImage);	//对image边缘添加一个像素
filterImage(cv::Mat& srcImage);	//滤波
KMeansSegImage(cv::Mat& srcImage);	//K-means分割
doExpandAndErosionImage(cv::Mat& srcBuildingImage);	//对image进行腐蚀膨胀操作
doLaplaceKernel();	//laplace卷积
doCannyKernel();	//canny
doSharpenKernel();	//图片锐化
makeBgBlack(int i, int j, int cPointB, int cPointG, int cPointR);	//使图片背景变为黑色
makeBgBlackAndForeWhite(int i, int j, int cPointB, int cPointG, int cPointR);	//目标物为白色背景为黑色
makeObjectDepartBgBlack(int i, int j, int c_pixel_b, int c_pixel_g, int c_pixel_r);	//背景为黑，目标物赋予不同的颜色
doSplitBuildAndGrasslake(cv::Mat srcGrassAndLakeImage, cv::Mat srcBuildingImage, int i, int j, int cPointB,int cPointG, int cPointR);	//建筑物和草地湖水分离
doMergeTwoImages(cv::Mat srcGrassAndLakeImage, cv::Mat srcBuildingImage, float weigthOne = 1.0, float weightTwo = 1.0);	//融合两张照片
getContoursByCplus(cv::Mat& srcImage, double minarea, double whRatio);	//利用findcontours方法进行边界提取
reserveInterestRegion(int i, int j, int cPointB, int cPointG, int cPointR);	
init();	//图片处理入口
//grid map
doGridMap(cv::Mat& srcImage);	//将图片网格化，建筑物草地湖水有不同的权重
biInterpolation(cv::Mat& src);	//双线性插值
```

### fermat

生成fermat线的代码

利用`getContoursByCplus`获取到的contours进行fermat化，具体实现方法见**Connected fermat spirals for layered fabrication**这篇文章