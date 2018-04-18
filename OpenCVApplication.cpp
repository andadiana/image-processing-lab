// OpenCVApplication.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "common.h"
#include <random>
#include <fstream>

using namespace std;

void testOpenImage()
{
	char fname[MAX_PATH];
	while (openFileDlg(fname))
	{
		Mat src;
		src = imread(fname);
		imshow("image", src);
		waitKey();
	}
}

void testOpenImagesFld()
{
	char folderName[MAX_PATH];
	if (openFolderDlg(folderName) == 0)
		return;
	char fname[MAX_PATH];
	FileGetter fg(folderName, "bmp");
	while (fg.getNextAbsFile(fname))
	{
		Mat src;
		src = imread(fname);
		imshow(fg.getFoundFileName(), src);
		if (waitKey() == 27) //ESC pressed
			break;
	}
}


void testResize()
{
	char fname[MAX_PATH];
	while (openFileDlg(fname))
	{
		Mat src;
		src = imread(fname);
		Mat dst1, dst2;
		//without interpolation
		resizeImg(src, dst1, 320, false);
		//with interpolation
		resizeImg(src, dst2, 320, true);
		imshow("input image", src);
		imshow("resized image (without interpolation)", dst1);
		imshow("resized image (with interpolation)", dst2);
		waitKey();
	}
}


void testVideoSequence()
{
	VideoCapture cap("Videos/rubic.avi"); // off-line video from file
										  //VideoCapture cap(0);	// live video from web cam
	if (!cap.isOpened()) {
		printf("Cannot open video capture device.\n");
		waitKey(0);
		return;
	}

	Mat edges;
	Mat frame;
	char c;

	while (cap.read(frame))
	{
		Mat grayFrame;
		cvtColor(frame, grayFrame, CV_BGR2GRAY);
		imshow("source", frame);
		imshow("gray", grayFrame);
		c = cvWaitKey(0);  // waits a key press to advance to the next frame
		if (c == 27) {
			// press ESC to exit
			printf("ESC pressed - capture finished\n");
			break;  //ESC pressed
		};
	}
}


void testSnap()
{
	VideoCapture cap(0); // open the deafult camera (i.e. the built in web cam)
	if (!cap.isOpened()) // openenig the video device failed
	{
		printf("Cannot open video capture device.\n");
		return;
	}

	Mat frame;
	char numberStr[256];
	char fileName[256];

	// video resolution
	Size capS = Size((int)cap.get(CV_CAP_PROP_FRAME_WIDTH),
		(int)cap.get(CV_CAP_PROP_FRAME_HEIGHT));

	// Display window
	const char* WIN_SRC = "Src"; //window for the source frame
	namedWindow(WIN_SRC, CV_WINDOW_AUTOSIZE);
	cvMoveWindow(WIN_SRC, 0, 0);

	const char* WIN_DST = "Snapped"; //window for showing the snapped frame
	namedWindow(WIN_DST, CV_WINDOW_AUTOSIZE);
	cvMoveWindow(WIN_DST, capS.width + 10, 0);

	char c;
	int frameNum = -1;
	int frameCount = 0;

	for (;;)
	{
		cap >> frame; // get a new frame from camera
		if (frame.empty())
		{
			printf("End of the video file\n");
			break;
		}

		++frameNum;

		imshow(WIN_SRC, frame);

		c = cvWaitKey(10);  // waits a key press to advance to the next frame
		if (c == 27) {
			// press ESC to exit
			printf("ESC pressed - capture finished");
			break;  //ESC pressed
		}
		if (c == 115) { //'s' pressed - snapp the image to a file
			frameCount++;
			fileName[0] = NULL;
			sprintf(numberStr, "%d", frameCount);
			strcat(fileName, "Images/A");
			strcat(fileName, numberStr);
			strcat(fileName, ".bmp");
			bool bSuccess = imwrite(fileName, frame);
			if (!bSuccess)
			{
				printf("Error writing the snapped image\n");
			}
			else
				imshow(WIN_DST, frame);
		}
	}

}

void MyCallBackFunc(int event, int x, int y, int flags, void* param)
{
	//More examples: http://opencvexamples.blogspot.com/2014/01/detect-mouse-clicks-and-moves-on-image.html
	Mat* src = (Mat*)param;
	if (event == CV_EVENT_LBUTTONDOWN)
	{
		printf("Pos(x,y): %d,%d  Color(RGB): %d,%d,%d\n",
			x, y,
			(int)(*src).at<Vec3b>(y, x)[2],
			(int)(*src).at<Vec3b>(y, x)[1],
			(int)(*src).at<Vec3b>(y, x)[0]);
	}
}

void testMouseClick()
{
	Mat src;
	// Read image from file 
	char fname[MAX_PATH];
	while (openFileDlg(fname))
	{
		src = imread(fname);
		//Create a window
		namedWindow("My Window", 1);

		//set the callback function for any mouse event
		setMouseCallback("My Window", MyCallBackFunc, &src);

		//show the image
		imshow("My Window", src);

		// Wait until user press some key
		waitKey(0);
	}
}

/* Histogram display function - display a histogram using bars (simlilar to L3 / PI)
Input:
name - destination (output) window name
hist - pointer to the vector containing the histogram values
hist_cols - no. of bins (elements) in the histogram = histogram image width
hist_height - height of the histogram image
Call example:
showHistogram ("MyHist", hist_dir, 255, 200);
*/
void showHistogram(const std::string& name, int* hist, const int  hist_cols, const int hist_height)
{
	Mat imgHist(hist_height, hist_cols, CV_8UC3, CV_RGB(255, 255, 255)); // constructs a white image

																		 //computes histogram maximum
	int max_hist = 0;
	for (int i = 0; i<hist_cols; i++)
		if (hist[i] > max_hist)
			max_hist = hist[i];
	double scale = 1.0;
	scale = (double)hist_height / max_hist;
	int baseline = hist_height - 1;

	for (int x = 0; x < hist_cols; x++) {
		Point p1 = Point(x, baseline);
		Point p2 = Point(x, baseline - cvRound(hist[x] * scale));
		line(imgHist, p1, p2, CV_RGB(255, 0, 255)); // histogram bins colored in magenta
	}

	imshow(name, imgHist);
}

// Lab 1
void negative_image() {
	Mat img = imread("Images/cameraman.bmp",
		CV_LOAD_IMAGE_GRAYSCALE);
	for (int i = 0; i<img.rows; i++) {
		for (int j = 0; j<img.cols; j++) {
			img.at<uchar>(i, j) = 255 - img.at<uchar>(i, j);
		}
	}
	imshow("negative image", img);
	waitKey(0);
}

uchar const ADDITIVE_VAL = 20;

void gray_levels_additive() {
	Mat img = imread("Images/saturn.bmp",
		CV_LOAD_IMAGE_GRAYSCALE);
	imshow("Before", img);
	for (int i = 0; i<img.rows; i++) {
		for (int j = 0; j<img.cols; j++) {
			int new_val = img.at<uchar>(i, j) + ADDITIVE_VAL;
			if (new_val > 255) {
				img.at<uchar>(i, j) = 255;
			}
			else if (new_val < 0) {
				img.at<uchar>(i, j) = 0;
			}
			else {
				img.at<uchar>(i, j) = new_val;
			}
		}
	}
	imshow("After", img);
	waitKey(0);
}

uchar const MULTIPLICATIVE_VAL = 2;

void gray_levels_multiplicative() {
	Mat img = imread("Images/cameraman.bmp",
		CV_LOAD_IMAGE_GRAYSCALE);
	imshow("Before", img);
	for (int i = 0; i<img.rows; i++) {
		for (int j = 0; j<img.cols; j++) {
			int new_val = img.at<uchar>(i, j) * MULTIPLICATIVE_VAL;
			if (new_val > 255)
				img.at<uchar>(i, j) = 255;
			else if (new_val < 0)
				img.at<uchar>(i, j) = 0;
			else
				img.at<uchar>(i, j) = new_val;
		}
	}
	imshow("After", img);
	imwrite("cameraman_mul.bmp", img);
	waitKey(0);
}

void create_color_image() {
	Mat img(256, 256, CV_8UC3);
	for (int i = 0; i<img.rows; i++) {
		for (int j = 0; j<img.cols; j++) {
			Vec3b pixel;
			if (i < 128 && j < 128) {
				pixel[0] = 255;
				pixel[1] = 255;
				pixel[2] = 255;
			}
			else if (i < 128 && j >= 128) {
				pixel[0] = 0;
				pixel[1] = 0;
				pixel[2] = 255;
			}
			else if (i >= 128 && j < 128) {
				pixel[0] = 0;
				pixel[1] = 255;
				pixel[2] = 0;
			}
			else {
				pixel[0] = 0;
				pixel[1] = 255;
				pixel[2] = 255;
			}
			img.at<Vec3b>(i, j) = pixel;
		}
	}
	imshow("Color image", img);
	waitKey(0);
}

void horizontal_flip() {

	Mat img = imread("Images/cameraman.bmp",
		CV_LOAD_IMAGE_GRAYSCALE);
	imshow("Before flipping", img);
	for (int i = 0; i<img.rows; i++) {
		for (int j = 0; j<img.cols / 2; j++) {
			uchar aux = img.at<uchar>(i, j);
			img.at<uchar>(i, j) = img.at<uchar>(i, img.cols - j - 1);
			img.at<uchar>(i, img.cols - j - 1) = aux;
		}
	}
	imshow("After flipping", img);
	waitKey(0);
}

void vertical_flip() {

	Mat img = imread("Images/cameraman.bmp",
		CV_LOAD_IMAGE_GRAYSCALE);
	imshow("Before flipping", img);
	for (int i = 0; i<img.rows / 2; i++) {
		for (int j = 0; j<img.cols; j++) {
			uchar aux = img.at<uchar>(i, j);
			img.at<uchar>(i, j) = img.at<uchar>(img.rows - i - 1, j);
			img.at<uchar>(img.rows - i - 1, j) = aux;
		}
	}
	imshow("After flipping", img);
	waitKey(0);
}

void center_crop(int height, int width) {
	Mat img = imread("Images/kids.bmp",
		CV_LOAD_IMAGE_COLOR);
	imshow("Original", img);
	Mat cropped(img.rows, img.cols, CV_8UC3);
	if (height > img.rows || width > img.cols) {
		printf("Incorrect crop size! Please input another");
		return;
	}
	int init_r_pos = img.rows / 2 - height / 2;
	int init_c_pos = img.cols / 2 - width / 2;
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			cropped.at<Vec3b>(i, j) = img.at<Vec3b>(init_r_pos + i, init_c_pos + j);
		}
	}
	imshow("Cropped", cropped);
	waitKey(0);
}

void resize_image(int new_height, int new_width) {
	Mat img = imread("Images/cameraman.bmp",
		CV_LOAD_IMAGE_COLOR);
	imshow("Original", img);
	printf("original size is: rows: %d, cols: %d\n", img.rows, img.cols);
	Mat resized(new_height, new_width, CV_8UC3);
	double height_ratio = double(img.rows) / new_height;
	double width_ratio = double(img.cols) / new_width;
	printf("height ratio is %.2f\n", height_ratio);
	printf("width ratio is %.2f\n", width_ratio);
	for (int i = 0; i < new_height; i++) {
		for (int j = 0; j < new_width; j++) {
			int nearest_i = round(i * height_ratio);
			int nearest_j = round(j * width_ratio);

			//printf("Nearest_i is %d\n", nearest_i);
			//printf("Nearest_j is %d\n", nearest_j);
			if (nearest_i >= img.rows) { nearest_i = img.rows - 1; }
			if (nearest_j >= img.cols) { nearest_j = img.cols - 1; }
			resized.at<Vec3b>(i, j) = img.at<Vec3b>(nearest_i, nearest_j);
			
		}
	}
	imshow("Resized", resized);
	waitKey(0);
}

// Lab 2
void separate_RGB() {
	Mat img = imread("Images/flowers_24bits.bmp",
		CV_LOAD_IMAGE_COLOR);
	imshow("Original", img);
	Mat imgR(img.rows, img.cols, CV_8UC1);
	Mat imgG(img.rows, img.cols, CV_8UC1);
	Mat imgB(img.rows, img.cols, CV_8UC1);

	for (int i = 0; i < img.rows; i++) {
		for (int j = 0; j < img.cols; j++) {
			imgB.at<uchar>(i, j) = img.at<Vec3b>(i, j)[0];
			imgG.at<uchar>(i, j) = img.at<Vec3b>(i, j)[1];
			imgR.at<uchar>(i, j) = img.at<Vec3b>(i, j)[2];
		}
	}

	imshow("R channel", imgR);
	imshow("G channel", imgG);
	imshow("B channel", imgB);
	waitKey(0);
}

void RGB_to_grayscale() {
	Mat img = imread("Images/flowers_24bits.bmp",
		CV_LOAD_IMAGE_COLOR);
	imshow("Original", img);
	Mat img_grayscale(img.rows, img.cols, CV_8UC1);
	for (int i = 0; i < img.rows; i++) {
		for (int j = 0; j < img.cols; j++) {
			uchar blue = img.at<Vec3b>(i, j)[0];
			uchar green = img.at<Vec3b>(i, j)[1];
			uchar red = img.at<Vec3b>(i, j)[2];
			uchar gray_value = (red + green + blue) / 3;
			img_grayscale.at<uchar>(i, j) = gray_value;
		}
	}
	imshow("Grayscale", img_grayscale);
	waitKey(0);
}

void grayscale_to_binary(unsigned int threshold) {
	Mat img = imread("Images/eight.bmp",
		CV_LOAD_IMAGE_GRAYSCALE);
	imshow("Original", img);
	Mat img_binary(img.rows, img.cols, CV_8UC1);
	for (int i = 0; i < img.rows; i++) {
		for (int j = 0; j < img.cols; j++) {
			uchar value = img.at<uchar>(i, j);
			if (value < threshold) {
				img_binary.at<uchar>(i, j) = 0;
			}
			else {
				img_binary.at<uchar>(i, j) = 255;
			}
		}
	}
	imshow("Grayscale", img_binary);
	waitKey(0);
}

void RGB_to_HSV() {
	Mat img = imread("Images/Lena_24bits.bmp",
		CV_LOAD_IMAGE_COLOR);
	imshow("Original", img);

	Mat imgH(img.rows, img.cols, CV_8UC1);
	Mat imgS(img.rows, img.cols, CV_8UC1);
	Mat imgV(img.rows, img.cols, CV_8UC1);

	for (int i = 0; i < img.rows; i++) {
		for (int j = 0; j < img.cols; j++) {
			float B = img.at<Vec3b>(i, j)[0];
			float G = img.at<Vec3b>(i, j)[1];
			float R = img.at<Vec3b>(i, j)[2];

			// normalize RGB values
			float r = R / 255;
			float g = G / 255;
			float b = B / 255;

			float max = max(max(r, g), b);
			float min = min(min(r, g), b);
			float C = max - min;

			//Value
			float V = max;

			//Saturation
			float S;
			if (V != 0) {
				S = C / V;
			} 
			else {
				S = 0; //grayscale
			}

			//Hue
			float H;
			if (C != 0.0) {
				if (max == r) H = 60 * (g - b) / C;
				if (max == g) H = 120 + 60 * (b - r) / C;
				if (max == b) H = 240 + 60 * (r - g) / C;
			}
			else {
				//grayscale
				H = 0;
			}
			if (H < 0.0) {
				H = H + 360;
			}

			//normalize HSV values
			H = (H * 255) / 360;
			S = S * 255;
			V = V * 255;

			imgH.at<uchar>(i, j) = H;
			imgS.at<uchar>(i, j) = S;
			imgV.at<uchar>(i, j) = V;
		}
	}
	imshow("H", imgH);
	imshow("S", imgS);
	imshow("V", imgV);
	waitKey(0);
}

void detect_traffic_sign() {
	Mat img = imread("Images/traffic_sign.png",
		CV_LOAD_IMAGE_COLOR);
	imshow("Original", img);
	Mat hsvImg;
	cv::cvtColor(img, hsvImg, CV_BGR2HSV);
	imshow("HSV image", hsvImg);

	Mat hsvChannels[3];
	cv::split(hsvImg, hsvChannels);

	Mat mask;
	cv::inRange(hsvChannels[0], 0, 15, mask);

	imshow("Mask", mask);
	waitKey(0);
}

void RGB_mask() {
	Mat img = imread("Images/traffic_sign.png",
		CV_LOAD_IMAGE_COLOR);
	imshow("Original", img);
	Mat hsvImg;
	cv::cvtColor(img, hsvImg, CV_BGR2HSV);
	imshow("HSV image", hsvImg);

	Mat hsvChannels[3];
	cv::split(hsvImg, hsvChannels);

	Mat redMask;
	cv::inRange(hsvChannels[0], 0, 7, redMask);

	Mat greenMask;
	cv::inRange(hsvChannels[0], 45, 65, greenMask);

	Mat blueMask;
	cv::inRange(hsvChannels[0], 90, 125, blueMask);

	imshow("Red mask", redMask);
	imshow("Green mask", greenMask);
	imshow("Blue mask", blueMask);
	waitKey(0);
}

int get_obj_area(Mat* img, uchar R, uchar G, uchar B) {
	int area = 0;
	for (int i = 0; i < (*img).rows; i++) {
		for (int j = 0; j < (*img).cols; j++) {
			uchar imgR = (*img).at<Vec3b>(i, j)[2];
			uchar imgG = (*img).at<Vec3b>(i, j)[1];
			uchar imgB = (*img).at<Vec3b>(i, j)[0];
			if (imgR == R && imgG == G && imgB == B) {
				area++;
			}

		}
	}
	return area;
}

int* get_obj_center_of_mass(Mat* img, uchar R, uchar G, uchar B) {
	int area = get_obj_area(img, R, G, B);
	int center_r = 0;
	int center_c = 0;
	for (int i = 0; i < (*img).rows; i++) {
		for (int j = 0; j < (*img).cols; j++) {
			uchar imgR = (*img).at<Vec3b>(i, j)[2];
			uchar imgG = (*img).at<Vec3b>(i, j)[1];
			uchar imgB = (*img).at<Vec3b>(i, j)[0];
			if (imgR == R && imgG == G && imgB == B) {
				center_r += i;
				center_c += j;
			}

		}
	}
	center_r = center_r / area;
	center_c = center_c / area;
	int center[2] = { center_r, center_c };
	return center;
}

double get_axis_of_elongation(Mat* img, uchar R, uchar G, uchar B) {
	int *center_of_mass;
	center_of_mass = get_obj_center_of_mass(img, R, G, B);
	double nominator = 0;
	double denominator = 0;
	int center_r = center_of_mass[0];
	int center_c = center_of_mass[1];
	for (int i = 0; i < (*img).rows; i++) {
		for (int j = 0; j < (*img).cols; j++) {
			uchar imgR = (*img).at<Vec3b>(i, j)[2];
			uchar imgG = (*img).at<Vec3b>(i, j)[1];
			uchar imgB = (*img).at<Vec3b>(i, j)[0];
			if (imgR == R && imgG == G && imgB == B) {
				nominator += ((i - center_r) * (j - center_c));
				denominator += (j -center_c) * (j - center_c) - (i - center_r) * (i - center_r);
			}
		}
	}
	nominator *= 2;
	double alpha = (atan2(nominator, denominator)) / 2;
	double angle_deg = (alpha * 180) / PI;
	return angle_deg;
}

boolean white_pixel(Vec3b pixel) {
	if (pixel[0] == 255 && pixel[1] == 255 && pixel[2] == 255)
		return true;
	return false;
}

int get_perimeter(Mat *img, uchar R, uchar G, uchar B) {
	int perimeter = 0;
	for (int i = 1; i < (*img).rows - 1; i++) {
		for (int j = 1; j < (*img).cols - 1; j++) {
			uchar imgR = (*img).at<Vec3b>(i, j)[2];
			uchar imgG = (*img).at<Vec3b>(i, j)[1];
			uchar imgB = (*img).at<Vec3b>(i, j)[0];
			if (imgR == R && imgG == G && imgB == B) {
				if (white_pixel((*img).at<Vec3b>(i - 1, j - 1)) ||
					white_pixel((*img).at<Vec3b>(i - 1, j)) ||
					white_pixel((*img).at<Vec3b>(i - 1, j + 1)) ||
					white_pixel((*img).at<Vec3b>(i, j - 1)) ||
					white_pixel((*img).at<Vec3b>(i, j + 1)) ||
					white_pixel((*img).at<Vec3b>(i + 1, j - 1)) ||
					white_pixel((*img).at<Vec3b>(i + 1, j)) ||
					white_pixel((*img).at<Vec3b>(i + 1, j + 1)))
				perimeter++;
			}
		}
	}
	return perimeter * (PI / 4);
}

void draw_contour(Mat *img, uchar R, uchar G, uchar B, Mat *newimg) {
	//Mat newimg(img->rows, img->cols, CV_8UC3);
	int perimeter = 0;
	
	for (int i = 1; i < (*img).rows - 1; i++) {
		for (int j = 1; j < (*img).cols - 1; j++) {
			uchar imgR = (*img).at<Vec3b>(i, j)[2];
			uchar imgG = (*img).at<Vec3b>(i, j)[1];
			uchar imgB = (*img).at<Vec3b>(i, j)[0];
			if (imgR == R && imgG == G && imgB == B) {
				if (white_pixel((*img).at<Vec3b>(i - 1, j - 1)) ||
					white_pixel((*img).at<Vec3b>(i - 1, j)) ||
					white_pixel((*img).at<Vec3b>(i - 1, j + 1)) ||
					white_pixel((*img).at<Vec3b>(i, j - 1)) ||
					white_pixel((*img).at<Vec3b>(i, j + 1)) ||
					white_pixel((*img).at<Vec3b>(i + 1, j - 1)) ||
					white_pixel((*img).at<Vec3b>(i + 1, j)) ||
					white_pixel((*img).at<Vec3b>(i + 1, j + 1))) {
					perimeter++;
					(*newimg).at<Vec3b>(i, j) = Vec3b(0, 0, 0);
				}
			}
		}
	}
	//imshow("Clone", newimg);
}

float get_thinness_ratio(int area, int perimeter) {
	return (4 * PI * area) / pow(perimeter, 2);
}

float get_aspect_ratio(Mat *img, uchar R, uchar G, uchar B) {
	int rmin = (*img).rows;
	int rmax = 0;
	int cmin = (*img).cols;
	int cmax = 0;
	for (int i = 0; i < (*img).rows; i++) {
		for (int j = 0; j < (*img).cols; j++) {
			uchar imgR = (*img).at<Vec3b>(i, j)[2];
			uchar imgG = (*img).at<Vec3b>(i, j)[1];
			uchar imgB = (*img).at<Vec3b>(i, j)[0];
			if (imgR == R && imgG == G && imgB == B) {
				if (i < rmin) { rmin = i; }
				if (i > rmax) { rmax = i; }
				if (j < cmin) { cmin = j; }
				if (j > cmax) { cmax = j; }
			}
		}
	}
	float asp_ratio = (float)(cmax - cmin + 1) / (rmax - rmin + 1);
	return asp_ratio;
}

void object_info_callback(int event, int x, int y, int flags, void* param)
{
	//More examples: http://opencvexamples.blogspot.com/2014/01/detect-mouse-clicks-and-moves-on-image.html
	Mat* src = (Mat*)param;
	//Mat clone((*src).rows, (*src).cols, CV_8UC3);
	if (event == CV_EVENT_LBUTTONDOWN)
	{
		Mat newimg = (*src).clone();
		uchar R = (*src).at<Vec3b>(y, x)[2];
		uchar G = (*src).at<Vec3b>(y, x)[1];
		uchar B = (*src).at<Vec3b>(y, x)[0];
		printf("Pos(x,y): %d,%d  Color(RGB): %d,%d,%d\n", x, y, R, G, B);
		int area = get_obj_area(src, R, G, B);
		printf("Object area is: %d\n", area);
		Mat mass_center(src->rows, src->cols, CV_8UC3);
		int *center_of_mass;
		center_of_mass = get_obj_center_of_mass(src, R, G, B);
		int center_r = center_of_mass[0];
		int center_c = center_of_mass[1];
		printf("Rows %d cols %d\n", newimg.rows, newimg.cols);
		printf("Object center of mass: row %d, column %d\n", center_r, center_c);
		double angle = get_axis_of_elongation(src, R, G, B);
		printf("Angle of the elongation axis: %lf\n", angle);
		int perimeter = get_perimeter(src, R, G, B);
		printf("Perimeter: %d\n", perimeter);
		printf("Thinness ratio: %.2f\n", get_thinness_ratio(area, perimeter));
		printf("Aspect ratio: %.2f\n", get_aspect_ratio(src, R, G, B));

		// Draw contour
		draw_contour(src, R, G, B, &newimg);
		newimg.at<Vec3b>(center_r, center_c) = Vec3b(0, 0, 0);
		Point line_start = Point(center_c, center_r);
		float rad_angle = ((angle + 180) * PI) / 180;
		int endx = line_start.x + sin(rad_angle) * 50;
		int endy = line_start.y + cos(rad_angle) * 50;
		Point line_end = Point(endy, endx);
		line(newimg, line_start, line_end, 1, 1);
		imshow("Clone", newimg);
	}
	
}

void get_object_info()
{
	Mat src;
	// Read image from file 
	char fname[MAX_PATH];
	while (openFileDlg(fname))
	{
		src = imread(fname);
		//Create a window
		namedWindow("Objects", 1);

		//set the callback function for any mouse event
		setMouseCallback("Objects", object_info_callback, &src);

		//show the image
		imshow("Objects", src);

		// Wait until user press some key
		waitKey(0);
	}
}

std::vector<Point2i> get_neighbors(int img_height, int img_width, int i, int j) {
	std::vector<Point2i> neighbors;
	neighbors.push_back(Point2i(i - 1, j - 1));
	neighbors.push_back(Point2i(i - 1, j));
	neighbors.push_back(Point2i(i - 1, j + 1));
	neighbors.push_back(Point2i(i, j - 1));
	neighbors.push_back(Point2i(i, j + 1));
	neighbors.push_back(Point2i(i + 1, j - 1));
	neighbors.push_back(Point2i(i + 1, j));
	neighbors.push_back(Point2i(i + 1, j + 1));
	int index = 0;
	for (Point2i p: neighbors) {
		if (p.x < 0 || p.y < 0 || p.x >= img_height || p.y >= img_width) {
			neighbors.erase(neighbors.begin() + index);
		}
		else {
			i++;
		}
	}
	return neighbors;
}

void BFS_labeling() {
	Mat img = imread("Images/labeling2.bmp",
		CV_LOAD_IMAGE_GRAYSCALE);
	imshow("Original", img);
	Mat newimg(img.rows, img.cols, CV_8UC3);
	uchar label = 0;
	uchar R = 0;
	uchar G = 0;
	uchar B = 0;

	Mat labels(img.rows, img.cols, CV_8UC1);
	for (int i = 0; i < img.rows; i++) {
		for (int j = 0; j < img.cols; j++) {
			labels.at<uchar>(i, j) = 0;
		}
	}
	
	for (int i = 0; i < img.rows; i++) {
		for (int j = 0; j < img.cols; j++) {
			if (img.at<uchar>(i, j) == 0 && labels.at<uchar>(i, j) == 0) {
				label++;
				std::queue<Point2i> q;
				labels.at<uchar>(i, j) = label;
				newimg.at<Vec3b>(i, j)[0] = B;
				newimg.at<Vec3b>(i, j)[1] = G;
				newimg.at<Vec3b>(i, j)[2] = R;

				q.push(Point2i(i, j));
				while (!q.empty()) {
					Point2i p = q.front();
					q.pop();
					std::vector<Point2i> neighbors = get_neighbors(img.rows, img.cols, p.x, p.y);
					for (Point2i p: neighbors) {
						int x = p.x;
						int y = p.y;
						if (img.at<uchar>(x, y) == 0 && labels.at<uchar>(x, y) == 0) {
							labels.at<uchar>(x, y) = label;
							q.push(p);
						}
					}
				}
			}
		}
	}

	Vec3b colors[50] = { Vec3b(0,0,0) };
	for (int i = 0; i < labels.rows; i++) {
		for (int j = 0; j < labels.cols; j++) {
			int label = labels.at<uchar>(i, j);
			if (label == 0) {
				newimg.at<Vec3b>(i, j) = Vec3b(255, 255, 255);
			}
			else {
				if (colors[label] == Vec3b(0, 0, 0)) {
					colors[label][0] = rand() % 256;
					colors[label][1] = rand() % 256;
					colors[label][2] = rand() % 256;
					printf("label %d %d %d %d\n", label, colors[label][0], colors[label][1], colors[label][2]);
					

				}
				newimg.at<Vec3b>(i, j) = colors[label];
			}
		}
	}

	imshow("Labels", newimg);
	waitKey(0);
}

std::vector<Point2i> get_prev_neighb(int img_height, int img_width, int i, int j) {
	std::vector<Point2i> neighbors;
	neighbors.push_back(Point2i(i - 1, j - 1));
	neighbors.push_back(Point2i(i - 1, j));
	neighbors.push_back(Point2i(i - 1, j + 1));
	neighbors.push_back(Point2i(i, j - 1));
	int index = 0;

	std::vector<Point2i>::iterator it = neighbors.begin();
	while (it != neighbors.end()) {
		Point2i p = *it;
		if (p.x < 0 || p.y < 0 || p.x >= img_height || p.y >= img_width) {
			it = neighbors.erase(it);
		}
		else ++it;
	}
	return neighbors;
}

void two_pass_label() {
	Mat img = imread("Images/labeling2.bmp",
		CV_LOAD_IMAGE_GRAYSCALE);
	imshow("Original", img);

	Mat newimg(img.rows, img.cols, CV_8UC3);
	uchar label = 0;
	Mat labels(img.rows, img.cols, CV_8UC1);

	for (int i = 0; i < img.rows; i++) {
		for (int j = 0; j < img.cols; j++) {
			labels.at<uchar>(i, j) = 0;
		}
	}

	int nrAccess = 0;

	std::vector<std::vector<int>> edges(1000);
	for (int i = 0; i < img.rows; i++) {
		for (int j = 0; j < img.cols; j++) {
			if (img.at<uchar>(i, j) == 0 && labels.at<uchar>(i, j) == 0) {
				std::vector<uchar> L;
				std::vector<Point2i> prev_n = get_prev_neighb(img.rows, img.cols, i, j);
				for (Point2i n : prev_n) {
					if (labels.at<uchar>(n.x, n.y) > 0) {
						L.push_back(labels.at<uchar>(n.x, n.y));
					}
				}
				if (L.size() == 0) {
					label++;
					labels.at<uchar>(i, j) = label;
				}
				else {
					uchar x = min_element(L.begin(), L.end())[0];
					labels.at<uchar>(i, j) = x;
					for (uchar y : L) {
						if (x != y) {
							//printf("nr access: %d\n", nrAccess++);
							//printf("x: %d\n", x);
							//printf("y: %d\n", y);
							edges[x].push_back(y);
							edges[y].push_back(x);
						}
					}
				}
				
			}
		}
	}
	
	uchar newlabel = 0;
	uchar *newlabels = (uchar*)malloc(img.rows * sizeof(uchar));
	for (int i = 0; i < img.rows; i++) {
		newlabels[i] = 0;
		
	}

	for (int i = 1; i <= label; i++) {
		if (newlabels[i] == 0) {
			newlabel++;
			std::queue<uchar> Q;
			newlabels[i] = newlabel;
			Q.push(i);
			while (!Q.empty()) {
				int x = Q.front();
				Q.pop();
				for (int y : edges[x]) {
					if (newlabels[y] == 0) {
						newlabels[y] = newlabel;
						Q.push(y);
					}
				}
			}
		}
	}

	for (int i = 0; i < img.rows; i++) {
		for (int j = 0; j < img.cols; j++) {
			labels.at<uchar>(i, j) = newlabels[labels.at<uchar>(i, j)];
		}
	}

	imshow("labels", labels);

	Vec3b colors[50] = { Vec3b(0,0,0) };
	for (int i = 0; i < labels.rows; i++) {
		for (int j = 0; j < labels.cols; j++) {
			int label = labels.at<uchar>(i, j);
			if (label == 0) {
				newimg.at<Vec3b>(i, j) = Vec3b(255, 255, 255);
			}
			else {
				if (colors[label] == Vec3b(0, 0, 0)) {
					colors[label][0] = rand() % 256;
					colors[label][1] = rand() % 256;
					colors[label][2] = rand() % 256;
					printf("label %d %d %d %d\n", label, colors[label][0], colors[label][1], colors[label][2]);


				}
				newimg.at<Vec3b>(i, j) = colors[label];
			}
		}
	}
	imshow("New image: ", newimg);
	
	waitKey(0);
}

Point2i getNextPixel(int dir, int i, int j) {
	Point2i nextP;
	switch (dir) {
	case 0: nextP.x = i;
		nextP.y = j + 1;
		break;
	case 1:
		nextP.x = i - 1;
		nextP.y = j + 1;
		break;
	case 2: nextP.x = i - 1;
		nextP.y = j;
		break;
	case 3: nextP.x = i - 1;
		nextP.y = j - 1;
		break;
	case 4: nextP.x = i;
		nextP.y = j - 1;
		break;
	case 5: nextP.x = i + 1;
		nextP.y = j - 1;
		break;
	case 6: nextP.x = i + 1;
		nextP.y = j;
		break;
	case 7: nextP.x = i + 1;
		nextP.y = j + 1;
		break;
	}
	return nextP;
}
void border_tracing() {
	Mat src;
	// Read image from file 
	char fname[MAX_PATH];
	while (openFileDlg(fname))
	{
		src = imread(fname, CV_LOAD_IMAGE_GRAYSCALE);
		//Create a window
		namedWindow("Objects", 1);

		//show the image
		imshow("Objects", src);

		//border tracing algorithm

		//find start point of object
		Point2i startPoint;
		for (int i = 0; i < src.rows; i++) {
			for (int j = 0; j < src.cols; j++) {
				if (src.at<uchar>(i, j) == 0) {
					startPoint = Point2i(i, j);
					i = src.rows;
					j = src.cols;
				}
			}
		}

		printf("starting point: %d, %d\n", startPoint.x, startPoint.y);
		int dir = 7;
		int n = 1;
		std::vector<Point2i> border;
		std::vector<int> chainCode;
		Point2i currentPoint = startPoint;
		border.push_back(startPoint);
		do {
			if (dir % 2 == 0) {
				dir = (dir + +7) % 8;
			}
			else {
				dir = (dir + 6) % 8;
			}
			
			Point2i nextP = getNextPixel(dir, currentPoint.x, currentPoint.y);

			//what happens if no neighbors are black?
			while (src.at<uchar>(nextP.x, nextP.y) != 0) {
				dir = (dir + 1) % 8;
				nextP = getNextPixel(dir, currentPoint.x, currentPoint.y);
			}
			
			if (src.at<uchar>(nextP.x, nextP.y) == 0) {
				border.push_back(nextP);
				n = border.size();
				currentPoint = nextP;
				chainCode.push_back(dir);
			}
		} while (border.size() <= 2 || (border.at(0) != border.at(n - 2) && border.at(1) != border.at(n - 1)));

		//build new image
		Mat borderImg(src.rows, src.cols, CV_8UC1);
		//init with black
		for (int i = 0; i < borderImg.rows; i++) {
			for (int j = 0; j < borderImg.cols; j++) {
				borderImg.at<uchar>(i, j) = 0;
			}
		}
		for (int k = 0; k < border.size(); k++) {
			Point2i p = border.at(k);
			borderImg.at<uchar>(p.x, p.y) = 255;
		}

		printf("Chain code: \n", chainCode.size());
		for (int i = 0; i < chainCode.size(); i++) {
			printf("%d ", chainCode.at(i));
		}

		printf("Derivative chain code: \n");
		for (int i = 1; i < chainCode.size(); i++) {
			int val = (chainCode.at(i) - chainCode.at(i - 1) + 8) % 8;
			printf("%d ", val);
		}

		imshow("Border image", borderImg);

		waitKey(0);
	}
}

void reconstruct_image() {
	ifstream f;
	f.open("Images/reconstruct.txt");
	if (!f) {
		cerr << "Unable to open file reconstruct.txt";
		exit(1);   // call system to stop
	}

	int x, y;
	f >> x;
	f >> y;
	Point2i startPoint = Point2i(x, y);
	int nrChainCode;
	f >> nrChainCode;
	vector<int> chainCode;
	for (int i = 0; i < nrChainCode; i++) {
		f >> x;
		chainCode.push_back(x);
	}
	f.close();

	Mat img = imread("Images/gray_background.bmp",
		CV_LOAD_IMAGE_GRAYSCALE);
	img.at<uchar>(startPoint.x, startPoint.y) = 0;
	Point2i currentPoint = startPoint;
	for (int k = 0; k < chainCode.size(); k++) {
		int dir = chainCode.at(k);
		Point2i nextP = getNextPixel(dir, currentPoint.x, currentPoint.y);
		img.at<uchar>(currentPoint.x, currentPoint.y) = 0;
		currentPoint = nextP;
	}
	imshow("Newimg", img);
	waitKey(0);
}

Mat dilation(Mat* img) {
	// image dilation
	Mat copy;
	img->copyTo(copy);
	
	for (int i = 1; i < img->rows - 1; i++) {
		for (int j = 1; j < img->cols - 1; j++) {
			if (img->at<uchar>(i, j) == 0) {
				copy.at<uchar>(i - 1, j - 1) = 0;
				copy.at<uchar>(i - 1, j) = 0;
				copy.at<uchar>(i - 1, j + 1) = 0;
				copy.at<uchar>(i, j - 1) = 0;
				copy.at<uchar>(i, j + 1) = 0;
				copy.at<uchar>(i + 1, j - 1) = 0;
				copy.at<uchar>(i + 1, j) = 0;
				copy.at<uchar>(i + 1, j + 1) = 0;
			}
		}
	}
	return copy;
}

Mat erosion(Mat* img) {
	// image erosion
	Mat copy;
	img->copyTo(copy);

	for (int i = 1; i < img->rows - 1; i++) {
		for (int j = 1; j < img->cols - 1; j++) {
			if (img->at<uchar>(i, j) == 0) {
				if (img->at<uchar>(i - 1, j - 1) == 255 ||
					img->at<uchar>(i - 1, j) == 255 ||
					img->at<uchar>(i - 1, j + 1) == 255 ||
					img->at<uchar>(i, j - 1) == 255 ||
					img->at<uchar>(i, j + 1) == 255 ||
					img->at<uchar>(i + 1, j - 1) == 255 ||
					img->at<uchar>(i + 1, j) == 255 ||
					img->at<uchar>(i + 1, j + 1) == 255) {
					copy.at<uchar>(i, j) = 255;
				}
			}
		}
	}
	return copy;
}

void dilate_ntimes(int n) {
	Mat img = imread("Images/Morphological_Op_Images/1_Dilate/reg1neg1_bw.bmp", CV_LOAD_IMAGE_GRAYSCALE);
	
	Mat copy;
	img.copyTo(copy);

	//show the image
	imshow("Objects", img);

	for (int i = 0; i < n; i++) {
		copy = dilation(&copy);
	}
	imshow("Copy", copy);
	waitKey(0);
}

void erode_ntimes(int n) {
	Mat img = imread("Images/Morphological_Op_Images/1_Dilate/reg1neg1_bw.bmp", CV_LOAD_IMAGE_GRAYSCALE);
	//Mat img = imread("Images/Morphological_Op_Images/1_Dilate/mon1thr1_bw.bmp", CV_LOAD_IMAGE_GRAYSCALE);
	Mat copy;
	img.copyTo(copy);

	//show the image
	imshow("Objects", img);

	for (int i = 0; i < n; i++) {
		copy = erosion(&copy);
	}
	imshow("Copy", copy);
	waitKey(0);
}

void opening() {
	Mat img = imread("Images/Morphological_Op_Images/3_Open/cel4thr3_bw.bmp", CV_LOAD_IMAGE_GRAYSCALE);
	//erosion followed by dilation
	Mat copy;
	img.copyTo(copy);
	imshow("Original", img);
	copy = erosion(&copy);
	copy = dilation(&copy);

	imshow("Copy", copy);
	waitKey(0);
}

void closing() {
	Mat img = imread("Images/Morphological_Op_Images/4_Close/phn1thr1_bw.bmp", CV_LOAD_IMAGE_GRAYSCALE);
	//erosion followed by dilation
	Mat copy;
	img.copyTo(copy);
	imshow("Original", img);
	copy = dilation(&copy);
	copy = erosion(&copy);

	imshow("Copy", copy);
	waitKey(0);
}

void boundary_extraction() {
	Mat img = imread("Images/Morphological_Op_Images/5_BoundaryExtraction/wdg2thr3_bw.bmp", CV_LOAD_IMAGE_GRAYSCALE);
	imshow("Original", img);
	Mat copy;
	img.copyTo(copy);
	copy = erosion(&copy);
	Mat boundary(img.rows, img.cols, CV_8UC1);
	boundary = img + (255 - copy);
	imshow("Boundary", boundary);
	waitKey(0);
}

Mat inters(Mat* mat1, Mat* mat2) {
	Mat inters(mat1->rows, mat1->cols, CV_8UC1);
	for (int i = 0; i < mat1->rows; i++) {
		for (int j = 0; j < mat1->cols; j++) {
			if (mat1->at<uchar>(i, j) == mat2->at<uchar>(i, j) && mat1->at<uchar>(i,j) == 0) {
				inters.at<uchar>(i, j) = 0;
			}
			else {
				inters.at<uchar>(i, j) = 255;
			}
		}
	}
	return inters;
}

Mat unionImg(Mat* mat1, Mat* mat2) {
	Mat inters(mat1->rows, mat1->cols, CV_8UC1);
	for (int i = 0; i < mat1->rows; i++) {
		for (int j = 0; j < mat1->cols; j++) {
			if (mat1->at<uchar>(i, j) = 0 || mat2->at<uchar>(i, j) == 0) {
				inters.at<uchar>(i, j) = 0;
			}
			else {
				inters.at<uchar>(i, j) = 255;
			}
		}
	}
	return inters;
}

bool areEqual(Mat *mat1, Mat *mat2) {
	for (int i = 0; i < mat1->rows; i++) {
		for (int j = 0; j < mat1->cols; j++) {
			if (mat1->at<uchar>(i, j) != mat2->at<uchar>(i, j)) {
				return false;
			}
		}
	}
	return true;
}

void region_filling() {
	Mat img = imread("Images/Morphological_Op_Images/6_RegionFilling/reg1neg1_bw.bmp", CV_LOAD_IMAGE_GRAYSCALE);
	imshow("Original", img);
	Mat copy(img.rows, img.cols, CV_8UC1);
	Mat complement;
	complement = (255 - img);
	
	int starti = img.rows / 2;
	int startj = img.cols / 2;
	printf("Starting region fill at %d %d\n", starti, startj);
	for (int i = 0; i < copy.rows; i++) {
		for (int j = 0; j < copy.cols; j++) {
			copy.at<uchar>(i, j) = 255;
		}
	}
	copy.at<uchar>(starti, startj) = 0;
	
	Mat copyNew;
	copyNew = dilation(&copy);
	copyNew = inters(&copyNew, &complement);
	int i = 0;
	while (!areEqual(&copyNew, &copy)) {
		copy = copyNew;
		copyNew = dilation(&copy);
		copyNew = inters(&copyNew, &complement);
	}
	copy = unionImg(&copy, &img);
	imshow("Copy", copy);
	waitKey(0);
}

std::vector<int> compute_histogram(Mat *img) {
	std::vector<int> histogram = std::vector<int>(256);
	//initialize vec
	for (int i = 0; i < 256; i++) {
		histogram.at(i) = 0;
	}

	for (int i = 0; i < img->rows; i++) {
		for (int j = 0; j < img->cols; j++) {
			uchar val = img->at<uchar>(i, j);
			histogram.at(val)++;
		}
	}
	return histogram;
}

float compute_mean(int nr_pixels, std::vector<int> hist) {
	float mean = 0;
	for (int g = 0; g < 256; g++) {
		mean += g * hist.at(g);
	}
	return mean / nr_pixels;
}

float standard_deviation(int nr_pixels, float mean, std::vector<int> hist) {
	float std = 0;
	for (int g = 0; g < 256; g++) {
		float pdf = float(hist.at(g)) / nr_pixels;
		std += pow(g - mean, 2) * pdf;
	}
	return sqrt(std);
}

void mean_std() {
	Mat src;
	// Read image from file 
	char fname[MAX_PATH];
	while (openFileDlg(fname))
	{
		src = imread(fname, CV_LOAD_IMAGE_GRAYSCALE);
		//Create a window
		namedWindow("Objects", 1);

		//show the image
		imshow("Objects", src);

		//compute histogram
		std::vector<int> hist = compute_histogram(&src);
		
		int nr_pixels = src.rows * src.cols;
		float mean = compute_mean(nr_pixels, hist);
		printf("Mean: %f\n", mean);
		float std = standard_deviation(nr_pixels, mean, hist);
		printf("Standard dev: %f\n", std);

		waitKey(0);
	}
}

void show_histogram() {
	Mat src;
	// Read image from file 
	char fname[MAX_PATH];
	while (openFileDlg(fname))
	{
		src = imread(fname, CV_LOAD_IMAGE_GRAYSCALE);
		//Create a window
		namedWindow("Objects", 1);

		//show the image
		imshow("Objects", src);

		//compute histogram
		std::vector<int> hist = compute_histogram(&src);
		
		int* h = &hist[0];

		showHistogram("Histogram", h, 256, 200);

		int nr_pixels = src.rows * src.cols;

		std::vector<float> pdf = std::vector<float>(256);
		for (int i = 0; i < 256; i++) {
			pdf.at(i) = float(hist.at(i)) / nr_pixels;
		}

		waitKey(0);
	}
}

void basic_global_threshold() {
	Mat src;
	// Read image from file 
	char fname[MAX_PATH];
	while (openFileDlg(fname))
	{
		src = imread(fname, CV_LOAD_IMAGE_GRAYSCALE);
		//Create a window
		namedWindow("Original", 1);

		//show the image
		imshow("Original", src);

		//compute histogram
		std::vector<int> hist = compute_histogram(&src);

		int* h = &hist[0];
		showHistogram("Histogram", h, 256, 200);

		int min = 255;
		int max = 0;
		for (int i = 0; i < src.rows; i++) {
			for (int j = 0; j < src.cols; j++) {
				uchar val = src.at<uchar>(i, j);
				if (val < min) {
					min = val;
				}
				if (val > max) {
					max = val;
				}
			}
		}

		printf("Min val: %d\n", min);
		printf("Max val: %d\n", max);

		
		float t = (max + min) / 2;
		float t1 = t;
		float error = 0.1;

		do {
			t = t1;
			int n1 = 0;
			float m1 = 0;
			for (int i = min; i <= t; i++) {
				n1 += hist.at(i);
				m1 += i * hist.at(i);
			}

			m1 = m1 / n1;

			int n2 = 0;
			float m2 = 0;
			for (int i = t + 1; i <= max; i++) {
				n2 += hist.at(i);
				m2 += i * hist.at(i);
			}
			m2 = m2 / n2;
			
			t1 = (m1 + m2) / 2;
		} while (abs(t1 - t) >= error);

		printf("threshold is: %f\n", t1);

		for (int i = 0; i < src.rows; i++) {
			for (int j = 0; j < src.cols; j++) {
				if (src.at<uchar>(i, j) < t1) {
					src.at<uchar>(i, j) = 0;
				}
				else {
					src.at<uchar>(i, j) = 255;
				}
			}
		}

		imshow("Thresholding", src);

		waitKey(0);
	}
}

void hist_stretch_shrink(int goutmin, int goutmax) {
	Mat src;
	// Read image from file 
	char fname[MAX_PATH];
	while (openFileDlg(fname))
	{
		src = imread(fname, CV_LOAD_IMAGE_GRAYSCALE);
		//Create a window
		namedWindow("Photo", 1);

		//show the image
		imshow("Photo", src);

		//compute histogram
		std::vector<int> hist = compute_histogram(&src);

		int* h = &hist[0];
		showHistogram("Histogram", h, 256, 200);

		int ginmin = 0;
		int ginmax = 0;

		for (int i = 0; i < 256; i++) {
			if (hist.at(i) != 0) {
				ginmin = i;
				break;
			}
		}
		for (int i = 255; i >= 0; i--) {
			if (hist.at(i) != 0) {
				ginmax = i;
				break;
			}
		}

		printf("G in min: %d\n", ginmin);
		printf("G in max: %d\n", ginmax);

		
		std::vector<int> modif_hist = std::vector<int>(256);
		for (int i = 0; i < 256; i++) {
			modif_hist.at(i) = 0;
		}

		Mat newimg(src.rows, src.cols, CV_8UC1);

		for (int i = 0; i < src.rows; i++) {
			for (int j = 0; j < src.cols; j++) {
				int val = src.at<uchar>(i, j);
				int gout = goutmin + (val - ginmin) * (goutmax - goutmin) / (ginmax - ginmin);
				if (gout > 255) {
					gout = 255;
				}
				else if (gout < 0) {
					gout = 0;
				}
				modif_hist.at(gout) = hist.at(val);
				newimg.at<uchar>(i, j) = gout;
			}
		}

		int* nh = &modif_hist[0];
		showHistogram("New histogram", nh, 256, 200);
		imshow("New photo", newimg);

		waitKey(0);
	}
}

void gamma_correction(float gamma) {
	Mat src;
	// Read image from file 
	char fname[MAX_PATH];
	while (openFileDlg(fname))
	{
		src = imread(fname, CV_LOAD_IMAGE_GRAYSCALE);
		//Create a window
		namedWindow("Photo", 1);

		//show the image
		imshow("Photo", src);

		Mat newimg;
		src.copyTo(newimg);

		for (int i = 0; i < src.rows; i++) {
			for (int j = 0; j < src.cols; j++) {
				int val = src.at<uchar>(i, j);
				float L = 255.0f;
				int gout = L * pow((val / L), gamma);
				//printf("gout is: %f\n", gout);
				if (gout > 255) {
					gout = 255;
				}
				else if (gout < 0) {
					gout = 0;
				}
				newimg.at<uchar>(i, j) = gout;
			}
		}
		imshow("New photo", newimg);

		waitKey(0);
	}
}

void hist_equalization() {
	Mat src;
	// Read image from file 
	char fname[MAX_PATH];
	while (openFileDlg(fname))
	{
		src = imread(fname, CV_LOAD_IMAGE_GRAYSCALE);
		//Create a window
		namedWindow("Objects", 1);

		//show the image
		imshow("Objects", src);

		//compute histogram
		std::vector<int> hist = compute_histogram(&src);

		int* h = &hist[0];
		showHistogram("Histogram", h, 256, 200);

		std::vector<int> cumulative_hist = std::vector<int>(256);
		int s = 0;
		for (int i = 0; i < 256; i++) {
			s += hist.at(i);
			cumulative_hist.at(i) = s;
		}

		int* ch = &cumulative_hist[0];
		showHistogram("Cumulative Histogram", ch, 256, 200);

		int nr_pixels = src.rows * src.cols;
		std::vector<float> cpdf = std::vector<float>(256);
		for (int i = 0; i < 256; i++) {
			cpdf.at(i) = float(cumulative_hist.at(i)) / nr_pixels;
		}

		Mat newimg;
		src.copyTo(newimg);
		for (int i = 0; i < src.rows; i++) {
			for (int j = 0; j < src.cols; j++) {
				int val = src.at<uchar>(i, j);
				int gout = 255 * cpdf.at(val);
				newimg.at<uchar>(i, j) = gout;
			}
		}

		imshow("New image", newimg);
		waitKey(0);
	}
}

int main()
{
	int op;
	do
	{
		system("cls");
		destroyAllWindows();
		printf("Menu:\n");
		printf(" 1 - Open image\n");
		printf(" 2 - Open BMP images from folder\n");
		printf(" 3 - Resize image\n");
		printf(" 4 - Process video\n");
		printf(" 5 - Snap frame from live video\n");
		printf(" 6 - Mouse callback demo\n");
		// Lab 1
		printf(" 7 - L1 Negative Image \n");
		printf(" 8 - Change gray levels by additive value\n");
		printf(" 9 - Change gray levels by multiplicative value\n");
		printf(" 10 - Create color image\n");
		printf(" 11 - Horizontal flip\n");
		printf(" 12 - Vertical flip \n");
		printf(" 13 - Center crop\n");
		printf(" 14 - Resize image\n");
		// Lab 2
		printf(" 15 - Separate RGB\n");
		printf(" 16 - RGB to grayscale\n");
		printf(" 17 - Grayscale to binary\n");
		printf(" 18 - RGB to HSV\n");
		printf(" 19 - Detect traffic sign\n");
		printf(" 20 - Segment by RGB\n"); //TODO
		// Lab 3
		printf(" 21 - Object information\n");
		//Lab 4
		printf(" 22 - BFS labeling\n");
		printf(" 23 - Two pass labeling\n");
		// Lab 5
		printf(" 24 - Border tracing (with chain codes and derivative chain code\n");
		printf(" 25 - Reconstruct image\n");

		// Lab 6
		printf(" 26 - Dilation\n");
		printf(" 27 - Erosion\n");
		printf(" 28 - Open\n");
		printf(" 29 - Close\n");
		printf(" 30 - Boundary extraction\n");
		printf(" 31 - Region filling\n");

		// Lab 7
		printf(" 32 - Mean and standard deviation\n");
		printf(" 33 - Show histogram\n");
		printf(" 34 - Global thresholding\n");
		printf(" 35 - Histogram stretch/shrink\n");
		printf(" 36 - Gamma correction\n");
		printf(" 37 - Histogram equalization\n");

		printf(" 0 - Exit\n\n");
		printf("Option: ");
		scanf("%d", &op);
		switch (op)
		{
		case 1:
			testOpenImage();
			break;
		case 2:
			testOpenImagesFld();
			break;
		case 3:
			testResize();
			break;
		case 4:
			testVideoSequence();
			break;
		case 5:
			testSnap();
			break;
		case 6:
			testMouseClick();
			break;
		case 7:
			negative_image();
			break;
		case 8:
			gray_levels_additive();
			break;
		case 9:
			gray_levels_multiplicative();
			break;
		case 10:
			create_color_image();
			break;
		case 11:
			horizontal_flip();
			break;
		case 12:
			vertical_flip();
			break;
		case 13:
			printf("Input crop size:\n");
			int height, width;
			scanf("%d %d", &height, &width);
			center_crop(height, width);
			break;
		case 14:
			printf("Input new size:\n");
			scanf("%d %d", &height, &width);
			resize_image(height, width);
			break;
		case 15:
			separate_RGB();
			break;
		case 16:
			RGB_to_grayscale();
			break;
		case 17:
			unsigned int threshold;
			while (1) {
				printf("Input threshold (between 0 and 255):\n");
				scanf("%u", &threshold);
				if (threshold >= 0 && threshold <= 255) {
					break;
				}
			}
			grayscale_to_binary(threshold);
			break;
		case 18:
			RGB_to_HSV();
			break;
		case 19:
			detect_traffic_sign();
			break;
		case 20:
			RGB_mask();
			break;
		case 21:
			get_object_info();
			break;
		case 22:
			BFS_labeling();
			break;
		case 23:
			two_pass_label();
			break;
		case 24:
			border_tracing();
			break;
		case 25:
			reconstruct_image();
			break;
		case 26:
			printf("Input number of times to apply dilation:\n");
			int n;
			scanf("%d", &n);
			dilate_ntimes(n);
			break;
		case 27:
			printf("Input number of times to apply erosion:\n");
			scanf("%d", &n);
			erode_ntimes(n);
			break;
		case 28:
			opening();
			break;
		case 29:
			closing();
			break;
		case 30:
			boundary_extraction();
			break;
		case 31:
			region_filling();
			break;
		case 32:
			mean_std();
			break;
		case 33:
			show_histogram();
			break;
		case 34:
			basic_global_threshold();
			break;
		case 35:
			
			int min;
			int max;
			printf("Input min range value:\n");
			scanf("%d", &min);
			printf("Input max range value:\n");
			scanf("%d", &max);
			hist_stretch_shrink(min, max);
			break;
		case 36:
			float gamma;
			printf("Input gamma:\n");
			std::cin >> gamma;
			gamma_correction(gamma);
			break;
		case 37:
			hist_equalization();
			break;
		}
		
	} while (op != 0);
	return 0;
}