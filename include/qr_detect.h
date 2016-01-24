#ifndef INCLUDE_QR_DETECT_H_
#define INCLUDE_QR_DETECT_H_
#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <vector>

namespace qr {
// Type Definitions //
typedef std::vector< std::vector<cv::Point> > Contour_t;

std::vector<std::vector<cv::Point2f> > QRDetectIdentifiers(cv::Mat image, Contour_t *ids);
}
#endif  // INCLUDE_QR_DETECT_H_
