#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <opencv2/opencv.hpp>
#include "log.h"
#include "qr_detect.h"

#define FREE_CHECK_NULL(x) do {                                                \
  if (!x) {                                                                    \
    free(x);                                                                   \
    x = NULL;                                                                  \
  }                                                                            \
} while(0)

int main(int argc, char *argv[]) {


  // Detect QR Corners Identifiers
  qr::Contour_t corners;
  qr::QRDetectIdentifiers(in_img, &corners);

  // cv::imshow("QR Test", in_img);
  cv::waitKey(0);
  FREE_CHECK_NULL(input_img);
  FREE_CHECK_NULL(output_img);
  return 0;
}