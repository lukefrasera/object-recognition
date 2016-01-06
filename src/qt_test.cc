#include "log.h"
#include "quad_tree.h"
#include <stdio.h>
#include <stdlib.h>
#include <opencv2/opencv.hpp>
#include <vector>
int main(int argc, char *argv[]) {
  printf("Testing Quad Tree...\n");
  cv::Mat image(100,100, CV_8UC3, cv::Scalar(0,0,0));
  utils::QuadTree qt(utils::AABB(utils::Point(50.0, 50.0), 50.0), 4);
  std::vector<utils::AABB> bounding_boxes;
  qt.GetBounds(&bounding_boxes);
  printf("number of Bounding Boxes: %d\n", bounding_boxes.size());
  float x = static_cast<float>(rand()) / static_cast<float>(RAND_MAX) * 10.0 - 5.0;
  float y = static_cast<float>(rand()) / static_cast<float>(RAND_MAX) * 10.0 - 5.0;
  qt.Insert(utils::Point(x, y));

  cv::Point x1, x2;
  for (int i = 0; i < bounding_boxes.size(); ++i) {
    x1 = cv::Point(
      bounding_boxes[i].center_.x_ - bounding_boxes[i].half_dimension_,
      bounding_boxes[i].center_.y_ - bounding_boxes[i].half_dimension_);
    printf("Top Left Corner: %d, %d\n", x1.x, x1.y);
    x2 = cv::Point(
      bounding_boxes[i].center_.x_ + bounding_boxes[i].half_dimension_,
      bounding_boxes[i].center_.y_ + bounding_boxes[i].half_dimension_);
    cv::rectangle(image, x1, x2, cv::Scalar(255,255,0), 2);
  }

  cv::imshow("QT_test", image);
  cv::waitKey(0);
  return 0;
}