#include "qr_detect.h"
#include <stdint.h>
#include <map>
#include "quad_tree.h"
#include "log.h"
#include <opencv2/opencv.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/flann/miniflann.hpp>
#include <iostream>
#include <cmath>

namespace qr {
struct ContourMoment {
  std::vector<cv::Point> contour;
  cv::Point2f moment;
};

void QRGetContoursFromImage(cv::Mat image, Contour_t *contours,
    std::vector<cv::Vec4i> *hierarchy) {
  cv::Mat gray(image.size(), CV_MAKETYPE(image.depth(), 1));
  cv::Mat edge(image.size(), CV_MAKETYPE(image.depth(), 1));

  cv::cvtColor(image, gray, CV_RGB2GRAY);
  cv::Canny(gray, edge, 100, 200, 3);

  cv::findContours(edge, *contours, *hierarchy, cv::RETR_TREE,
    cv::CHAIN_APPROX_SIMPLE);
}

void QRGetContoursMassCenters(const Contour_t *contours,
    std::vector<cv::Point2f> *points) {
  std::vector<cv::Moments> mu(contours->size());
  for (int i = 0; i < contours->size(); ++i) {
    mu[i] = moments((*contours)[i], false);
    points->push_back(cv::Point2f(mu[i].m10/mu[i].m00, mu[i].m01/mu[i].m00));
  }
}

void QRCullContourCriteriaLengthRatio(const Contour_t *contours,
    const std::vector<cv::Vec4i> hierarchy,
    int32_t *mask) {
  float xmin, xmax, ymin, ymax;
  float e = 0.5;
  for (int i = 0; i < contours->size(); ++i) {
    if (mask[i]) {
      // Find Xmin, Xmax, Ymin, Ymax
      xmin = (*contours)[i][0].x;xmax = xmin;
      ymin = (*contours)[i][0].y;ymax = ymin;
      for (int j = 1; j < (*contours)[i].size(); ++j) {
        int x = (*contours)[i][j].x;
        int y = (*contours)[i][j].y;
        if (x > xmax)
          xmax = x;
        if (x < xmin)
          xmin = x;
        if (y > ymax)
          ymax = y;
        if (y < ymin)
          ymin = y;
      }
      // Determin Length ratio
      float r = (xmax - xmin)/(ymax - ymin);
      // printf("R: %f\n", r);
      if (!(1 - e < r && r < 1 + e))
        mask[i] = 0;
    }
  }
}

void QRCullContourCriteriaLength(const Contour_t *contours,
    const std::vector<cv::Vec4i> hierarchy,
    int32_t *mask) {
  int t = 25;
  for (int i = 0; i < contours->size(); ++i) {
    if ((*contours)[i].size() < t) {
      mask[i] = 0;
    }
  }
}

uint32_t CountContourChildren(const std::vector<cv::Vec4i> *hierarchy,
    int idx) {
  if ((*hierarchy)[idx][2] == -1) {
    return 0;
  } else {
    return CountContourChildren(hierarchy, (*hierarchy)[idx][2]) + 1;
  }
}

void QRCullContourCriteriaInterior(const Contour_t *contours,
    const std::vector<cv::Vec4i> hierarchy,
    int32_t *mask) {
  for (int i = 0; i < contours->size(); ++i) {
    if (mask[i] == 1) {
      if (CountContourChildren(&hierarchy, i) < 3) {
        mask[i] = 0;
      }
    }
  }
}

void QRCullContourCriteriaOverlap(const Contour_t *contours,
    const std::vector<cv::Vec4i> hierarchy,
    int32_t *mask) {
  float dist = 2;
  cv::Mat features;
  utils::QuadTree qt(
    utils::AABB(utils::Point(1280.0/2.0, 720.0/2), 1280.0/2.0), 4);
  std::vector<cv::Point2f> moments;
  std::vector<cv::Point2f> feature_points;
  // Compute moments
  QRGetContoursMassCenters(contours, &moments);

  for (int i = 0; i < moments.size(); ++i) {
    if (mask[i] == 1) {
      if (!(std::isnan(moments[i].x) || std::isnan(moments[i].y))){
        feature_points.push_back(moments[i]);
      } else {
        mask[i] = 0;
      }
    }
  }
  // features = cv::Mat(feature_points.size(), 2, CV_32F);
  // float *data = (float*)features.data;
  // for (int i = 0; i < feature_points.size(); ++i) {
  //   data[2*i]     = feature_points[i].x;
  //   data[2*i + 1] = feature_points[i].y;
  // }

  // std::cout << features << std::endl;
  // cv::flann::KDTreeIndexParams params;
  // // Create Momoments matrix as samples
  // cv::flann::Index knn_object(features, params, cvflann::FLANN_DIST_EUCLIDEAN);

  // std::vector<std::vector<float> > range(5);
  // std::vector<std::vector<float> > distance(5);
  // std::vector<float> query(2);
  // for (int i = 0; i < moments.size(); ++i) {
  //   if (mask[i] == 1) {
  //     query[0] = moments[i].x;
  //     query[1] = moments[i].y;
  //     knn_object.radiusSearch(query, range, distance, 3.0, 5);
  //     std::cout << range[0][0] << std::endl;
  //   }
  // }

  for (int i = 0; i < contours->size(); ++i) {
    if (mask[i] == 1) {
      qt.Insert(utils::Point(moments[i]));
    }
  }
  // Find Nearest Neghbors to discover identification markers
  std::vector<utils::Point> range;
  for (int i = 0; i < contours->size(); ++i) {
    if (mask[i] == 1) {
      range = qt.QueryRange(utils::AABB(utils::Point(moments[i]), dist));
      if (range.size() < 3)
        mask[i] = 0;
    }
  }
}

void QRApplyContourMask(Contour_t *contours, const int32_t *mask) {
  std::vector<std::vector<cv::Point> >::iterator it;
  Contour_t culled;
  int i = 0;
  for (it = contours->begin(); it != contours->end(); ++it) {
    if (mask[i]) {
      culled.push_back(*it);
    }
    i++;
  }
  (*contours) = culled;
}

void QRCullContours(Contour_t *contours, std::vector<cv::Vec4i> hierarchy) {
  int32_t *mask = new int32_t[contours->size()];
  for (int i = 0; i < contours->size(); ++i) {mask[i] = 1;}
  QRCullContourCriteriaLengthRatio(contours, hierarchy, mask);
  QRCullContourCriteriaLength(contours, hierarchy, mask);
  QRCullContourCriteriaOverlap(contours, hierarchy, mask);
  QRCullContourCriteriaInterior(contours, hierarchy, mask);
  QRApplyContourMask(contours, mask);
  delete [] mask;
}

void QRDetectIdentifiers(cv::Mat image, Contour_t *ids) {
  std::vector<cv::Vec4i> hierarchy;

  QRGetContoursFromImage(image, ids, &hierarchy);
  QRCullContours(ids, hierarchy);
}
}