#include "qr_detect.h"
#include <stdint.h>
#include <map>
#include "quad_tree.h"
#include "log.h"
#include <opencv2/opencv.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/flann/miniflann.hpp>
#include <opencv2/flann/flann.hpp>
#include <iostream>
#include <cmath>
#include <queue>

namespace qr {
struct ContourMoment {
  std::vector<cv::Point> contour;
  cv::Point2f moment;
};

struct PointCoord {
  PointCoord(int v_, cv::Point p_) {
    value = v_;
    p = p_;
  }
  int value;
  cv::Point p;
};

bool operator>(const PointCoord &a, const PointCoord &c) {
  return a.value > c.value;
}

bool operator>=(const PointCoord &a, const PointCoord &c) {
  return a.value >= c.value;
}

bool operator<(const PointCoord &a, const PointCoord &c) {
  return a.value < c.value;
}

bool operator<=(const PointCoord &a, const PointCoord &c) {
  return a.value <= c.value;
}

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
      if (range.size() < 3) {
        mask[i] = 0;
      }
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

cv::Mat QRGeneratePoints(Contour_t *contours) {
  std::vector<cv::Point2f> points;
  QRGetContoursMassCenters(contours, &points);
  utils::QuadTree qt(
    utils::AABB(utils::Point(640.0/2.0, 480/2.0), 640.0/2.0), 4);

  for (int i = 0; i < points.size(); ++i) {
    if (!std::isnan(points[i].x))
      qt.Insert(utils::Point(points[i]));
  }

  std::vector<utils::Point> range;
  std::vector<cv::Point2f> features;
  for (int i = 0; i < points.size(); ++i) {
    range = qt.QueryRange(utils::AABB(utils::Point(points[i]), 3.0));
    // find mean range
    cv::Point2f avg_point;
    avg_point.x = 0;
    avg_point.y = 0;
    for (int j = 0; j < range.size(); ++j) {
      avg_point.x += range[j].x_;
      avg_point.y += range[j].y_;
    }
    avg_point.x /= range.size();
    avg_point.y /= range.size();

    // Check if point is in features yet
    bool within = false;
    for (int j = 0; j < features.size(); ++j) {
      cv::Point2f diff = features[j] - avg_point;
      float dist = std::sqrt(std::pow(diff.x, 2) + std::pow(diff.y, 2));
      if (dist < 3)
        within = true;
    }
    if (!within) {
      features.push_back(avg_point);
    }
  }
  cv::Mat_<float> output(features.size(), 2);
  float *data = (float*)output.data;
  for (int i = 0; i < features.size(); ++i) {
    data[2*i]    = features[i].x;
    data[2*i +1] = features[i].y;
  }
  return output;
}

cv::Point ComputeIncenter(
    cv::Point2f A,
    cv::Point2f B,
    cv::Point2f C) {
  float a,b,c;
  cv::Point2f ab = A - B;
  cv::Point2f bc = B - C;
  cv::Point2f ca = C - A;
  a = std::sqrt(pow(ab.x,2) + pow(ab.y,2));
  b = std::sqrt(pow(bc.x,2) + pow(bc.y,2));
  c = std::sqrt(pow(ca.x,2) + pow(ca.y,2));

  cv::Point center;
  center.x = (A.x + B.x + C.x) / (a + b + c);
  center.y = (A.y + B.y + C.y) / (a + b + c);

  return center;
}

std::vector<cv::Point> FindNMaxPoints(cv::Mat field, int num) {
  std::priority_queue<PointCoord, std::vector<PointCoord>, std::less<std::vector<PointCoord>::value_type> > queue;
  int *data = (int*)field.data;
  for (int i = 0; i < field.rows; ++i) {
    for (int j = 0; j < field.cols; ++j) {
      queue.push(PointCoord(data[field.cols * i + j], cv::Point(j, i)));
    }
  }
  std::vector<cv::Point> points;
  for (int i = 0; i < num; ++i) {
    points.push_back(queue.top().p);
    queue.pop();
  }
  return points;
}

std::vector<std::vector<cv::Point2f> > QRDetectIdentifiers(cv::Mat image, Contour_t *ids) {
  std::vector<cv::Vec4i> hierarchy;
  std::vector<std::vector<cv::Point2f> > points;
  QRGetContoursFromImage(image, ids, &hierarchy);
  QRCullContours(ids, hierarchy);

  // average contours into single mass centers
  cv::Mat feature_points = QRGeneratePoints(ids);

  if (feature_points.rows > 3) {
    if (!feature_points.empty()) {
      // Perform KNN to put points into groups of three
      int knn = 3;
      cv::Mat_<int> indices(feature_points.rows, knn);
      cv::Mat_<float> dists(feature_points.rows, knn);

      cv::flann::GenericIndex<cvflann::ChiSquareDistance<float> > index(feature_points,
        cvflann::KDTreeIndexParams());

      index.knnSearch(feature_points, indices, dists, knn, cvflann::SearchParams(feature_points.cols-1));

      // group points into 3's
      uint32_t num_groups = feature_points.rows / 3;
      LOG_INFO("HELLO");
      // Vote on Groups based on 3 point triangle incenter
      std::vector<cv::Point> incenter_vector(feature_points.rows);
      cv::Mat_<int> voting_field = cv::Mat::zeros(image.rows, image.cols, CV_32S);
      for (int j = 0; j < feature_points.rows; ++j) {
        // compute incenter
        cv::Point2f A,B,C;
        cv::Point incenter;
        int a_idx, b_idx, c_idx;
        a_idx = indices.at<int>(j, 0);
        b_idx = indices.at<int>(j, 1);
        c_idx = indices.at<int>(j, 2);

        A.x = feature_points.at<float>(a_idx, 0);
        A.y = feature_points.at<float>(a_idx, 1);

        B.x = feature_points.at<float>(b_idx, 0);
        B.y = feature_points.at<float>(b_idx, 1);

        C.x = feature_points.at<float>(c_idx, 0);
        C.y = feature_points.at<float>(c_idx, 1);
        incenter = ComputeIncenter(A, B, C);
        LOG_INFO("Point Centers: %f, %f", incenter.x, incenter.y);
        voting_field.at<uint32_t>(incenter.y, incenter.x) += 1;
        incenter_vector.push_back(incenter);
      }
      LOG_INFO("HELLO");

      std::vector<cv::Point> group_incenters = FindNMaxPoints(voting_field, num_groups);
      points = std::vector<std::vector<cv::Point2f> >(group_incenters.size());
      LOG_INFO("HELLO");
      for (int i = 0; i < group_incenters.size(); ++i) {
        for (int j = 0; j < feature_points.rows; ++j) {
          if (incenter_vector[i] == group_incenters[j]) {
            int a_idx = indices.at<int>(j,0);
            int b_idx = indices.at<int>(j,1);
            int c_idx = indices.at<int>(j,2);

            points[i].push_back(
                                cv::Point2f(
                                  feature_points.at<float>(a_idx, 0),
                                  feature_points.at<float>(a_idx, 1)));
            points[i].push_back(
                                cv::Point2f(
                                  feature_points.at<float>(b_idx, 0),
                                  feature_points.at<float>(b_idx, 1)));
            points[i].push_back(
                                cv::Point2f(
                                  feature_points.at<float>(c_idx, 0),
                                  feature_points.at<float>(c_idx, 1)));
            break;
          }
        }
      }
    }
  } else if (feature_points.rows == 3) {
    if (feature_points.rows > 0) {
      float *data = (float*)feature_points.data;
      std::vector<cv::Point2f> point;
      for (int i = 0; i < feature_points.rows; ++i) {
        point.push_back(cv::Point2f(data[2*i], data[2*i+1]));
      }
      points.push_back(point);
    }
  }
  return points;
}
}