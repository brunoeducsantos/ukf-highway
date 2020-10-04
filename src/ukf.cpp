#include "ukf.h"
#include <Eigen/Dense>
#include <iostream>
#include <yaml-cpp/yaml.h>
using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF(std::string filename)
{
  YAML::Node config = YAML::LoadFile(filename);
  
  // if this is false, laser measurements will be ignored (except during init)
  // use_laser_ = true;
  use_laser_= config["use_laser"].as<bool>();
  // if this is false, radar measurements will be ignored (except during init)
  // use_radar_ = true;
  use_radar_= config["use_radar"].as<bool>();
  
  // is_initialized_ = false;
  is_initialized_= config["is_initialized"].as<bool>();

  // initial state vector
  x_ = VectorXd(5);
  x_.setZero();
  // initial covariance matrix
  P_ = MatrixXd::Zero(5, 5);
  // Process noise standard deviation longitudinal acceleration in m/s^2
  // std_a_ = 3.5;
  std_a_= config["std_a"].as<double>();
  // Process noise standard deviation yaw acceleration in rad/s^2
  // std_yawdd_ = 0.8;
  std_yawdd_= config["std_yawdd"].as<double>();
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
   * End DO NOT MODIFY section for measurement noise values 
   */

  /**
   *  Initialization. 
   */
  n_aug_ = 7;
  n_x_ = 5;
  
  lambda_ = (3-n_aug_)*1. ;
  Xsig_pred_ = MatrixXd::Zero(n_x_, 2 * n_aug_ + 1);
  // set weights
  weights_ = VectorXd(2 * n_aug_ + 1);
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights_(0) = weight_0;
  for (int i = 1; i < 2 * n_aug_ + 1; ++i)
  { // 2n+1 weights
    double weight = 0.5 /(n_aug_ + lambda_);
    weights_(i) = weight;
  }
  //NIS metric
  epsilon=0.;
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package)
{


  if (!is_initialized_)
  {
    //TODO Fix state initialization
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    { 
      double rho=  meas_package.raw_measurements_(0);
      double phi=  meas_package.raw_measurements_(1);
      x_(0)= rho*cos(phi);
      x_(1)= rho*sin(phi);
      x_.tail(3) = meas_package.raw_measurements_;
      P_(2, 2) = std_radr_ * std_radr_;
      P_(3, 3) = std_radphi_ * std_radphi_;
      P_(4, 4) = std_radrd_ * std_radrd_;
      
    }
    else
    {
      x_.head(2) = meas_package.raw_measurements_;
      P_(0, 0) = std_laspx_ * std_laspx_;
      P_(1, 1) = std_laspy_ * std_laspy_;
      }
    //covariance initialization

    is_initialized_ = true;
    time_us_= meas_package.timestamp_;
    return ;
  }

  else
  {
    
    double dt = (meas_package.timestamp_ - time_us_) / pow(10,6);

    Prediction(dt);
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
    {
      UpdateRadar(meas_package);
    }
    else if (use_laser_ && (meas_package.sensor_type_ == MeasurementPackage::LASER))
    {
      UpdateLidar(meas_package);
    }
    time_us_= meas_package.timestamp_;
  }
}
void UKF::SigmaPoints(double delta_t)
{
  //Create augmented state vector
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.fill(0.0);
  VectorXd x_aug(n_aug_);
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //Augmented covariance matrix
  MatrixXd P_aug = MatrixXd::Zero(7, 7);
  P_aug.topLeftCorner(5, 5) = P_;
  MatrixXd Q(2, 2);
  Q = MatrixXd::Zero(2, 2);
  Q(1, 1) = std_yawdd_ * std_yawdd_;
  Q(0, 0) = std_a_ * std_a_;
  P_aug.bottomRightCorner(2, 2) = Q;

  // create square root matrix
  // Record start time
  MatrixXd L = P_aug.llt().matrixL();
  // create augmented sigma points
  Xsig_aug.col(0) = x_aug;

  for (int i = 0; i < n_aug_; ++i)
  { 
    Xsig_aug.col(i + 1) = x_aug + sqrt((lambda_ + n_aug_)*1.) * L.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt((lambda_ + n_aug_)*1.) * L.col(i);
  }

  // predict sigma points
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    // extract values for better readability
    double p_x = Xsig_aug(0, i);
    double p_y = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);
    // predicted state values
    double px_p, py_p;

    // avoid division by zero
    if (fabs(yawd) > 0.001)
    {
      px_p = p_x + (v / yawd) * (sin(yaw + yawd * delta_t) - sin(yaw));
      py_p = p_y + (v / yawd) * (cos(yaw) - cos(yaw + yawd * delta_t));
    }
    else
    {
      px_p = p_x + v * delta_t * cos(yaw);
      py_p = p_y + v * delta_t * sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    // add noise
    px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    v_p = v_p + nu_a * delta_t;

    yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p = yawd_p + nu_yawdd * delta_t;

    // write predicted sigma point into right column
    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;
  }
}

void UKF::Prediction(double delta_t)
{
  /**   

   *  Predict sigma points, the state, 
   * and the state covariance matrix.
   */
  //Predict sigma points
  SigmaPoints(delta_t);

  // predicted state mean
  VectorXd x_pred(n_x_);
  x_pred.fill(0.);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  { // iterate over sigma points
    x_pred += weights_(i) * Xsig_pred_.col(i);
  }

  // predicted state covariance matrix
  P_.fill(0.);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  { // iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    
    // angle normalization
    while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
   
  x_=x_pred;
}

void UKF::UpdateLidar(MeasurementPackage meas_package)
{
  /**
   *  Lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */

  int n_z = 2;
  Zsig_ = MatrixXd(n_z, 2 * n_aug_ + 1);
  Zsig_.fill(0.0);
  Eigen::VectorXd z(n_z);
  z.fill(0.0);
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);

  // transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  { // 2n+1 simga points
    // extract values for better readability
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);
    Zsig_(0, i) = p_x;
    Zsig_(1, i) = p_y;
  }
  // mean predicted measurement
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    z = z + weights_(i) * Zsig_.col(i);
  }
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  { // 2n+1 simga points
    // residual
    VectorXd z_diff = Zsig_.col(i) - z;

    // angle normalization

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  // add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_laspx_ * std_laspx_, 0,
      0, std_laspy_ * std_laspy_;
  S = S + R;

  Eigen::MatrixXd T(n_x_, n_z);
  T.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

    VectorXd z_diff = Zsig_.col(i) - z;
    // angle normalization

    T = T + weights_(i) * x_diff * z_diff.transpose();
  }
  MatrixXd K = T * S.inverse();

  // angle normalization
  VectorXd z_diff = meas_package.raw_measurements_ - z;
  //Mean update
  x_ = x_ + K * z_diff;
  //Covariance update
  P_ = P_ - (K * S * K.transpose());

  //NIS
  epsilon=z_diff.transpose()*S.inverse()*z_diff;
  //TODO accumulate into vector
}

void UKF::UpdateRadar(MeasurementPackage meas_package)
{
  /**
   * Update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  int n_z = 3;
  Zsig_ = MatrixXd(n_z, 2 * n_aug_ + 1);
  Zsig_.fill(0.0);
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);

  // transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  { // 2n+1 simga points
    // extract values for better readability
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    while (yaw > M_PI) yaw -= 2. * M_PI;
    while (yaw< -M_PI) yaw += 2. * M_PI;

    double v1 = cos(yaw) * v;
    double v2 = sin(yaw) * v;
    // measurement model
    Zsig_(0, i) = sqrt(p_x * p_x + p_y * p_y);                         // r
    Zsig_(1, i) = atan2(p_y, p_x);                                     // phi
    Zsig_(2, i) = (p_x * v1 + p_y * v2) / sqrt(p_x * p_x + p_y * p_y); // r_dot
  }
  // mean predicted measurement
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    z_pred = z_pred + weights_(i) * Zsig_.col(i);
  }

  // innovation covariance matrix S
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  { // 2n+1 simga points
    // residual
    VectorXd z_diff = Zsig_.col(i) - z_pred;

    // angle normalization
    while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  // add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_radr_ * std_radr_, 0, 0,
      0, std_radphi_ * std_radphi_, 0,
      0, 0, std_radrd_ * std_radrd_;
  S = S + R;

  //UKF radar update

  Eigen::MatrixXd T(n_x_, n_z);
  T.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

    VectorXd z_diff = Zsig_.col(i) - z_pred;

    // angle normalization
    while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;
    T = T + weights_(i) * x_diff * z_diff.transpose();
  }
  MatrixXd K = T * S.inverse();

  // angle normalization
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;
  while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
  while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;
  //Mean update
  x_ = x_ + K * z_diff;
  //Covariance update
  P_ = P_ - (K * S * K.transpose());
  //NIS
  epsilon=z_diff.transpose()*S.inverse()*z_diff;
  //TODO accumulate into vector
}