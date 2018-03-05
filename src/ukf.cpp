#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.3;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  is_initialized_ = false;
  P_ << 1, 0, 0, 0, 0,
	   0, 1, 0, 0, 0,
	   0, 0, 1, 0, 0,
	   0, 0, 0, 1, 0,
	   0, 0, 0, 0, 1;

  n_x_ = 5;
  n_aug_ = 7;

  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);
  weights_ = VectorXd(2*n_aug_+1);

  lambda_ = 3 - n_aug_;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
	if (!is_initialized_)
	{
		if (meas_package.sensor_type_ == MeasurementPackage::LASER)
		{
			//initialize state with first measurement
			x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
		}
		else
		{
			//RADAR measurement
			float ro = meas_package.raw_measurements_[0];
			float theta = meas_package.raw_measurements_[1];
			float rodot = meas_package.raw_measurements_[2];
			float v = sqrt(rodot*cos(theta)*rodot*cos(theta) + rodot*sin(theta)*rodot*sin(theta));
			x_ << ro*cos(theta), ro*sin(theta), v, 0, 0;
		}
		is_initialized_ = true;
		time_us_ = meas_package.timestamp_;
		return;
	}
	//call prediction function
	if ((meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) ||
		(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_))
	{
		double delta_t = (meas_package.timestamp_ - time_us_)/1000000.0;
		time_us_ = meas_package.timestamp_;
		Prediction(delta_t);
	}

	//call update function based on type of sensor
	if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
	{
		UpdateLidar(meas_package);
	}
	else if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
	{
		UpdateRadar(meas_package);
	}
}

void UKF::PredictSigmaPoints(VectorXd& x_aug, MatrixXd& P_aug, double delta_t) {
	  //sigma point matrix
	  MatrixXd X_sigma_aug = MatrixXd(n_aug_, 2*n_aug_+1);
	  MatrixXd A = P_aug.llt().matrixL();
	  double mult = (lambda_ + n_aug_);

	  X_sigma_aug.col(0) = x_aug;
	  VectorXd A_col_factor;
	  for (int i = 0; i < n_aug_; i++)
	  {
		  A_col_factor = sqrt(mult)*A.col(i);
		  X_sigma_aug.col(i+1) = x_aug + A_col_factor;
		  X_sigma_aug.col(n_aug_+i+1) = x_aug - A_col_factor;
	  }

	  //Prediction
	  Xsig_pred_.fill(0.0);
	  for (int i = 0; i < 2*n_aug_+1; i++)
	  {
		  float px = X_sigma_aug(0,i);
		  float py = X_sigma_aug(1,i);
		  float v = X_sigma_aug(2,i);
		  float yaw = X_sigma_aug(3,i);
		  float yawd = X_sigma_aug(4,i);
		  float nu_a = X_sigma_aug(5,i);
		  float nu_yawdd = X_sigma_aug(6,i);

		  if (fabs(yawd) > 0.001)
		  {
			  Xsig_pred_(0,i) = px + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
			  Xsig_pred_(1,i) = py + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
		  }
		  else
		  {
			  Xsig_pred_(0,i) = px + v*delta_t*cos(yaw);
			  Xsig_pred_(1,i) = py + v*delta_t*sin(yaw);
		  }
		  Xsig_pred_(0,i) += 0.5*delta_t*delta_t*cos(yaw)*nu_a;
		  Xsig_pred_(1,i) += 0.5*delta_t*delta_t*sin(yaw)*nu_a;
		  Xsig_pred_(2,i) = v + delta_t*nu_a;
		  Xsig_pred_(3,i) = yaw + yawd*delta_t + 0.5*delta_t*delta_t*nu_yawdd;
		  Xsig_pred_(4,i) = yawd + delta_t*nu_yawdd;
	  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_;

  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_,n_x_) = std_a_*std_a_;
  P_aug(n_x_+1,n_x_+1) = std_yawdd_*std_yawdd_;

  PredictSigmaPoints(x_aug, P_aug, delta_t); //populates Xsig_pred_

  //Predict state mean and covariance
  weights_(0) = lambda_/(lambda_+n_aug_);
  for (int i = 1; i < 2*n_aug_+1; i++)
  {
	  weights_(i) = 0.5/(lambda_+n_aug_);
  }

  	  x_.fill(0.0);
	  for (int i = 0; i < 2*n_aug_+1; i++)
	  {
		  x_ = x_ + weights_(i)*Xsig_pred_.col(i);
	  }

	  P_.fill(0.0);
	  for (int i = 0; i < 2*n_aug_+1; i++)
	  {
		  VectorXd xdiff = Xsig_pred_.col(i) - x_;
		  if (xdiff(3) > M_PI) xdiff(3) -= 2*M_PI;
		  else if (xdiff(3) < -M_PI) xdiff(3) += 2*M_PI;

		  P_ = P_ + weights_(i)*xdiff*xdiff.transpose();
	  }

  }

void UKF::measMeanAndCovariance(MatrixXd Zsig, MatrixXd R, VectorXd& z_pred, MatrixXd& S, bool isRadar) {
	z_pred.fill(0.0);
	S.fill(0.0);

	for (int i = 0; i < 2*n_aug_+1; i++)
	{
		z_pred = z_pred + weights_(i)*Zsig.col(i);
	}

	for (int i = 0; i < 2*n_aug_+1; i++)
	{
		VectorXd zdiff = Zsig.col(i) - z_pred;

		if (isRadar)
		{
			if (zdiff(1) > M_PI) zdiff(1) -= 2*M_PI;
			else if (zdiff(1) < -M_PI) zdiff(1) += 2*M_PI;
		}

		S = S + weights_(i)*zdiff*zdiff.transpose();
	}
	S += R;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
	int nx = 2;
	VectorXd z_meas = meas_package.raw_measurements_;
	MatrixXd R = MatrixXd(nx,nx);
	R << std_laspx_*std_laspx_, 0,
		 0, std_laspy_*std_laspy_;

	//transform sigma points to measurement space
	MatrixXd Zsig = MatrixXd(nx, 2*n_aug_+1);
	for (int i = 0; i < 2*n_aug_+1; i++)
	{
		Zsig.col(i) = Xsig_pred_.col(i).head(nx);
	}

	//compute mean and covariance of measurement
	VectorXd z_pred = VectorXd(nx);
	MatrixXd S = MatrixXd(nx, nx);
	measMeanAndCovariance(Zsig, R, z_pred, S, false/*isRadar*/);

	MatrixXd Tc = MatrixXd(n_x_, nx);
	Tc.fill(0.0);
	for (int i = 0; i < 2*n_aug_+1; i++)
	{
		VectorXd xdiff = Xsig_pred_.col(i) - x_;
		VectorXd zdiff = Zsig.col(i) - z_pred;
		Tc += weights_(i)*xdiff*zdiff.transpose();
	}

	VectorXd z_diff = z_meas - z_pred;
	MatrixXd K = Tc*S.inverse();
	x_ = x_ + K*z_diff;
	P_ = P_ - K*S*K.transpose();

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
	//compute predicted measurement mean z and covariance S
	VectorXd z_meas = meas_package.raw_measurements_;
	int n_z = 3; //radar measurement
	MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);
	MatrixXd R = MatrixXd(n_z,n_z);
	R << std_radr_*std_radr_, 0, 0,
		 0, std_radphi_*std_radphi_, 0,
		 0, 0, std_radrd_*std_radrd_;

	for (int i = 0; i < 2*n_aug_+1; i++)
	{
		double px = Xsig_pred_(0,i);
		double py = Xsig_pred_(1,i);
		double v = Xsig_pred_(2,i);
		double yaw = Xsig_pred_(3,i);

		double vx = v*cos(yaw);
		double vy = v*sin(yaw);

		Zsig(0,i) = sqrt(px*px + py*py);
		Zsig(1,i) = atan2(py,px);
		if (fabs(Zsig(0,i)) > 0.001)
		{
			Zsig(2,i) = (px*vx + py*vy)/Zsig(0,i);
		}
		else
			Zsig(2,i) = 0;
	}

	//compute mean and covariance of measurement prediction
	VectorXd z_pred = VectorXd(n_z);
	MatrixXd S = MatrixXd(n_z, n_z);
	measMeanAndCovariance(Zsig, R, z_pred, S, true/*isRadar*/);

	MatrixXd Tc = MatrixXd(n_x_, n_z);
	Tc.fill(0.0);
	for (int i = 0; i < 2*n_aug_+1; i++)
	{
		VectorXd xdiff = Xsig_pred_.col(i) - x_;
		VectorXd zdiff = Zsig.col(i) - z_pred;
		Tc += weights_(i)*xdiff*zdiff.transpose();
	}
	MatrixXd K = Tc*S.inverse();

	//update measurement mean and covariance
	VectorXd z_diff = z_meas - z_pred;
	x_ = x_ + K * z_diff;
	P_ = P_ - K * S * K.transpose();

}
