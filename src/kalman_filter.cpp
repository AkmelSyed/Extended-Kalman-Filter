#include "kalman_filter.h"
using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in, MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::SameUpdateStep(VectorXd y) {
  MatrixXd Ht = H_.transpose();
  // Calculating PHt before the calculation of S allows for less matrix multiplications
  // Computes 2 matrix multiplications instead of 3
  MatrixXd PHt = P_ * Ht;
	MatrixXd S = H_ * PHt + R_;
	MatrixXd Si = S.inverse();
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateLidar(const VectorXd &z) {
  VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	SameUpdateStep(y);
}

void KalmanFilter::UpdateRadar(const VectorXd &z) {
  double px = x_(0);
  double py = x_(1);
  double vx = x_(2);
  double vy = x_(3);

  // avoid initializing with nonsense data
  if ( fabs(px) < 0.001 && fabs(py) < 0.001 ) {
   px = 0.1;
   py = 0.1;
  }

  double rho = sqrt(px*px + py*py);
  double phi = atan2(py, px);
  double rho_dot;

  if (fabs(rho) < 0.001) { rho = 0.001; }

  rho_dot = (px*vx + py*vy)/rho;

  VectorXd z_pred = VectorXd(3);
  z_pred << rho, phi, rho_dot;

  VectorXd y = z - z_pred;

  bool positive = false;
  if(y(1)>0) { positive = true; }

  // make sure the phi in the y-vector is between -PI and PI
  y(1) = fabs(y(1));
  while(y(1)>M_PI) { y(1) -= 2*M_PI; }
  if(positive == false) { y(1) = -y(1); }

  SameUpdateStep(y);
}
