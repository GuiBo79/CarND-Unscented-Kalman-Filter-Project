#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
    
    
  ///* State dimension
  n_x_ = 5;
    
  ///* Augmented state dimension
  n_aug_ = 7;
    
  ///* Sigma point spreading parameter
  lambda_ = 3 - n_aug_;
    
  // initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;
    
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_.setIdentity();
    
  //predicted sigma points matrix
  Xsig_pred_ = Eigen::MatrixXd (n_x_ , 2*n_aug_ + 1);
    
  // Weights of sigma points
  weights_ = VectorXd(2*n_aug_ + 1);
    
  // time when the state is true, in us
  time_us_ = 0;
    
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.5;

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
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    
    
    
    if (is_initialized_ == false){
        time_us_ = meas_package.timestamp_;
        
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
            
            double ro = meas_package.raw_measurements_[0];
            double theta = meas_package.raw_measurements_[1];
            //double ro_dot = meas_package.raw_measurements_[2];
            
            double px = ro * cos(theta);
            double py = ro * sin(theta);
            x_ << px,py,0,0,0;
            
        }
        
        if (meas_package.sensor_type_ == MeasurementPackage::LASER){
            
            double px = meas_package.raw_measurements_[0];
            double py = meas_package.raw_measurements_[1];
            
            x_ << px,py,0,0,0;
            
            
        }
        
        is_initialized_ = true;
            
            
    }
    
    double  delta_t = (meas_package.timestamp_ - time_us_)/1000000.0;
    
    Prediction(delta_t);
    time_us_ = meas_package.timestamp_;
    
    
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR  && use_radar_ == true){
        
        UpdateRadar (meas_package);
        
    }
    
    if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_ == true){
        
        UpdateLidar (meas_package);
        
    }

    
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
    
    /******Prediction of the Sigma Points **********//////
    
    //create augmented mean vector
    VectorXd x_aug = VectorXd(7);
    
    //create augmented state covariance
    MatrixXd P_aug = MatrixXd(7, 7);
    
    //create sigma point matrix
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
    
    //create matrix with predicted sigma points as columns
    
    MatrixXd State_Aux = MatrixXd(5,1);
    
    x_aug << x_,0,0;
    
    // Q Matrix
    MatrixXd Q(2,2);
    Q << std_a_*std_a_,0,
    0, std_yawdd_*std_yawdd_;
    
    //P augmented Matrix initialization
    P_aug << P_,
    MatrixXd::Zero(n_x_, n_aug_ - n_x_),
    MatrixXd::Zero(n_aug_ - n_x_, n_x_),
    Q;
    
    
    
    //Square root of P_aug Matrix
    MatrixXd A = P_aug.llt().matrixL();
    
    //set first column of sigma point matrix
    Xsig_aug.col(0)  = x_aug;
    
    //set remaining sigma points
    for (int i = 0; i < n_aug_; i++)
    {
        Xsig_aug.col(i+1)     = x_aug + sqrt(lambda_+n_aug_) * A.col(i);
        Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * A.col(i);
    }
    
    
    for (int i=0 ; i < 2*n_aug_ +1 ;i++){
        
        //predict sigma points
        double px = Xsig_aug(0,i);
        double py = Xsig_aug(1,i);
        double v = Xsig_aug(2,i);
        double yaw = Xsig_aug(3,i);
        double yaw_dot = Xsig_aug(4,i);
        double nu_a = Xsig_aug(5,i);
        double nu_yaw_ddot = Xsig_aug(6,i);
        double p_px = 0;
        double p_py = 0;
        
        //Auxiliar State Vector
        State_Aux << 0.5*delta_t*delta_t*cos(yaw)*nu_a,
        0.5*delta_t*delta_t*sin(yaw)*nu_a,
        delta_t*nu_a,
        0.5*delta_t*delta_t*nu_yaw_ddot,
        delta_t*nu_yaw_ddot;
        
        
        
        
        //PX and PY Prediction Calculation Calculation
        if (fabs(yaw_dot) > 0.0001){
            
            p_px = px + v/yaw_dot*(sin(yaw+yaw_dot*delta_t)-sin(yaw));
            p_py = py + v/yaw_dot*(-cos(yaw+yaw_dot*delta_t)+cos(yaw));
            
        } else {
            
            p_px = px + v*cos(yaw)*delta_t;
            p_py = py + v*sin(yaw)*delta_t;
            
        }
        
        //Sigma Points Prediction Matrix
        Xsig_pred_.col(i) << p_px + State_Aux(0),
                             p_py + State_Aux(1),
                             v + State_Aux(2),
                             yaw + yaw_dot*delta_t + State_Aux(3),
                             yaw_dot + State_Aux(4);
        
    }
    
    
    
    /////****Prediction of Mean and covariance Matrix ****///////
    
    //set weights
    weights_(0) = lambda_/(lambda_+n_aug_);
    
    for (int i = 1; i < 2*n_aug_ + 1; i++){
        
        weights_(i) = 1/(2*(lambda_+n_aug_));
        
    }
    
    //predict state mean
    x_.fill(0.0);
    for (int i = 0; i < 2*n_aug_ + 1; i++){
        
        x_ = x_ + weights_(i)*Xsig_pred_.col(i);
        
    }
    
    //predict state covariance matrix
    P_.fill(0.0);
    for (int i = 0; i < 2*n_aug_ +1; i++){
        
        VectorXd diff = Xsig_pred_.col(i) - x_;
        //angle normalization
        while (diff(3)> M_PI) diff(3)-=2.*M_PI;
        while (diff(3)<-M_PI) diff(3)+=2.*M_PI;
        
        P_ = P_ + weights_(i)*diff*diff.transpose();
    }
    
    
    
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
    
    int n_z_ = 2;
    
    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);
    
    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z_);
    
    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z_,n_z_);
    
    //transform sigma points into measurement space
    Zsig = Xsig_pred_.block(0, 0, n_z_, 2*n_aug_+1);
    
    //calculate mean predicted measurement
    z_pred.fill(0.0);
    for (int i = 0; i < 2*n_aug_ + 1; i++){
        
        z_pred = z_pred + weights_(i)*Zsig.col(i);
        
    }
    
    //calculate measurement covariance matrix S
    S.fill(0.0);
    MatrixXd R = MatrixXd(2,2);
    R << std_laspx_*std_laspx_, 0,
         0, std_laspy_*std_laspy_;
    
    
    
    
    for (int i = 0; i < 2*n_aug_ + 1; i++){
        
        VectorXd Z_diff = Zsig.col(i) - z_pred;
        
        S = S + weights_(i)*Z_diff*Z_diff.transpose();
    }
    
    S = S + R;
    
    
    
    //Z vector from raw measurements
    VectorXd z = VectorXd(n_z_);
    z << meas_package.raw_measurements_[0],
         meas_package.raw_measurements_[1];
    
    
    
    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z_);
    MatrixXd K = MatrixXd(n_x_, n_z_);
    
    
    
    //calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < 2*n_aug_ + 1; i++){
        
        VectorXd X_diff = VectorXd(n_x_);
        VectorXd Z_diff = VectorXd(n_z_);
        
        X_diff = Xsig_pred_.col(i) - x_;
        //angle normalization
        while (X_diff(3)> M_PI) X_diff(3)-=2.*M_PI;
        while (X_diff(3)<-M_PI) X_diff(3)+=2.*M_PI;
        
        Z_diff = Zsig.col(i) - z_pred;
        //angle normalization
        while (Z_diff(1)> M_PI) Z_diff(1)-=2.*M_PI;
        while (Z_diff(1)<-M_PI) Z_diff(1)+=2.*M_PI;
        
        Tc = Tc + weights_(i)*X_diff*Z_diff.transpose();
        
    }
    
    
    //calculate Kalman gain K;
    K.fill(0.0);
    K = Tc*S.inverse();
    
    
    
    //update state mean and covariance matrix
    VectorXd Z_diff = z - z_pred;
    while (Z_diff(1)> M_PI) Z_diff(1)-=2.*M_PI;
    while (Z_diff(1)<-M_PI) Z_diff(1)+=2.*M_PI;
    
    
    
    x_ = x_ + K*(Z_diff);
    P_ = P_ - K*S*K.transpose();
    
    NIS_laser_ = Z_diff.transpose()*S.inverse()*Z_diff;
    cout << "NIS laser:" << NIS_laser_ << endl;
    
    
  /*
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
    
    int n_z_ = 3;
    
    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);
    
    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z_);
    
    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z_,n_z_);
    
    //transform sigma points into measurement space
    for (int i = 0; i < 2*n_aug_ + 1; i++){
        
        double px = Xsig_pred_.col(i)(0);
        double py = Xsig_pred_.col(i)(1);
        double v = Xsig_pred_.col(i)(2);
        double yaw = Xsig_pred_.col(i)(3);
        //double yaw_dot = Xsig_pred_.col(i)(4);
        
        double rho = sqrt(px*px+py*py);
        double phi = atan2(py,px);
        double rho_dot = (px*cos(yaw)*v + py*sin(yaw)*v)/rho;
        
        Zsig.col(i) << rho,
        phi,
        rho_dot;
        
        
    }
    
    //calculate mean predicted measurement
    z_pred.fill(0.0);
    for (int i = 0; i < 2*n_aug_ + 1; i++){
        
        z_pred = z_pred + weights_(i)*Zsig.col(i);
        
    }
    
    //calculate measurement covariance matrix S
    S.fill(0.0);
    MatrixXd R = MatrixXd(3,3);
    R << std_radr_*std_radr_, 0, 0,
    0, std_radphi_*std_radphi_, 0,
    0, 0, std_radrd_*std_radrd_;
    
    for (int i = 0; i < 2*n_aug_ + 1; i++){
        
        VectorXd Z_diff = Zsig.col(i) - z_pred;
        
        S = S + weights_(i)*Z_diff*Z_diff.transpose();
    }
    
    S = S + R;
    
    
    
    //Z vector from raw measurements
    VectorXd z = VectorXd(n_z_);
    z << meas_package.raw_measurements_[0],
         meas_package.raw_measurements_[1],
         meas_package.raw_measurements_[2];
    
    
    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z_);
    MatrixXd K = MatrixXd(n_x_, n_z_);
    
    
    
    //calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < 2*n_aug_ + 1; i++){
        
        VectorXd X_diff = VectorXd(n_x_);
        VectorXd Z_diff = VectorXd(n_z_);
        
        X_diff = Xsig_pred_.col(i) - x_;
        //angle normalization
        while (X_diff(3)> M_PI) X_diff(3)-=2.*M_PI;
        while (X_diff(3)<-M_PI) X_diff(3)+=2.*M_PI;
        
        Z_diff = Zsig.col(i) - z_pred;
        //angle normalization
        while (Z_diff(1)> M_PI) Z_diff(1)-=2.*M_PI;
        while (Z_diff(1)<-M_PI) Z_diff(1)+=2.*M_PI;
        
        Tc = Tc + weights_(i)*X_diff*Z_diff.transpose();
        
    }
    
    
    //calculate Kalman gain K;
    K.fill(0.0);
    K = Tc*S.inverse();
    
    
    
    //update state mean and covariance matrix
    VectorXd Z_diff = z - z_pred;
    while (Z_diff(1)> M_PI) Z_diff(1)-=2.*M_PI;
    while (Z_diff(1)<-M_PI) Z_diff(1)+=2.*M_PI;
    
    
    
    x_ = x_ + K*(Z_diff);
    P_ = P_ - K*S*K.transpose();
    
    NIS_radar_ = Z_diff.transpose()*S.inverse()*Z_diff;
    cout << "NIS radar:" << NIS_radar_ << endl;
    
    
    
  /*
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  
}
