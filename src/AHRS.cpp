#include "AHRS.h"

#include <math.h>

#include <Eigen/LU>

void AHRS::calcDCM(const Quaternion &X, RotationMatrix &DCM) {
    float_prec q0, q1, q2 ,q3;
    q0 = X(0);
    q1 = X(1);
    q2 = X(2);
    q3 = X(3);

    DCM(0,0) = q0*q0 + q1*q1 - q2*q2 - q3*q3;
    DCM(0,1) = 2 * (q1*q2 + q0*q3);
    DCM(0,2) = 2 * (q1*q3 - q0*q2);
    DCM(1,0) = 2 * (q1*q2 - q0*q3);
    DCM(1,1) = q0*q0 - q1*q1 + q2*q2 - q3*q3;
    DCM(1,2) = 2 * (q2*q3 + q0*q1);
    DCM(2,0) = 2 * (q1*q3 + q0*q2);
    DCM(2,1) = 2 * (q2*q3 - q0*q1);
    DCM(2,2) = q0*q0 - q1*q1 - q2*q2 + q3*q3;
}

void AHRS::calcDCM(const Euler& angles, RotationMatrix& DCM) {
    float_prec yaw, pitch, roll;
    roll = angles(0);
    pitch = angles(1);
    yaw = angles(2);

    DCM(0,0) = cos(yaw) * cos(pitch);
    DCM(0,1) = sin(yaw) * cos(pitch);
    DCM(0,2) = -sin(pitch);
    DCM(1,0) = sin(roll) * sin(pitch) * cos(yaw) - sin(yaw) * cos(roll);
    DCM(1,1) = sin(roll) * sin(yaw) * sin(pitch) + cos(roll) * cos(yaw);
    DCM(1,2) = sin(roll) * cos(pitch);
    DCM(2,0) = sin(roll) * sin(yaw) + sin(pitch) * cos(roll) * cos(yaw);
    DCM(2,1) = -sin(roll) * cos(yaw) + sin(yaw) * sin(pitch) * cos(roll);
    DCM(2,2) = cos(roll) * cos(pitch);
}

float_prec AHRS::getYaw(const Vec3D &A, const Vec3D &B) {
    float_prec roll = getRoll(A);
    float_prec pitch = getPitch(A);

    Euler angles;
    angles(0) = roll;
    angles(1) = pitch;
    angles(0) = 0.0;

    RotationMatrix DCM;
    calcDCM(angles, DCM);

    Vec3D BW = DCM.transpose() * B;

    return atan2(BW(1), BW(0));

}

void AHRS::quat2euler(const Quaternion& quat, Euler& angles) {
    float_prec q0, q1, q2, q3;
    q0 = quat(0);
    q1 = quat(1);
    q2 = quat(2);
    q3 = quat(3);

    angles(0) = atan2(2 * (q0*q1 + q2*q3), 1 - 2 * (q1*q1 + q2*q2));
    angles(1) = asin(2 * (q0*q2 - q3*q1));
    angles(2) = atan2(2 * (q0*q3 + q1*q2), 1 - 2 * (q2*q2 + q3*q3));
}

void AHRS::euler2quat(const Euler& angles, Quaternion& quat) {
    float_prec roll = angles(0);
    float_prec pitch = angles(1);
    float_prec yaw = angles(2);

    float_prec cy = cos(yaw * 0.5);
    float_prec sy = sin(yaw * 0.5);
    float_prec cp = cos(pitch * 0.5);
    float_prec sp = sin(pitch * 0.5);
    float_prec cr = cos(roll * 0.5);
    float_prec sr = sin(roll * 0.5);

    quat(0) = cr * cp * cy + sr * sp * sy;
    quat(1) = sr * cp * cy - cr * sp * sy;
    quat(2) = cr * sp * cy + sr * cp * sy;
    quat(3) = cr * cp * sy - sr * sp * cy;
}

void AHRS::Estimator::DEKF::init(const Quaternion& X0, float_prec dt) {
    this->X = X0;
    this->_dt = dt;
}

void AHRS::Estimator::DEKF::update(const inputVector& U, const outputVector& Y) {
    // * Predict the next state using the non-linear state function
    // 1. X[k|k-1] = F(X[k-1], U[k])
    UpdateX(X, U, Xp);

    // Update the Covariance matrix using state Jacobian
    // 2. dF = dF(X[t-1], U[t]) / dx
    JacobianF(X, U, dF);

    // 3. P[k|k-1] = dF @ P @ dF.T + Q
    P = (dF * (P * dF.transpose())) + Q;

    // Compute Jacobian of measurement function in order to correct the prediction
    // 4. dH = dH(X[t-1], U[t]) / dx
    JacobianH(X, U, dH);

    // Compute the Kalman gain
    // 5. S = dH @ P @ dH.T + R
    S = (dH * (P * dH.transpose())) + R;
    // 6. K = P @ dH.T @ S^(-1)
    K = P * (dH.transpose() * S.inverse());

    // Update measured state using the non-linear measurement function
    // 7. y[k] = H([X[k-1], U[k])
    UpdateY(X, U, Yp);

    // Correct the state prediction using the measurement and update the Covariance matrix
    // 8. X[k] = X[k|k-1] + (K @ (y[k] - y[k|k-1]))
    X = Xp +  (K * (Y - Yp));
    // 9. P[k] = (I - (K @ dH)) @ P
    P = (I - (K * dH)) * P;
}

void AHRS::Estimator::DEKF::UpdateX(const stateVector& Xprev, const inputVector& U, stateVector& Xpred) {
    float_prec p, q, r;
    p = U(0);
    q = U(1);
    r = U(2);

    float_prec q0, q1, q2, q3;
    q0 = Xprev(0);
    q1 = Xprev(1);
    q2 = Xprev(2);
    q3 = Xprev(3);

    Xpred(0) += (_dt/2) * (-p*q1 - q*q2 - r*q3);
    Xpred(1) += (_dt/2) * (p*q0 + r*q2 - q*q3);
    Xpred(2) += (_dt/2) * (q*q0 - r*q1 + p*q3);
    Xpred(3) += (_dt/2) * (r*q0 + q*q1 - p*q2);

    Xpred.normalize();
}

void AHRS::Estimator::DEKF::UpdateY(const stateVector& Xprev, const inputVector& U, outputVector& Yp) {
    float_prec q0, q1, q2, q3;
    q0 = Xprev(0);
    q1 = Xprev(1);
    q2 = Xprev(2);
    q3 = Xprev(3);

    float_prec q02 = q0 * q0;
    float_prec q12 = q1 * q1;
    float_prec q22 = q2 * q2;
    float_prec q32 = q3 * q3;

    Yp(0) = 2 * (q1*q3 - q0*q2);
    Yp(1) = 2 * (q2*q3 - q0*q1);
    Yp(2) = q02 - q12 - q22 + q32;
    Yp(3) = B0(0) * (q02 + q12 - q22 - q32) + B0(1) * (2 * (q1*q2 + q0*q3)) + B0(2) * (2 * (q1*q3 - q0*q2));
    Yp(4) = B0(0) * (2 * (q1*q2 - q0*q3)) + B0(1) * (q02 - q12 + q22 - q32) + B0(2) * (2 * (q2*q3 + q0*q1));
    Yp(5) = B0(0) * (2 * (q1*q3 + q0*q2)) + B0(1) * (2 * (q2*q3 - q0*q1)) + B0(2) * (q02 - q12 - q22 + q32);
}

void AHRS::Estimator::DEKF::JacobianF(const stateVector& Xprev, const inputVector& U, stateMatrix& dF) {
    float_prec p, q, r;
    p = U(0);
    q = U(1);
    r = U(2);

    dF(0,0) = 1.0;
    dF(0,1) = -0.5 * _dt * p;
    dF(0,2) = -0.5 * _dt * q;
    dF(0,3) = -0.5 * _dt * r;
    dF(1,0) = 0.5 * _dt * p;
    dF(1,1) = 1.0;
    dF(1,2) = 0.5 * _dt * r;
    dF(1,3) = -0.5 * _dt * q;
    dF(2,0) = 0.5 * _dt * q;
    dF(2,1) = -0.5 * _dt * r;
    dF(2,2) = 1.0;
    dF(2,3) = 0.5 * _dt * p;
    dF(3,0) = 0.5 * _dt * r;
    dF(3,1) = 0.5 * _dt * q;
    dF(3,2) = -0.5 * _dt * p;
    dF(3,3) = 1.0;
}

void AHRS::Estimator::DEKF::JacobianH(const stateVector& Xprev, const inputVector& U, outputMatrix& dH) {
    float_prec q0, q1, q2, q3;

    q0 = Xprev(0);
    q1 = Xprev(1);
    q2 = Xprev(2);
    q3 = Xprev(3);

    dH(0,0) = -2*q2;
    dH(0,1) =  2*q3;
    dH(0,2) = -2*q0;
    dH(0,3) =  2*q1;

    dH(1,0) = 2*q1;
    dH(1,1) = 2*q0;
    dH(1,2) = 2*q3;
    dH(1,3) = 2*q2;

    dH(2,0) =  2*q0;
    dH(2,1) = -2*q1;
    dH(2,2) = -2*q2;
    dH(2,3) =  2*q3;

    dH(3,0) =  2*q0*B0(0) + 2*q3*B0(1) - 2*q2*B0(2);
    dH(3,1) =  2*q1*B0(0) + 2*q2*B0(1) + 2*q3*B0(2);
    dH(3,2) = -2*q2*B0(0) + 2*q1*B0(1) - 2*q0*B0(2);
    dH(3,3) = -2*q3*B0(0) + 2*q0*B0(1) + 2*q1*B0(2);

    dH(4,0) = -2*q3*B0(0) + 2*q0*B0(1) + 2*q1*B0(2);
    dH(4,1) =  2*q2*B0(0) - 2*q1*B0(1) + 2*q0*B0(2);
    dH(4,2) =  2*q1*B0(0) + 2*q2*B0(1) + 2*q3*B0(2);
    dH(4,3) = -2*q0*B0(0) - 2*q3*B0(1) + 2*q2*B0(2);

    dH(5,0) =  2*q2*B0(0) - 2*q1*B0(1) + 2*q0*B0(2);
    dH(5,1) =  2*q3*B0(0) - 2*q0*B0(1) - 2*q1*B0(2);
    dH(5,2) =  2*q0*B0(0) + 2*q3*B0(1) - 2*q2*B0(2);
    dH(5,3) =  2*q1*B0(0) + 2*q2*B0(1) + 2*q3*B0(2);
}
