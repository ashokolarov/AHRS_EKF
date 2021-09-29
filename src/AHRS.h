#ifndef AHRS_H
#define AHRS_H

#include <Eigen.h>

#include "config.h"

namespace AHRS{
    using namespace Eigen;

    typedef Matrix<float_prec, 4, 1> Quaternion;
    typedef Matrix<float_prec, 1, 3> Euler;
    typedef Matrix<float_prec, 3, 3> RotationMatrix;
    typedef Matrix<float_prec, 3, 1> Vec3D;

    /**********************************************************
     * Calculate pitch based on set of accelerometer readings.
     * @param A -> Vec3D containing X,Y and Z accelerometer readings
     * @return  -> float_prec computed pitch in [rad].
    **********************************************************/
    inline float_prec getPitch(const Vec3D& A){
        return atan2(-A(0), A(2));
    }

    /**********************************************************
     * Calculate roll based on set of accelerometer readings.
     * @param A -> Vec3D containing X,Y and Z accelerometer readings
     * @return  -> float_prec computed roll in [rad].
    **********************************************************/
    inline float_prec getRoll(const Vec3D& A){
        return atan2(A(1), sqrt(A(0)*A(0) + A(2)*A(2)));
    }

    /***********************************************************
     * Compute yaw by transforming magnetometer measurements
     * using a DCM computed using pitch and roll from
     * accelerometer
     * @param A -> Vec3D containing X,Y and Z accelerometer readings
     * @param B -> Vec3D containing X,Y and Z magnetometer readings
     * @return  -> float_prec yaw in [rad]
     */
    float_prec getYaw(const Vec3D& A, const Vec3D& B);

    /**********************************************************
     * Compute Direction Cosine Matrix to transform from World
     * to Body reference frame.
     * @param X -> 4x1 Eigen matrix holding a quaternion
     * @param DCM -> 3x3 Eigen matrix to output the DCM
    **********************************************************/
    void calcDCM(const Quaternion& X, RotationMatrix& DCM);

    /**********************************************************
     * Compute Direction Cosine Matrix to transform from World
     * to Body reference frame.
     * @param angles -> 3x1 Eigen matrix holding a Euler angles
     * @param DCM -> 3x3 Eigen matrix to output the DCM
    **********************************************************/
    void calcDCM(const Euler& angles, RotationMatrix& DCM);

    /*************************************************************
     * Convert a set of Euler angles [roll, pitch, yaw] into a
     * unit quaternion
     * @param angles -> Euler vector holding the roll, pitch and yaw
     * angles in [rad]
     * @param quat -> Empty quaternion to write result to
    **************************************************************/
    void euler2quat(const Euler& angles, Quaternion& quat);

    /***********************************************************************
     * Convert a unit quaternion to a set of Euler angles [roll, pitch, yaw]
     * @param quat -> Quaternion holding the orientation
     * @param angles -> Euler vector holding the roll, pitch and yaw
     * angles in [rad]
    ***********************************************************************/
    void quat2euler(const Quaternion& quat, Euler& angles);

    namespace Estimator{
        static constexpr uint8_t dimX = 4;
        static constexpr uint8_t dimU = 3;
        static constexpr uint8_t dimY = 6;

        typedef Matrix<float_prec, dimX, 1> stateVector;
        typedef Matrix<float_prec, dimX, dimX> stateMatrix;
        typedef Matrix<float_prec, dimU, 1> inputVector;
        typedef Matrix<float_prec, dimY, 1> outputVector;
        typedef Matrix<float_prec, dimY, dimX> outputMatrix;

        static const Quaternion X0 = {1.0, 0.0, 0.0, 0.0}; // Default initial state estimate
        static Vec3D B0 = {-1.0, 0.0, 0.0};

        // Discrete Extended Kalman Filter class used to estimate attitude in the form of a Qaternion using measurements
        // from a 9-DOF IMU.
        class DEKF{

        public:

            /*****************************************************
             * Initialize the EKF with an initial state guess and
             * system time step
             * @param X0 -> Quaternion representing initial state
             * @param dt -> timestep used in [s]
            *****************************************************/
            void init(const Quaternion& X0, float_prec dt);

            /*************************************************************************************
            * Compute the state of the system at current iteration
            * @param Xprev -> state of the system at previous iteration
            * @param U -> Input vector holding gyroscope measurements
            * @param Y -> Measurement vector holding accelerometer and magnetometer measurements.
            *************************************************************************************/
            void update(const inputVector& U, const outputVector& Y);

            /***************************************************
             * Obtain current state estimate
             * @return -> quaternion representing current state
            ***************************************************/
            inline Quaternion getX() { return X; }

            /***************************************************
             * Obtain current measurement prediction
             * @return -> output vector with accelerometer
             * and magnetometer predicted readings
            ****************************************************/
            inline outputVector getY() { return Yp; }

        private:
            constexpr static float_prec P_init = 100.0;
            constexpr static float_prec Q_init = 1E-5;
            constexpr static float_prec R_init = 0.5;

            stateMatrix P = stateMatrix::Identity(dimX, dimX) * P_init; // Covariance Matrix
            stateMatrix Q = stateMatrix::Identity(dimX, dimX) * Q_init; // Process noise
            Matrix<float_prec, dimY, dimY> R = MatrixXf::Identity(dimY, dimY) * R_init; // Measurement noise

            float_prec _dt;

            stateVector X;
            stateVector Xp;
            outputVector Yp;
            stateMatrix dF;
            outputMatrix dH;
            Matrix<float_prec, dimY, dimY> S;
            Matrix<float_prec, dimX, dimY> K;
            Matrix<float_prec, dimX, dimX> I = MatrixXf::Identity(dimX, dimX);

            /******************************************************************************
             * Predict state at k+1 using non-linear state update function F(X[k-1], U[k])
             * @param Xprev -> quaternion holding previous state estimate
             * @param U     -> inputVector holding gyroscope readings
             * @param Xpred -> quaternion to write the predicted state to
            *******************************************************************************/
            void UpdateX(const stateVector& Xprev, const inputVector& U, stateVector& Xpred);

            /*********************************************************************************************
             * Predict measurement vector using non-linear function H(X[k-1], U[k])
             * @param Xprev -> quaternion holding previous state estimate
             * @param Yp    -> outputVector to write the predicted accelerometer and magnetometer readings to
             * @param U     -> inputVector holding gyroscope readings
            **********************************************************************************************/
            void UpdateY(const stateVector& Xprev, const inputVector& U, outputVector& Yp);

            /********************************************************************
             * Compute the Jacobian of the non-linear function F w.r.p. to the state
             * @param Xprev -> quaternion holding previous state estimate
             * @param U     -> inputVector holding gyroscope readings
             * @param dF    -> stateMatrix holding the numerical values
            *********************************************************************/
            void JacobianF(const stateVector& Xprev, const inputVector& U, stateMatrix& dF);

            /********************************************************************
             * Compute the Jacobian of the non-linear function H w.r.p. to the state
             * @param Xprev -> quaternion holding previous state estimate
             * @param U     -> inputVector holding gyroscope readings
             * @param dH    -> outputMatrix holding the numerical values
            *********************************************************************/
            void JacobianH(const stateVector& Xprev, const inputVector& U, outputMatrix& dH);
        };
    }
}

#endif
