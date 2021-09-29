#include "config.h"
#include "mpu9250.h"
#include "AHRS.h"
#include "tools.h"

// containers for IMU data
AHRS::Vec3D A;
AHRS::Vec3D W;
AHRS::Vec3D B;

// Attitude estimation state and EFK filter
AHRS::Quaternion X_cur;
AHRS::Vec3D U;
AHRS::Estimator::outputVector Y;
AHRS::Estimator::DEKF attitudeFilter;

// Declare IMU object
SimpleMPU9250 IMU(Wire, 0x68);

unsigned long previousMicros;

void setup() {
  Serial.begin(115200);

  // Initialize IMU
  int IMUstatus = IMU.begin();   /* start communication with IMU */
    if (IMUstatus < 0) {
        Serial.println("IMU initialization unsuccessful");
        Serial.println("Check IMU wiring or try cycling power");
        Serial.print("Status: ");
        Serial.println(IMUstatus);
        while(1) {}
    }

  // set Acceleromter and Gyro range as well as the Digital low-pass filter
  IMU.setAccelRange(IMU.ACCEL_RANGE_2G);
  IMU.setGyroRange(IMU.GYRO_RANGE_500DPS);
  IMU.setDlpfBandwidth(IMU.DLPF_BANDWIDTH_92HZ);

  attitudeFilter.init(AHRS::Estimator::X0, SYSTEM::dt);

  previousMicros = micros();
}

void loop() {
  if ((micros() - previousMicros) > SYSTEM::dtMicros) {

    // Read IMU
    IMU.readSensor();
    
    A(0) = IMU.ax;
    A(1) = IMU.ay;
    A(2) = IMU.az;
    A.normalize();
    
    B(0) = IMU.hx;
    B(1) = IMU.hy;
    B(2) = IMU.hz;
    B.normalize();

    W(0) = IMU.gx;
    W(1) = IMU.gy;
    W(2) = IMU.gz;

    // Update Kalman Filter
    U(0) = W(0);
    U(1) = W(1);
    U(2) = W(2);

    Y(0) = A(0);
    Y(1) = A(1);
    Y(2) = A(2);
    Y(3) = B(0);
    Y(4) = B(1);
    Y(5) = B(2);

    attitudeFilter.update(U, Y);
    X_cur = attitudeFilter.getX();
    print_visual(X_cur);
   
    previousMicros = micros();
  }

  
}
