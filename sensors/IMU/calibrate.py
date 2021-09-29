import serial
import numpy as np
import time 

def readData(ser, duration):
    start = time.time()
    container = []
    
    while time.time() - start < duration:
        try:
            ser_bytes = ser.readline()
            data = ser_bytes.decode("utf-8").strip().split(' ')
            data = np.asarray(data, dtype=np.float32)
            container.append(data)
            
        except Exception as e:
            print(e)
            break
        
    return container

def calibrateGyro():
    print(f'Put IMU into a static position and press Enter')
    input('')
    duration = 7
    data = readData(ser, duration)
    means = np.mean(data, axis=0)
    return means

def calibrateMagnetometer():
    lam = 0.99
    RLS_theta = np.zeros((4,1))
    RLS_P = np.eye(4) * 1000
    RLS_in = np.zeros((4,1))
    RLS_out = 0
    RLS_gain = np.zeros((4,1))

    R_MAG = np.dot(RLS_P.transpose(), RLS_P)[0,0]

    print("Pick up the IMU and start rotating it in every direction")
    time.sleep(3)

    i = 0
    
    while R_MAG > 1e-6:
        try:
            ser_bytes = ser.readline()
            data = ser_bytes.decode("utf-8").strip().split(' ')
            data = np.asarray(data, dtype=np.float32)
        except Exception as e:
            print(e)
            break

        if i%10 == 0:
            RLS_in[0,0] = data[0]
            RLS_in[1,0] = data[1]
            RLS_in[2,0] = data[2]
            RLS_in[3,0] = 1
            RLS_out = np.linalg.norm(data) ** 2

            err = RLS_out - RLS_in.transpose() @ RLS_theta
            RLS_gain = RLS_P @ RLS_in / (lam + RLS_in.transpose() @ RLS_P @ RLS_in)
            RLS_P = (RLS_P - RLS_gain @ RLS_in.transpose() @ RLS_P) / lam

            RLS_theta = RLS_theta + err * RLS_gain

            R_MAG = np.dot(RLS_P.transpose(), RLS_P)[0,0]
            print(f"Covariange norm: {R_MAG}")
        i += 1 
  

    HARD_IRON_BIAS = RLS_theta/2.0
    return HARD_IRON_BIAS[:3]

def gatherAccelData():
    container = np.zeros((9,3))
    print("To calibrate the accelerometer, data from 9 static positions is needed")
    print("Put the IMU in one of the 9 orientations and press Enter")
    input()
    for i in range(9):
        try:
            ser.flushInput()
            ser_bytes = ser.readline()
            data = ser_bytes.decode("utf-8").strip().split(' ')
            data = np.asarray(data, dtype=np.float32)
            container[i,:] = data
            if i != 8:
                print("Switch orientation and press Enter again")
                input()
            else:
                print("Measurement procedure finished")
        except Exception as e:
            print(e)
            break

    np.savetxt('accelData.txt', container)


def recordVariance():
    print("Put the IMU into a stationary position and press Enter")
    input()
    dur = 10
    data = readData(ser, dur)
    var = np.var(data, axis=0)
    return var
    

if __name__ == "__main__":
    ser = serial.Serial('COM3')
    ser.flushInput()
    #gyroBias = calibrateGyro()
    #magBias = calibrateMagnetometer()
    #gatherAccelData()
    var = recordVariance()
    
