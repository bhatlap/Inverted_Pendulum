import serial
import time
import struct
import numpy as np

class SerialCom:
    def __init__(self, port, baud):
        self.ser = serial.Serial(port, baud)  # open serial port 115200/250000
        self.t=0

    def send(self, u, check):
        self.t = time.time()
        print('')
        inpuT = struct.pack('f', u)
        print('input:')
        print(u)
        self.ser.write(inpuT)
        self.ser.write(b"\n")

    def receive(self):
        data1 = self.ser.read(4)
        data1 = struct.unpack_from('f', data1, offset=0)[0]
        #print('theta:')
        #print(data - check)
        # if abs(data-check)>=0.5:
        #     print('error')
        #     break

        data2 = self.ser.read(4)
        data2 = struct.unpack_from('f', data2, offset=0)[0]
        #print('x:')
        #print(data1)
        
        data3 = self.ser.read(4)
        data3 = struct.unpack_from('f', data3, offset=0)[0]
        
        data4 = self.ser.read(4)
        data4 = struct.unpack_from('f', data4, offset=0)[0]

        elapsed = time.time() - self.t
        time.sleep(0.07 - elapsed)
        print('freq:')
        print(1 / (time.time() - self.t))
        x = np.array([data1,data2,data3,data4])
        return x



with open('initialguessU.txt') as f:
    for line in f:
        currentline = line.split(",")
input = [float(i) for i in currentline]

with open('initialguessX.txt') as f:
    for line in f:
        currentline = line.split(",")
inputX = [float(i) for i in currentline]

ser = SerialCom('COM7', 250000)
for i in range(len(input)):
    ser.send(input[i], inputX[i])
    x9 = ser.receive()
    print(x9)
    
