import numpy as np
import serial
import time
import struct


class SerialCom:
    """
        Class to interact with the ESP32 microcontroller
    """
    def __init__(self, port, baud):
        self.ser = serial.Serial(port, baud)  # open serial port 115200/256000
        time.sleep(3)
        self.u_prev = 0
        self.t = time.perf_counter()
        self.u_prev_dec = 0
        self.x = np.array([0, 0, 0, 0]).reshape(4, 1)

    def apply(self, u):
        """
        Function to apply the current input to the pendulum
        :param: (float) u - Input for the current time step
        """

        self.ser.write(struct.pack('f', u))
        self.ser.write(b"\n")

    def measure(self):
        """
        Function to read the states of the pendulum sent by the ESP32 microcontroller
        :return: x - (float- 4 x 1) - Pendulum states at current time
        """
        buffer = self.ser.read_all()
        buffer = buffer[-42:]
        print("current buffer len: ", len(buffer))
        if len(buffer) <= 21:
            x = self.x
        else:
            for i in range(21, -1, -1):
                byte = buffer[i:i + 1]
                if byte == b"?":
                    print("message buffer len: ", len(buffer[i + 1:i + 21]))

                    # Input
                    u_rec = buffer[i + 1:i + 5]
                    if u_rec != self.u_prev:
                        print("time taken to update u ", time.perf_counter() - self.t)
                        self.t = time.perf_counter()
                    self.u_prev = u_rec
                    u_rec = struct.unpack_from('f', u_rec, offset=0)[0]
                    self.u_prev_dec = u_rec
                    print(f'u is {u_rec}')

                    # Angular Position
                    data1 = buffer[i + 5:i + 9]
                    data1 = struct.unpack_from('f', data1, offset=0)[0]

                    # Angular Velocity
                    data2 = buffer[i + 9:i + 13]
                    data2 = struct.unpack_from('f', data2, offset=0)[0]

                    # Linear Position
                    data3 = buffer[i + 13:i + 17]
                    data3 = struct.unpack_from('f', data3, offset=0)[0]

                    # Linear Velocity
                    data4 = buffer[i + 17:i + 21]
                    data4 = struct.unpack_from('f', data4, offset=0)[0]

                    x = np.array([data1, data2, data3, data4]).reshape(4, 1)
                    self.x = x
                    break
            if 'x' in locals():
                pass
            else:
                print(buffer)
        return x