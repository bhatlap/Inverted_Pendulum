import serial
import time
import struct
import numpy as np
import matplotlib.pyplot as plt


class SerialCom:
    def __init__(self, port, baud):
        self.ser = serial.Serial(port, baud)  # open serial port 115200/256000
        self.t = 0

    def send(self, u):
        self.t = time.time()
        # print('')
        inpuT = struct.pack('f', u)
        # print('input:')
        # print(inpuT.hex())
        self.ser.write(inpuT)
        self.ser.write(b"\n")

    def receive(self):
        data1 = self.ser.read(4)
        # print(data1.hex())
        data1 = struct.unpack_from('f', data1, offset=0)[0]
        # print('u receive:')
        # print(data1)

        data1 = self.ser.read(4)
        # print(data1.hex())
        data1 = struct.unpack_from('f', data1, offset=0)[0]
        # print(data1)

        data2 = self.ser.read(4)
        # print(data2.hex())
        data2 = struct.unpack_from('f', data2, offset=0)[0]

        # print('x:')
        # print(data1)
        data3 = self.ser.read(4)
        # print(data3.hex())
        data3 = struct.unpack_from('f', data3, offset=0)[0]

        data4 = self.ser.read(4)
        data4 = struct.unpack_from('f', data4, offset=0)[0]
        x = np.array([data1, data2, data3, data4])
        # print(data4.hex())

        # elapsed = time.time() - self.t
        # time.sleep(0.5)
        # time.sleep(0.08 - elapsed)
        # print('freq:')
        # print(1 / (time.time() - self.t))
        return x

def plot_trajectories(x, u, dt, title=""):
    time = np.arange(0, x.shape[1]) * dt
    fig, axs = plt.subplots(5, 1)
    fig.suptitle(title)
    axs[0].plot(time, x[0, :], label='$\Theta$')
    axs[1].plot(time, x[1, :], label='$\dot{\Theta}$')
    axs[2].plot(time, x[2, :], label='$x$')
    axs[3].plot(time, x[3, :], label='$\dot{x}$')
    axs[4].plot(time[:-1], u, label='$u$')  #

    for ax in axs:
        ax.set_xlabel("$t$ [s]")
        ax.grid()
        ax.legend()
    print(f"theta at end of {title} = {x[0, -1]} ")


def measure(ser):
    x = ser.receive()
    return x


def apply(ser, u):
    ser.send(u)


with open('initialguessU.txt') as f:
    for line in f:
        currentline = line.split(",")
input = [float(i) for i in currentline]

# with open('initialguessX.txt') as f:
#     for line in f:
#         currentline = line.split(",")
# inputX = [float(i) for i in currentline]

ser = SerialCom('COM3', )
x=[]
u=[]
for i in range(len(input)-1000):
    t_ca = time.perf_counter()
    apply(ser, input[i])
    while  (0.08 - time.perf_counter() + t_ca)>=0:
        x_0 = measure(ser).reshape(4, 1)
        x.append(x_0)
        u.append(input[i])

x = np.concatenate(x, axis=1)
u = np.array(u[0:x.shape[1]-1])
plot_trajectories(x, u, 0.008, title="MPC")
plt.show()
