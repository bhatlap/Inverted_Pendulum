import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from casadi import *
import time
from statistics import mean
import scipy.linalg as lin
from scipy import signal
import serial
from pySerialTransfer import pySerialTransfer as txfer

import struct

import psutil
import os

p = psutil.Process(os.getpid())
p.nice(psutil.HIGH_PRIORITY_CLASS)

class sendStruct(object):
    u = 0.0
    x1 = 0.0
    x2 = 0.0
    x3 = 0.0
    x4 = 0.0

class SerialCom:
    def __init__(self, port):
        self.link = txfer.SerialTransfer(port)
        self.link.open()
        self.t = 0

    def send(self, u):
        send_size = 0
        self.t = time.time()
        ###################################################################
        # Send a float
        ###################################################################
        float_ = u.item()
        send_size = self.link.tx_obj(float_, start_pos=send_size)
        ###################################################################
        # Transmit all the data to send in a single packet
        ###################################################################
        self.link.send(send_size)

    def receive(self):
        ###################################################################
        # Wait for a response and report any errors while receiving packets
        ###################################################################
        while not self.link.available():
            if self.link.status < 0:
                if self.link.status == txfer.CRC_ERROR:
                    print('ERROR: CRC_ERROR')
                elif self.link.status == txfer.PAYLOAD_ERROR:
                    print('ERROR: PAYLOAD_ERROR')
                elif self.link.status == txfer.STOP_BYTE_ERROR:
                    print('ERROR: STOP_BYTE_ERROR')
                else:
                    print('ERROR: {}'.format(self.link.status))

        ###################################################################
        # Parse response float
        ###################################################################
        recSize = 0
        sendStruct = struct
        sendStruct.u = self.link.rx_obj(obj_type='f', start_pos=recSize)
        recSize += txfer.STRUCT_FORMAT_LENGTHS['f']

        sendStruct.x1 = self.link.rx_obj(obj_type='f', start_pos=recSize)
        recSize += txfer.STRUCT_FORMAT_LENGTHS['f']

        sendStruct.x2 = self.link.rx_obj(obj_type='f', start_pos=recSize)
        recSize += txfer.STRUCT_FORMAT_LENGTHS['f']

        sendStruct.x3 = self.link.rx_obj(obj_type='f', start_pos=recSize)
        recSize += txfer.STRUCT_FORMAT_LENGTHS['f']

        sendStruct.x4 = self.link.rx_obj(obj_type='f', start_pos=recSize)
        recSize += txfer.STRUCT_FORMAT_LENGTHS['f']

        x = np.array([sendStruct.x1, sendStruct.x2, sendStruct.x3, sendStruct.x4])

        # elapsed = time.time() - self.t
        # time.sleep(0.5)
        # time.sleep(0.08 - elapsed)
        # print('freq:')
        # print(1 / (time.time() - self.t))
        return x


def measure(ser):
    x = ser.receive()
    return x


def apply(ser, u):
    ser.send(u)


