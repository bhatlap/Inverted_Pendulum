#include <ESP32Encoder.h>
ESP32Encoder MotorEncoder;
ESP32Encoder PendulumEncoder;
#include "driver/ledc.h"

// CONFIGURE COMMUNICATION PARAMETERS FOR PWM CHANNEL
const int FREQUENCY = 2000000;
const int CHANNEL = 0;
const int RESOLUTION = 8;

// INITIALIZE MOTOR DRIVER PINS
const int PULSE_PIN = 17;
const int DIRECTION_PIN = 16;
const int ENABLE_PIN = 4;

// INITIALIZE LEFT- AND RIGHT-SWITCH PINS
const int LSW = 14;
const int RSW = 15;

// INITIALIZE MODEL PARAMETERS FOR SIMULATION
const float m = 0.3;      // mass of the pendulum [m]
const float g = 9.80665;  // gravitational acceleration [m / s^2]
const float l = 0.4;      // length of the pendulum [m]
const float b = 0.016;            // coefficient of friction

// INITIALIZE CONTROLLER PARAMETERS
float u_k = 0;  // current controller output  (m/s^2)
float a = 0.0;   // current constrained controller output (m / s^2)
uint32_t u = 0;  // current constrained controller output (Hz)

// INITIALIZE ENCODER VARIABLES
float pend_pulses = 0.0;  // pulses from the pendulum encoder
float mot_pulses = 0.0;   // pulses from the motor encoder

// INITIALIZE SYSTEM VARIABLES
float x_k = 0.0;   // current position of the cart   (m)
float x_k1 = 0.0;  // Previous position of the cart (m)

float th_k = 0.0;   // current angle of the pendulum  (rad)
float th_k1 = 0.0;  // previous angle of the pendulum (rad)

float xd_k = 0.0;   // current velocity of the cart   (mm / s)
float xd_k1 = 0.0;  // previous velocity of the cart  (mm / s)

float v_k = 0.0000;  // current cart velocity (m / s)

float thd_k = 0.0;  // current angular velocity of the pendulum (rad / s)
float thd_k1 = 0.0;  // first previous angular velocity of the pendulum (rad / s)
float thd_k2 = 0.0;  // second previous angular velocity of the pendulum (rad / s)
float thd_k3 = 0.0;  // third previous angular velocity of the pendulum (rad / s)

// INITIALIZE RK4 SIMULATATED SYSTEM 
float f1 = 0.0000;
float f2 = 0.0000;
float f3 = 0.0000;
float f4 = 0.0000;
float g1 = 0.0000;
float g2 = 0.0000;
float g3 = 0.0000;
float g4 = 0.0000;
float h1 = 0.0000;
float h2 = 0.0000;
float h3 = 0.0000;
float h4 = 0.0000;
float k1 = 0.0000;
float k2 = 0.0000;
float k3 = 0.0000;
float k4 = 0.0000;
float x1 = 3.1412;  
float x2 = 0.0000;
float x3 = 0.0000;
float x4 = 0.0000;

// INITIALIZE ITERATION TIMER
unsigned long start_timer;  
unsigned long current_timer;

// INITIALIZE CONTROL SEQUENCE FLAGS
boolean controller_on = false;      // flag to turn the controller on
boolean system_on = true;           // flag to turn the system on

// INITIALIZE TO RECEIVE INPUT
union ArrayToInteger 
{
  byte array[4];
  float ffloat;
};
ArrayToInteger input;
boolean newData = false;
const byte numChars = 4;
char receivedChars[numChars];  // an array to store the received data

// SIMULATION FLAG
const boolean Simulation = false;

void setup() {
  Serial.begin(115200);
  delay(20);
  Serial.setTimeout(5);

  input.ffloat = 0;
  u_k = 0;
  if (Simulation == false) {

    // ATTACH PWM CHANNEL TO PWM PIN
    ledcSetup(CHANNEL, FREQUENCY, RESOLUTION);
    ledcAttachPin(PULSE_PIN, CHANNEL);

    // SETUP ENCODERS                             => LOOKING FROM THE PENDULUM'S SIDE
    MotorEncoder.attachFullQuad(27, 26);     // 27-26 ==> RIGHT => positive, LEFT => negative
    PendulumEncoder.attachFullQuad(33, 32);  // 33-32 ==> CCW => degrees increase, CW => degrees decrease

    // SETUP MODES FOR THE DIRECTION AND ENABLE PINS
    pinMode(DIRECTION_PIN, OUTPUT);
    pinMode(ENABLE_PIN, OUTPUT);

    // WRITE INITIAL VALUES TO DIRECTION AND ENABLE PINS
    digitalWrite(DIRECTION_PIN, HIGH);
    digitalWrite(ENABLE_PIN, LOW);

    if (system_on) {
      delay(1000);

      controller_on = true;
      system_on = false;

      thd_k1 = 0;
      thd_k2 = 0;
      thd_k3 = 0;

      pend_pulses = PendulumEncoder.getCount();
      mot_pulses = MotorEncoder.getCount();

      th_k = ((float)pend_pulses + 800) * 2.0 * 3.14159 / 1600.0;
      th_k1 = ((float)pend_pulses + 800) * 2.0 * 3.14159 / 1600.0;
      
      x_k = ((float)mot_pulses) / 20 / 1000;
      x_k1 = ((float)mot_pulses) / 20 / 1000;

      xd_k = 0;
      xd_k1 = 0;


      start_timer = esp_timer_get_time();
    }
  }
  else 
  {
    th_k = x1;
    thd_k = x2;
    x_k = x3;
    xd_k = x4;
    v_k = x4;

    start_timer = esp_timer_get_time();
  }

}  


void loop()
{

  read();

  if (Simulation == false) 
  {
    // COUNT PULSES FROM ENCODERS
    pend_pulses = PendulumEncoder.getCount();
    mot_pulses = MotorEncoder.getCount();

    // CALCULATE CURRENT PENDULUM ANGLE AND CART POSITION
    th_k = ((float)pend_pulses + 800) * 2.0 * 3.14159 / 1600.0;
    x_k = ((float)mot_pulses) / 20 / 1000;
  }
  
  current_timer = esp_timer_get_time();
  if ((current_timer - start_timer) <= (5000))
  {
    current_timer = esp_timer_get_time();  // WAITING
  } 
  else 
  {
    if (Simulation == false) 
    {
      // PENDULUM ANGULAR VELOCITY
      thd_k = (th_k - th_k1) / 0.005;
      th_k1 = th_k;

      thd_k = (thd_k + thd_k1 + thd_k2 + thd_k3) / 4;
      thd_k1 = thd_k;
      thd_k2 = thd_k1;
      thd_k3 = thd_k2;

      // PENDULUM LINEAR VELOCITY
      v_k = (x_k - x_k1)/0.005 + a*0.005/2;
      x_k1 = x_k;
    }
    sendStates();

    if (Simulation == false) 
    {
      real_dynamics();
    }
    else 
    {
      simulated_dynamics();
    }

    start_timer = current_timer;
  }
}


void recvWithEndMarker() 
{

  static byte ndx = 0;
  char endMarker = '\n';
  char rc;

  while (Serial.available() > 0 && newData == false) 
  {
    rc = Serial.read();
    if (rc != endMarker) 
    {
      receivedChars[ndx] = rc;
      ndx++;
    } 
    else 
    {
      receivedChars[ndx + 1] = '\0';  // terminate the string
      ndx = 0;
      newData = true;
    }
  }
}

void storeNewData() 
{
  if (newData == true) 
  {
    input.array[0] = receivedChars[0];  //save something to each byte in the array
    input.array[1] = receivedChars[1];  //save something to each byte in the array
    input.array[2] = receivedChars[2];  //save something to each byte in the array
    input.array[3] = receivedChars[3];  //save something to each byte in the array
    newData = false;
    u_k = input.ffloat;
  }
}

void read()
{
  recvWithEndMarker();
  storeNewData();
}

void sendStates() 
{
  ArrayToInteger backconverter_u;
  ArrayToInteger backconverter_x1;
  ArrayToInteger backconverter_x2;
  ArrayToInteger backconverter_x3;
  ArrayToInteger backconverter_x4;

  Serial.write('?');
  backconverter_u.ffloat = u_k;
  Serial.write(backconverter_u.array, 4);

  backconverter_x1.ffloat = th_k;
  Serial.write(backconverter_x1.array, 4);

  backconverter_x2.ffloat = thd_k;
  Serial.write(backconverter_x2.array, 4);

  backconverter_x3.ffloat = x_k;
  Serial.write(backconverter_x3.array, 4);

  backconverter_x4.ffloat = v_k; //Alterantive velocity = xd_k/(Simulation ? 1 : 1000)
  Serial.write(backconverter_x4.array, 4);
}

void simulated_dynamics() 
{
  float dT = 0.005; 
  x1 = th_k;
  x2 = thd_k;
  x3 = x_k;
  x4 = xd_k;

  f1 = x2;
  f2 = (3 * g) / (2 * l) * sin(x1) - (3 * b) / (m * l * l) * x2 + (3) / (2 * l) * cos(x1) * u_k;
  f3 = x4;
  f4 = u_k;

  g1 = x2 + dT / 2 * f2;
  g2 = (3 * g) / (2 * l) * sin(x1 + dT / 2 * f1) - (3 * b) / (m * l * l) * (x2 + dT / 2 * f2) + (3) / (2 * l) * cos(x1 + dT / 2 * f1) * u_k;
  g3 = x4 + dT / 2 * f4;
  g4 = u_k;

  h1 = x2 + dT / 2 * g2;
  h2 = (3 * g) / (2 * l) * sin(x1 + dT / 2 * g1) - (3 * b) / (m * l * l) * (x2 + dT / 2 * g2) + (3) / (2 * l) * cos(x1 + dT / 2 * g1) * u_k;
  h3 = x4 + dT / 2 * g4;
  h4 = u_k;

  k1 = x2 + dT * h2;
  k2 = (3 * g) / (2 * l) * sin(x1 + dT * h1) - (3 * b) / (m * l * l) * (x2 + dT * h2) + (3) / (2 * l) * cos(x1 + dT * h1) * u_k;
  k3 = x4 + dT * h4;
  k4 = u_k;

  x1 = x1 + dT / 6 * (f1 + 2 * g1 + 2 * h1 + k1);
  x2 = x2 + dT / 6 * (f2 + 2 * g2 + 2 * h2 + k2);
  x3 = x3 + dT / 6 * (f3 + 2 * g3 + 2 * h3 + k3);
  x4 = x4 + dT / 6 * (f4 + 2 * g4 + 2 * h4 + k4);

  th_k = x1;
  thd_k = x2;
  x_k = x3;
  xd_k = x4;
  v_k = x4;
}

void real_dynamics() 
{
  if (controller_on)
  {
    if (!digitalRead(LSW) && !digitalRead(RSW)) 
    {
      
      // COMPUTE CURRENT CONSTRAINED CONTROLLER OUTPUT (AS ACCELERATION)
      
      a = constrain(u_k, -25, 25);
      
      xd_k = a * 0.005 * 1000 + xd_k1;
      xd_k1 = xd_k;


      // WRITE ON DIRECTION PIN TO UPDATE IF VELOCITY HAS CHANGED
      digitalWrite(DIRECTION_PIN, (xd_k > 0 ? 0 : 1));
      
  
      uint32_t u1 = u;  // PREVIOUS INPUT VALUE FOR COMPARISON

      // COMPUTE CURRENT CONSTRAINED CONTROLLER OUTPUT (AS FREQUENCY)
      u = constrain(abs(round(xd_k * 20)), 0, FREQUENCY);  //2000000);
      
      if ((u1-u)>0)     
      {
        ledcSetup(CHANNEL, u, 9); //ledcSetup sets the frequency of the PWM signal 
        ledcWrite(CHANNEL, 0xFF); //ledcWrite controlls the duty cycle of the PWM signal and start the signal
      }
    }
  }
}