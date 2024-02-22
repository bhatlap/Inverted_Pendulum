#define sign(x) ((x) < 0 ? -1 : ((x) > 0 ? 1 : 0))

#include <ESP32Encoder.h>
ESP32Encoder MotorEncoder;
ESP32Encoder PendulumEncoder;

// CREATE A HARDWARE TIMER  (CLOCK FREQUENCY => 40 MHZ)
hw_timer_t * timer = NULL;                                                  

// CONFIGURE COMMUNICATION PARAMETERS FOR PWM CHANNEL 
const int FREQUENCY  = 2000;
const int CHANNEL    = 0;
const int RESOLUTION = 8;

// INITIALIZE MOTOR DRIVER PINS
const int PULSE_PIN     = 17;
const int DIRECTION_PIN = 16;
const int ENABLE_PIN    = 4;

// INITIALIZE LEFT- AND RIGHT-SWITCH PINS
const int LSW = 14;
const int RSW = 15;

// INITIALIZE MODEL PARAMETERS
const float m = 0.3;                // mass of the pendulum [m]
const float g = 9.80665;            // gravitational acceleration [m / s^2]
const float l = 0.4;                // length of the pendulum [m]

const float b = 0.016;              // coefficient of friction 
const float J = m*(l^2)/3;          // moment of inertia of the pendulum

const float E_star = m*g*l;         // reference energy for the upright position
const int x_star = 0;               // position to be tracked for the position controller

// INITIALIZE CONTROLLER PARAMETERS
const float dT   = 0.008; // sampling time (s)
const int KE     = 3;     // controller gain for energy error
const int KC     = 20;    // controller gain for position error

float u_k  = 0;       // current controller output  (m/s^2)

float a    = 0.0;     // current constrained controller output (m / s^2)
uint32_t u = 0;       // current constrained controller output (Hz)

// INITIALIZE SYSTEM VARIABLES
float x_k    = 0.0;   // current position of the cart   (m)

float th_k   = 0.0;   // current angle of the pendulum  (rad)
float th_k1  = 0.0;   // previous angle of the pendulum (rad)

float xd_k   = 0.0;   // current velocity of the cart   (mm / s)
float xd_k1  = 0.0;   // previous velocity of the cart  (mm / s)

float thd_k  = 0.0;   // current angular velocity of the pendulum

float thd_k1 = 0.0;   // first previous angular velocity of the pendulum
float thd_k2 = 0.0;   // second previous angular velocity of the pendulum
float thd_k3 = 0.0;   // third previous angular velocity of the pendulum
float thd_k4 = 0.0;   // fourth previous angular velocity of the pendulum

// INITIALIZE ENCODER VARIABLES
float pend_pulses = 0.0;  // pulses from the pendulum encoder
float mot_pulses = 0.0;   // pulses from the motor encoder

// INITIALIZE CONTROL SEQUENCE FLAGS
boolean controller_on = false;     // flag to turn the controller on
boolean system_on = true;          // flag to turn the system on
volatile boolean timer_on = false; // flag to turn the timer on

void IRAM_ATTR onTimer(){
  timer_on = !timer_on;
}

void setup() {
  Serial.begin(115200);
  delay(20);
  Serial.setTimeout(5);

  // SETUP HARDWARE TIMER
  timer = timerBegin(0, 40, true);
  timerAttachInterrupt(timer, &onTimer, true);
  timerAlarmWrite(timer, 8000, true);
  timerAlarmEnable(timer);

  // ATTACH PWM CHANNEL TO PWM PIN
  ledcSetup(CHANNEL, FREQUENCY, RESOLUTION);
  ledcAttachPin(PULSE_PIN, CHANNEL);

  // SETUP ENCODERS                             => LOOKING FROM THE PENDULUM'S SIDE
  MotorEncoder.attachFullQuad(27, 26);          // 27-26 ==> RIGHT => positive, LEFT => negative 
  PendulumEncoder.attachFullQuad(33, 32);       // 33-32 ==> CCW => degrees increase, CW => degrees decrease

  // SETUP MODES FOR THE DIRECTION AND ENABLE PINS
  pinMode(DIRECTION_PIN, OUTPUT);
  pinMode(ENABLE_PIN, OUTPUT);

  // WRITE INITIAL VALUES TO DIRECTION AND ENABLE PINS
  digitalWrite(DIRECTION_PIN, HIGH);
  digitalWrite(ENABLE_PIN, LOW);
  
}

void loop() {

    // COUNT PULSES FROM ENCODERS
    pend_pulses = PendulumEncoder.getCount();
    mot_pulses  = MotorEncoder.getCount();
    
    // CALCULATE CURRENT PENDULUM ANGLE AND CART POSITION
    th_k = ((float)pend_pulses + 800)*2.0*3.14159/1600.0;   
    x_k  = ((float)mot_pulses)/20/1000;                       

    Serial.println("th: " + String(th_k) + " x: " + String(x_k));
    if(system_on){
      delay(3000);
      
      controller_on = true;
      system_on = false;

      thd_k1 = 0;
      thd_k2 = 0;
      thd_k3 = 0;
      thd_k4 = 0;
      
      th_k1  = 0;
      
      xd_k   = 0;
      xd_k1  = 0;
      
    }
    
    if(timer_on and controller_on){
      if( !digitalRead(LSW) && !digitalRead(RSW) ){     
          if(abs(th_k) > 0.1745){
            
            // COMPUTE PENDULUM ANGULAR VELOCITY BY AVERAGING OVER DERIVATIVE APPROXIMATIONS
            thd_k = (th_k - th_k1) / dT;
            thd_k = (thd_k + thd_k1 + thd_k2 + thd_k3 + thd_k4) / 5;

            thd_k1 = thd_k;
            thd_k2 = thd_k1;
            thd_k3 = thd_k2;
            thd_k4 = thd_k3;

            // COMPUTE PENDULUM ENERGY AND ERRORS
            E_pendulum = 0.5*J*(thd_k^2) + m*g*(l/2)*(cos(th_k) + 1);
            EE = E_star - E_pendulum;
            PE = x_star - x_k;

            // COMPUTE CURRENT CONTROLLER OUTPUT (AS ACCELERATION)
            if EE*thd_k*cos(th_k) == 0
              u_k = KE + KC*PE;
            else
              u_k = KE*sign(EE*thd*cos(th_k)) + KC*PE;
            end        

            // COMPUTE CURRENT CONSTRAINED CONTROLLER OUTPUT (AS ACCELERATION)
            th_k1 = th_k;
            a = constrain(u_k, -25, 25);
            
            // COMPUTE CURRENT CART VELOCITY TO BE APPLIED AS INPUT
            xd_k = a*dT*1000 + xd_k1;
            
            // WRITE ON DIRECTION PIN TO UPDATE IF VELOCITY HAS CHANGED 
            digitalWrite(DIRECTION_PIN, (xd_k>0?0:1));
            xd_k1 = xd_k;
            
            // COMPUTE CURRENT CONSTRAINED CONTROLLER OUTPUT (AS FREQUENCY)
            u = constrain(abs(round(xd_k*40)), 0, 2000000);

            // APPLY INPUT
            ledcWriteTone(CHANNEL, u);   // ledcWriteTone sets the frequency of the PWM signal
                                         // ledcWrite controlls the duty cycle of the PWM signal
            Serial.println(u);
            
            // DEBUG LOOP
            /*
            Serial.println("th: " + String(th_k) + " thd: " + String(thd_k) + " x: " + String(x_k) + " xd: " + String(xd_k));
            Serial.println("thd_k1: " + String(thd_k1) + " thd_k2: " + String(thd_k2) + " thd_k3: " + String(thd_k3) + " thd_k4: " + String(thd_k4));
            Serial.println("u_k: " + String(u_k) + " u: " + String(u));
            */
            
          }else{
            controller_on = false;
          }
      }else{
        ledcWriteTone(CHANNEL, 0);
        digitalWrite(PULSE_PIN, LOW);
      }
    
    timer_on = false;
  }
 
}