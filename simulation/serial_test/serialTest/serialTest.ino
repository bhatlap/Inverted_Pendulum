const byte numChars = 4;
char receivedChars[numChars];   // an array to store the received data

float u=0.0000;
float m = 0.3000;
float l = 0.4000;
float g = 9.8100;
float b = 0.0160;
float f1=0.0000;
float f2=0.0000;
float f3=0.0000;
float f4=0.0000;
float g1=0.0000;
float g2=0.0000;
float g3=0.0000;
float g4=0.0000;
float h1=0.0000;
float h2=0.0000;
float h3=0.0000;
float h4=0.0000;
float k1=0.0000;
float k2=0.0000;
float k3=0.0000;
float k4=0.0000;
float x1=3.1412;
float x2=0.0000;
float x3=0.0000;
float x4=0.0000;
double dt=0.0080;

boolean newData = false;
uint32_t x = 0;
float y =0.00;


union ArrayToInteger {
  byte array[4];
  uint32_t integer;
  float ffloat;
};

void setup() {
    Serial.begin(256000);
    //Serial.begin(250000);

}

void loop() {
    recvWithEndMarker();
    showNewData();
}

void recvWithEndMarker() {
    static byte ndx = 0;
    char endMarker = '\n';
    char rc;
    
    while (Serial.available() > 0 && newData == false) {
        rc = Serial.read();

        if (rc != endMarker || (rc == endMarker && ndx < numChars-1)) {
            receivedChars[ndx] = rc;
            ndx++;
            //if (ndx >= numChars) {
            //    ndx = numChars - 1;
            //}
        }
        else {
            receivedChars[ndx+1] = '\0'; // terminate the string
            ndx = 0;
            newData = true;
        }
    }
}

void showNewData() {
    if (newData == true) {

        ArrayToInteger converter; //Create a converter
//        converter.array[0] = receivedChars[0]; //save something to each byte in the array
//        converter.array[1] = receivedChars[1]; //save something to each byte in the array
//        converter.array[2] = receivedChars[2]; //save something to each byte in the array
//        converter.array[3] = receivedChars[3]; //save something to each byte in the array
        
        converter.array[0] = 0x00;//receivedChars[0]; //save something to each byte in the array
        converter.array[1] = 0x00;//receivedChars[1]; //save something to each byte in the array
        converter.array[2] = 0xb8;//receivedChars[2]; //save something to each byte in the array
        converter.array[3] = 0xc1;//receivedChars[3]; //save something to each byte in the array

        u = converter.ffloat; 
        f1 = x2;
        f2 = (3*g)/(2*l)*sin(x1)-(3*b)/(m*l*l)*x2+(3)/(2*l)*cos(x1)*u;
        f3 = x4;
        f4 = u;

        g1 = x2+dt/2*f2;
        g2 = (3*g)/(2*l)*sin(x1+dt/2*f1)-(3*b)/(m*l*l)*(x2+dt/2*f2)+(3)/(2*l)*cos(x1+dt/2*f1)*u;
        g3 = x4+dt/2*f4;
        g4 = u;

        h1 = x2+dt/2*g2;
        h2 = (3*g)/(2*l)*sin(x1+dt/2*g1)-(3*b)/(m*l*l)*(x2+dt/2*g2)+(3)/(2*l)*cos(x1+dt/2*g1)*u;
        h3 = x4+dt/2*g4;
        h4 = u;

        k1 = x2+dt*h2;
        k2 = (3*g)/(2*l)*sin(x1+dt*h1)-(3*b)/(m*l*l)*(x2+dt*h2)+(3)/(2*l)*cos(x1+dt*h1)*u;
        k3 = x4+dt*h4;
        k4 = u;
    
        x1= x1 + dt/6*(f1+2*g1+2*h1+k1);
        x2= x2 + dt/6*(f2+2*g2+2*h2+k2);
        x3= x3 + dt/6*(f3+2*g3+2*h3+k3);
        x4= x4 + dt/6*(f4+2*g4+2*h4+k4);

        ArrayToInteger backconverter_u; 
        ArrayToInteger backconverter_x1;
        ArrayToInteger backconverter_x2;
        ArrayToInteger backconverter_x3;
        ArrayToInteger backconverter_x4;

        backconverter_u.ffloat = u;
        Serial.write(backconverter_u.array,4);
//                Serial.print('\n'); 

        backconverter_x1.ffloat = x1;
        Serial.write(backconverter_x1.array,4);
 //          Serial.print('\n'); 
               backconverter_x2.ffloat = x2;
        Serial.write(backconverter_x2.array,4);
  //                   Serial.print('\n'); 

        backconverter_x3.ffloat = x3;
        // Serial.write(backconverter_x.array,4);
        Serial.write(backconverter_x3.array[0]);
        Serial.write(backconverter_x3.array[1]);
        Serial.write(backconverter_x3.array[2]);
        Serial.write(backconverter_x3.array[3]);
     //      Serial.print('\n'); 

        backconverter_x4.ffloat = x4;
        Serial.write(backconverter_x4.array,4);
        // Serial.print('\n'); 
          //         Serial.print('\n'); 

        newData = false;
    }
}
