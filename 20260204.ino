#include <Wire.h>
#include <SD.h>      
#include <math.h>

// ================================================================
// [SECTION 0] 사용자 설정 및 튜닝된 계수
// ================================================================
const int buzzerPin = 37;    
const float LAUNCH_THRESHOLD_G = 3.0f; 

// [튜닝된 필터 계수]
#define MADGWICK_BETA 0.3f 
float Q_pos = 0.0003800563f;
float Q_vel = 0.0380056300f;
float R_meas = 0.0755295250f;

// [보정 임계값]
const int descent_check_limit = 500; 
const int landing_check_limit = 2000; 
const float MAX_BARO_INNOVATION = 15.0f; // 화약 사출 시 15m 이상의 기압 튐 무시

// = [구조체] SD카드 저장 데이터 형식 =
struct LogData 
{
  uint32_t time_ms;
  float ax, ay, az;      
  float gx, gy, gz;      
  float baro_alt;        
  float est_alt;         
  float est_vel;         
  float vertical_acc;    
  uint8_t event_status;  // 0:보정, 1:준비, 2:상승, 3:하강, 4:종료
};

// ================================================================
// [SECTION 1] Filter Class & Global variations
// ================================================================
class RocketMadgwick 
{
public:
    float q0, q1, q2, q3; 
    float beta;
    float invSampleFreq;

    RocketMadgwick() : q0(1.0f), q1(0.0f), q2(0.0f), q3(0.0f), beta(MADGWICK_BETA), invSampleFreq(0.001f) {}
    void begin(float sampleFrequency) 
    {
        invSampleFreq = 1.0f / sampleFrequency;
    }
 
    float invSqrt(float x) 
    {
        if (x == 0.0f) return 0.0f;
        float halfx = 0.5f * x;
        float y = x;
        long i = *(long*)&y;
        i = 0x5f3759df - (i >> 1);
        y = *(float*)&i;
        y = y * (1.5f - (halfx * y * y));
        return y;
    }
    
    void update(float gx, float gy, float gz, float ax, float ay, float az) 
    {
        float recipNorm;
        float s0, s1, s2, s3;
        float qDot1, qDot2, qDot3, qDot4;
        float _2q0, _2q1, _2q2, _2q3, _4q0, _4q1, _4q2 ,_8q1, _8q2, q0q0, q1q1, q2q2, q3q3;

        gx *= 0.0174533f; 
        gy *= 0.0174533f; 
        gz *= 0.0174533f;
        qDot1 = 0.5f * (-q1 * gx - q2 * gy - q3 * gz);
        qDot2 = 0.5f * (q0 * gx + q2 * gz - q3 * gy);
        qDot3 = 0.5f * (q0 * gy - q1 * gz + q3 * gx);
        qDot4 = 0.5f * (q0 * gz + q1 * gy - q2 * gx);

        if(!((ax == 0.0f) && (ay == 0.0f) && (az == 0.0f))) 
        {
            recipNorm = invSqrt(ax * ax + ay * ay + az * az);
            ax *= recipNorm; 
            ay *= recipNorm; 
            az *= recipNorm;   
            _2q0 = 2.0f * q0; 
            _2q1 = 2.0f * q1; 
            _2q2 = 2.0f * q2; 
            _2q3 = 2.0f * q3;
            _4q0 = 4.0f * q0; 
            _4q1 = 4.0f * q1; 
            _4q2 = 4.0f * q2;
            _8q1 = 8.0f * q1; 
            _8q2 = 8.0f * q2;
            q0q0 = q0 * q0; 
            q1q1 = q1 * q1; 
            q2q2 = q2 * q2; 
            q3q3 = q3 * q3;

            s0 = _4q0 * q2q2 + _2q2 * ax + _4q0 * q1q1 - _2q1 * ay;
            s1 = _4q1 * q3q3 - _2q3 * ax + 4.0f * q0q0 * q1 - _2q0 * ay - _4q1 + _8q1 * q1q1 + _8q1 * q2q2 + _4q1 * az;
            s2 = 4.0f * q0q0 * q2 + _2q0 * ax + _4q2 * q3q3 - _2q3 * ay - _4q2 + _8q2 * q1q1 + _8q2 * q2q2 + _4q2 * az;
            s3 = 4.0f * q1q1 * q3 - _2q1 * ax + 4.0f * q2q2 * q3 - _2q2 * ay;
            recipNorm = invSqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3);
            s0 *= recipNorm; 
            s1 *= recipNorm; 
            s2 *= recipNorm; 
            s3 *= recipNorm;
            qDot1 -= beta * s0; 
            qDot2 -= beta * s1; 
            qDot3 -= beta * s2; 
            qDot4 -= beta * s3;
        }
        q0 += qDot1 * invSampleFreq; 
        q1 += qDot2 * invSampleFreq;
        q2 += qDot3 * invSampleFreq; 
        q3 += qDot4 * invSampleFreq;
        recipNorm = invSqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
        q0 *= recipNorm; 
        q1 *= recipNorm; 
        q2 *= recipNorm; 
        q3 *= recipNorm;
    }
};

RocketMadgwick filter;
File dataFile;
char fileName[16];

unsigned long LOOP_INTERVAL_US = 1000; 
unsigned long nextLoopTime_us = 0;
unsigned long startTime_ms = 0;
unsigned long lastFlush_ms = 0;

float est_altitude = 0.0f;
float est_velocity = 0.0f; 
float P[2][2] = {{1, 0}, {0, 1}}; 
float gyro_offset[3] = {0,0,0};
float filt_acc[3] = {0,0,0};

#define ICM_Addr 0x68
#define MS5611_Addr 0x77
#define SEALEVEL 1013.25
uint16_t C[7]; 
const unsigned long BARO_INTERVAL_US = 50000; 
unsigned long nextBaroTime_us = 0; 
unsigned long baroWaitStart_us = 0;
enum BaroState { BARO_IDLE, BARO_WAIT_D1, BARO_WAIT_D2 };
BaroState baroState = BARO_IDLE;
uint32_t D1 = 0, D2 = 0;
float current_baro_alt = 0.0f;
bool newBaroData = false;

float ground_alt = 0.0f;
float acc_bias = 0.0f;
bool isCalibrated = false;
float altSum = 0, accSum = 0;
int calibSamples = 0;

uint8_t current_status = 0;
bool isLaunched = false;
unsigned long launchStartTime = 0;
float max_altitude = -9999.0f;
int descent_counter = 0;
int landing_counter = 0;

// ================================================================
// [SECTION 2] Helper Functions
// ================================================================

void writeReg(uint8_t reg, uint8_t data) 
{
    Wire.beginTransmission(ICM_Addr);
    Wire.write(reg); 
    Wire.write(data);
    Wire.endTransmission();
}

void setBank(uint8_t bank) 
{
    writeReg(0x7F, bank << 4);
}

void ms5611_cmd(uint8_t cmd) 
{
    Wire1.beginTransmission(MS5611_Addr);
    Wire1.write(cmd); 
    Wire1.endTransmission();
}

uint16_t ms5611_readCoeff(uint8_t idx) 
{
    Wire1.beginTransmission(MS5611_Addr);
    Wire1.write(0xA0 + idx * 2); 
    Wire1.endTransmission(false);
    Wire1.requestFrom(MS5611_Addr, 2, true);
    if(Wire1.available()>=2)
    {
        return (Wire1.read()<<8)|Wire1.read();
    }
    return 0;
}

uint32_t ms5611_readADC() 
{
    Wire1.beginTransmission(MS5611_Addr);
    Wire1.write(0x00); 
    Wire1.endTransmission(false);
    Wire1.requestFrom(MS5611_Addr, 3, true);
    if(Wire1.available()>=3) 
    {
        uint32_t h=Wire1.read(); 
        uint32_t m=Wire1.read(); 
        uint32_t l=Wire1.read();
        return (h<<16)|(m<<8) | l;
    }
    return 0;
}

void updateBarometer() 
{
    unsigned long now = micros();
    if (baroState == BARO_IDLE && now >= nextBaroTime_us) 
    {
        ms5611_cmd(0x48); 
        baroWaitStart_us = now; 
        baroState = BARO_WAIT_D1;
        nextBaroTime_us = now + BARO_INTERVAL_US;
    }
    if (baroState == BARO_WAIT_D1 && (now - baroWaitStart_us >= 9500)) 
    {
        D1 = ms5611_readADC(); 
        ms5611_cmd(0x58); 
        baroWaitStart_us = now; 
        baroState = BARO_WAIT_D2;
    }
    else if (baroState == BARO_WAIT_D2 && (now - baroWaitStart_us >= 9500)) 
    {
        D2 = ms5611_readADC();
        int64_t dT = D2 - ((int64_t)C[5] << 8);
        int64_t TEMP = 2000 + ((dT * C[6]) >> 23);
        int64_t OFF = ((int64_t)C[2] << 16) + ((C[4] * dT) >> 7);
        int64_t SENS = ((int64_t)C[1] << 15) + ((C[3] * dT) >> 8);
        int32_t P_val = (((D1 * SENS) >> 21) - OFF) >> 15;
        current_baro_alt = 44330.0f * (1.0f - pow(((float)P_val/100.0f) / SEALEVEL, 0.190295f));
        newBaroData = true; 
        baroState = BARO_IDLE;
    }
}

void kalmanUpdate(float acc_z, float meas_alt, float dt) 
{
    est_altitude += est_velocity * dt + 0.5f * acc_z * dt * dt;
    est_velocity += acc_z * dt;
    float dt2 = dt * dt;
    P[0][0] += dt * (P[1][0] + P[0][1]) + dt2 * P[1][1] + Q_pos;
    P[0][1] += dt * P[1][1]; 
    P[1][0] += dt * P[1][1]; 
    P[1][1] += dt2 * Q_vel;
    if (newBaroData) 
    {
        float y = meas_alt - est_altitude;
        if (abs(y) > MAX_BARO_INNOVATION) {
            y = (y > 0) ? MAX_BARO_INNOVATION : -MAX_BARO_INNOVATION;
        }
        float S = P[0][0] + R_meas;
        float K0 = P[0][0] / S; 
        float K1 = P[1][0] / S;
        est_altitude += K0 * y; 
        est_velocity += K1 * y;
        float P00 = P[0][0]; 
        float P01 = P[0][1];
        P[0][0] -= K0 * P00; 
        P[0][1] -= K0 * P01;
        P[1][0] -= K1 * P00; 
        P[1][1] -= K1 * P01;
        newBaroData = false;
    }
}

void gyroCalib() 
{
    digitalWrite(buzzerPin, HIGH); 
    delay(100); 
    digitalWrite(buzzerPin, LOW); 
    delay(500);
    long s[3] = {0,0,0};
    for(int i=0; i<1000; i++) 
    {
        Wire.beginTransmission(ICM_Addr); 
        Wire.write(0x33); 
        Wire.endTransmission(false);
        Wire.requestFrom(ICM_Addr, 6, true);
        s[0] += (int16_t)(Wire.read()<<8|Wire.read());
        s[1] += (int16_t)(Wire.read()<<8|Wire.read());
        s[2] += (int16_t)(Wire.read()<<8|Wire.read());
        delay(1);
    }
    for(int i=0; i<3; i++)
    {
        gyro_offset[i] = (float)s[i]/1000.0/8.2;
    }
    digitalWrite(buzzerPin, HIGH); 
    delay(100); 
    digitalWrite(buzzerPin, LOW);
}

// ================================================================
// [SECTION 3] Setup
// ================================================================
void setup() 
{
    pinMode(buzzerPin, OUTPUT);

    Wire.begin(); 
    Wire.setClock(400000);
    Wire1.begin(); 
    Wire1.setClock(400000);

    setBank(0); 
    delay(10); 
    writeReg(0x06, 0x01); 
    delay(10); 
    setBank(2); 
    delay(10);
    writeReg(0x14, 0x1F); 
    writeReg(0x01, 0x1B); 
    setBank(0);

    ms5611_cmd(0x1E); 
    delay(100);
    for(int i=1;i<=6;i++) 
    {
        C[i] = ms5611_readCoeff(i);
        if(C[i] == 0 || C[i] == 0xFFFF) 
        {
            while(1) 
            { 
                digitalWrite(buzzerPin, HIGH); delay(50); 
                digitalWrite(buzzerPin, LOW); delay(50); 
            }
        }
    }

    if(!SD.begin(BUILTIN_SDCARD)) 
    { 
        while(1) 
        { 
            digitalWrite(buzzerPin, HIGH); delay(100); 
            digitalWrite(buzzerPin, LOW); delay(100); 
        } 
    }
    
    for (int i = 0; i < 100; i++) 
    {
        sprintf(fileName, "LOG_%02d.BIN", i);
        if (!SD.exists(fileName)) break;
    }
    dataFile = SD.open(fileName, FILE_WRITE);

    gyroCalib();
    filter.begin(1000); 

    for(int i=0; i<3; i++) 
    { 
        digitalWrite(buzzerPin, HIGH); delay(50); 
        digitalWrite(buzzerPin, LOW); delay(50); 
    }
    
    startTime_ms = millis();
    nextLoopTime_us = micros();
}

// ================================================================
// [SECTION 4] Main Loop
// ================================================================
void loop() 
{
    unsigned long now_us = micros();
    if (now_us < nextLoopTime_us) return;
    nextLoopTime_us += LOOP_INTERVAL_US;

    // 1. 센서 데이터 취득
    Wire.beginTransmission(ICM_Addr);
    Wire.write(0x2D);
    Wire.endTransmission(false);
    Wire.requestFrom(ICM_Addr, 12, true);
    int16_t ax_raw = (Wire.read()<<8)|Wire.read();
    int16_t ay_raw = (Wire.read()<<8)|Wire.read();
    int16_t az_raw = (Wire.read()<<8)|Wire.read();
    int16_t gx_raw = (Wire.read()<<8)|Wire.read();
    int16_t gy_raw = (Wire.read()<<8)|Wire.read();
    int16_t gz_raw = (Wire.read()<<8)|Wire.read();

    updateBarometer();

    // 2. 데이터 보정 및 필터 업데이트
    float ax = (float)ax_raw / 1024.0f;
    float ay = (float)ay_raw / 1024.0f;
    float az = (float)az_raw / 1024.0f;
    filt_acc[0] = filt_acc[0]*0.8f + ax*0.2f;
    filt_acc[1] = filt_acc[1]*0.8f + ay*0.2f;
    filt_acc[2] = filt_acc[2]*0.8f + az*0.2f;

    float gx_f = ((float)gx_raw/8.2f) - gyro_offset[0];
    float gy_f = ((float)gy_raw/8.2f) - gyro_offset[1];
    float gz_f = ((float)gz_raw/8.2f) - gyro_offset[2];

    filter.update(gx_f, gy_f, gz_f, filt_acc[0], filt_acc[1], filt_acc[2]);
    
    float vz_x = 2.0f*(filter.q1*filter.q3 - filter.q0*filter.q2);
    float vz_y = 2.0f*(filter.q0*filter.q1 + filter.q2*filter.q3);
    float vz_z = filter.q0*filter.q0 - filter.q1*filter.q1 - filter.q2*filter.q2 + filter.q3*filter.q3;
    float final_acc_g = (filt_acc[0]*vz_x + filt_acc[1]*vz_y + filt_acc[2]*vz_z) - 1.0f;
    float final_acc_ms2 = final_acc_g * 9.80665f;

    // 3. 지상 영점 보정 로직
    unsigned long elapsed = millis() - startTime_ms;
    if (!isCalibrated) 
    {
        if (elapsed > 8000 && elapsed < 10000) 
        {
            altSum += current_baro_alt;
            accSum += final_acc_ms2;
            calibSamples++;
        } 
        else if (elapsed >= 10000) 
        {
            ground_alt = altSum / (float)calibSamples;
            acc_bias = accSum / (float)calibSamples;
            isCalibrated = true;
            current_status = 1;
            for(int abc=0;abc<10;abc++)
            {
                digitalWrite(buzzerPin, HIGH); delay(50);
                digitalWrite(buzzerPin, LOW); delay(50);
            }
            digitalWrite(buzzerPin, HIGH); delay(500); 
            digitalWrite(buzzerPin, LOW);
        }
    }

    // 4. 비행 상태 관리 (State Machine)
    float corrected_acc = 0, corrected_alt = 0;
    if (isCalibrated) 
    {
        corrected_acc = final_acc_ms2 - acc_bias;
        corrected_alt = current_baro_alt - ground_alt;
        kalmanUpdate(corrected_acc, corrected_alt, 0.001f);
        
        unsigned long current_time = millis();

        // [State 1 -> 2] 발사 감지
        if (current_status == 1 && corrected_acc > (LAUNCH_THRESHOLD_G * 9.80665f)) 
        {
            isLaunched = true;
            launchStartTime = current_time;
            current_status = 2;
        } 
        // [State 2 -> 3] 정점 및 하강 감지 (이지미니 화약 교란 방어형)
        else if (current_status == 2) 
        {
            if (est_altitude > max_altitude) max_altitude = est_altitude;
            
            bool isSafetyOff = (current_time - launchStartTime) > 2000;
            if (isSafetyOff) 
            {
                bool isPastApogee = est_altitude < (max_altitude - 2.0f);
                bool isDescending = est_velocity < -2.0f;
                
                if (isPastApogee || isDescending) {
                    descent_counter++;
                } else {
                    if (descent_counter > 0) descent_counter--; 
                }

                if (descent_counter >= descent_check_limit) {
                    current_status = 3;
                    digitalWrite(buzzerPin, HIGH); delay(100); digitalWrite(buzzerPin, LOW);
                }
            }
        }
        // [State 3 -> 4] 착륙 판정
        else if (current_status == 3) 
        {
            if (est_altitude < 15.0f && abs(est_velocity) < 1.0f) {
                landing_counter++;
            } else {
                landing_counter = 0;
            }

            if (landing_counter >= landing_check_limit) {
                current_status = 4;
            }
        }
    }

    // 5. 로깅 및 사후 처리
    if (dataFile && current_status <= 4) 
    {
        LogData log;
        log.time_ms = millis() - startTime_ms;
        log.ax = ax; log.ay = ay; log.az = az;
        log.gx = gx_f; log.gy = gy_f; log.gz = gz_f;
        log.baro_alt = current_baro_alt;
        log.est_alt = est_altitude;
        log.est_vel = est_velocity;
        log.vertical_acc = corrected_acc;
        log.event_status = current_status;

        dataFile.write((const uint8_t *)&log, sizeof(LogData));

        if (current_status == 4) 
        {
            dataFile.flush();
            dataFile.close();
            while(1) {
                digitalWrite(buzzerPin, HIGH); delay(500);
                digitalWrite(buzzerPin, LOW); delay(500);
            }
        }
    }

    if (dataFile && (millis() - lastFlush_ms >= 1000)) 
    {
        dataFile.flush();
        lastFlush_ms = millis();
    }
}