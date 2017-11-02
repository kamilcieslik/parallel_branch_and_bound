//
// Created by mrfarinq on 22.05.17.
//

#ifndef SDIZO_2_TIMEMEASUREMENT_H
#define SDIZO_2_TIMEMEASUREMENT_H

#include <iostream>
#include <chrono>

typedef std::chrono::high_resolution_clock Clock;


class TimeMeasurement {
private:
    Clock::time_point timeStart;
    Clock::time_point timeStop;

public:
    TimeMeasurement();
    
    ~TimeMeasurement();
    
    void TimeStart();
    
    void TimeStop();
    
    double GetTimeInSeconds();
    
    const std::string currentDateTime();
};


#endif //SDIZO_2_TIMEMEASUREMENT_H
