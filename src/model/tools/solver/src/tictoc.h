#ifndef TICTOC_H
#define TICTOC_H

#include <iostream>
#include <chrono>

typedef std::chrono::high_resolution_clock Clock;
typedef std::chrono::milliseconds milliseconds;
typedef std::chrono::microseconds microseconds;
typedef std::chrono::nanoseconds nanoseconds;
static Clock::time_point t0 = Clock::now();

using namespace std::chrono;
using namespace std;

void tic()
{
 t0 = Clock::now();
 //cout.sync_with_stdio(false);
}

void toc(float tsim = 1.0f)
{
    Clock::time_point t1 = Clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t1 - t0);

    cout <<"Elapsed time:  "<< time_span.count() << " s\n";
    cout <<"Sampling rate: "<< tsim/(time_span.count()) << " Hz\n";
}

#endif