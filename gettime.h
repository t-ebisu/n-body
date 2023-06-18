#include <stdio.h>
#include <time.h>

double getTime(double *old_time) {
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    double now_time = ts.tv_nsec * 1e-9;
    double diff_time = now_time - *old_time;
    *old_time = now_time;

    return diff_time;
}