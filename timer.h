#ifndef TYMER_H
#define TYMER_H

#include <stdio.h>
#include <stdlib.h>
#include <chrono>

typedef struct {
  double time;
  std::chrono::system_clock::time_point start;
  std::chrono::system_clock::time_point end;
} Timer;

void timer_malloc(Timer **timer);
void timer_init(Timer *timer);
void timer_start(Timer *timer);
void timer_stop(Timer *timer);
void timer_free(Timer *timer);
double timer_get(Timer *timer);
void timer_print(FILE *fp, const char *s, Timer *timer);
void timer_create(Timer **timer);

#define TCREATE(t) Timer *t; timer_create(&t)
#define TSTART(t) timer_start(t)
#define TSTOP(t) timer_stop(t)
#define TPRINT(t, outfile) outfile << #t << ": " << timer_get(t) << std::endl
#define TFREE(t) timer_free(t)

#define TICK(t) TCREATE(t); TSTART(t)
#define TACK(t, outfile) TSTOP(t); TPRINT(t, outfile); TFREE(t)

#endif // TIMER_H
