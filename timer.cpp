#include "timer.h"

/*
  See 
  https://en.cppreference.com/w/cpp/chrono/time_point
  http://vivi.dyndns.org/tech/cpp/timeMeasurement.html
*/


void timer_malloc(Timer **timer)
{
  *timer = (Timer *)malloc(sizeof(Timer));
  if (*timer == NULL) {
    fprintf(stderr, "### %s: fail to allocate. Exit.\n", __FUNCTION__);
    exit(EXIT_FAILURE);
  }
}

void timer_init(Timer *timer)
{
  timer->time = 0.0;
}  

void timer_start(Timer *timer)
{
  timer->start = std::chrono::system_clock::now();
}

void timer_stop(Timer *timer)
{
  timer->end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = timer->end - timer->start;
  timer->time += elapsed_seconds.count();
}

void timer_free(Timer *timer)
{
  free(timer);
}

double timer_get(Timer *timer)
{
  return timer->time;
}

void timer_print(FILE *fp, const char *s, Timer *timer)
{
  fprintf(fp, "# %s: %s = %f\n", __FUNCTION__, s, timer->time);
}

void timer_create(Timer **timer)
{
  timer_malloc(timer);
  timer_init(*timer);
}
