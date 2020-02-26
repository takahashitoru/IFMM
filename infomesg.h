#ifndef INFOMESG_H
#define INFOMESG_H

#include <stdio.h>
#include <stdlib.h>
#include <omp.h> // 180619

#ifndef INFO
#if defined(_OPENMP)
#define INFO(fmt, ...) fprintf(stderr, "# TH%02d %s: " fmt, omp_get_thread_num(), __FUNCTION__, ## __VA_ARGS__)
#else
#define INFO(fmt, ...) fprintf(stderr, "# %s: " fmt, __FUNCTION__, ## __VA_ARGS__)
#endif
#else
#error Already defined.
#endif

#endif /* INFOMESG_H */
