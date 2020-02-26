#ifndef DEBUGMACROS_H
#define DEBUGMACROS_H

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <omp.h>

#if defined(MYDEBUG)

/* Note: No separation between '"xxx"' and 'fmt'. */

#if defined(_OPENMP)
#define DBG(fmt, ...) fprintf(stderr, "### TH%02d %s: " fmt, omp_get_thread_num(), __FUNCTION__, ## __VA_ARGS__)
#define DBG0(fmt, ...) fprintf(stderr, fmt, ## __VA_ARGS__)
#else // ! _OPENMP
#define DBG(fmt, ...) fprintf(stderr, "### %s: " fmt, __FUNCTION__, ## __VA_ARGS__)
#define DBG0(fmt, ...) fprintf(stderr, fmt, ## __VA_ARGS__)
#endif // ! _OPENMP

#else // ! MYDEBUG

#define DBG(fmt, ...)
#define DBG0(fmt, ...)
//#define DBG(fmt, ...) fprintf(stderr, "### %s: " fmt, __FUNCTION__, ## __VA_ARGS__)
//#define DBG0(fmt, ...) fprintf(stderr, fmt, ## __VA_ARGS__)

#endif // ! MYDEBUG

#if defined(_DEBUG) || defined(MYDEBUG)
#if defined(_OPENMP) // 180619
#define ASSERT(c) {                             \
    if (!(c)) {                                 \
      fprintf(stderr, "# TH%02d %s: '%s' is violated at line %d. Exit.\n", omp_get_thread_num(), __FUNCTION__, #c, __LINE__); \
      exit(EXIT_FAILURE);                       \
    }                                           \
  }
#else
#include <assert.h>
#define ASSERT(s) assert(s)
#endif
#else
#define ASSERT(s)
#endif

#define ASSUME(c) {                             \
    if (!(c)) {                                 \
    INFO("'%s' is violated at line %d. Exit.\n", #c, __LINE__);		\
      exit(EXIT_FAILURE);                       \
    }                                           \
  }


#endif /* DEBUGMACROS_H */
