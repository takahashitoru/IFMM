#ifndef MATHMACROS_H
#define MATHMACROS_H

//150706#ifndef PI
//150706#define PI (3.1415926535897932384626433832795029)
//150706#else
//150706#error Redefined.
//150706#endif

//150706#ifndef PI2
//150706#define PI2 (6.2831853071795864769252867665590058)
//150706#else
//150706#error Redefined.
//150706#endif

//150708#ifndef PI4
//150708#define PI4 (12.5663706143591729538505735331180116)
//150708#else
//150708#error Redefined.
//150708#endif

//150708#ifndef PI4I
//150708#define PI4I (7.957747154594767e-02)
//150708#else
//150708#error Redefined.
//150708#endif

#include <math.h>


#ifndef FLOORUP
#define FLOORUP(n, m) (((n) & ((m) - 1)) ? ((((n) / (m)) + 1) * (m)) : (n))
#else
#error Redefined.
#endif

#ifndef MAX
#define MAX(x, y) (( (x) > (y) ) ? (x) : (y))
#else
#error Redefined.
#endif

#ifndef MIN
#define MIN(x, y) (( (x) < (y) ) ? (x) : (y))
#else
#error Redefined.
#endif

#ifndef SQUARE
#define SQUARE(n) ((n) * (n))     /* n^2 */
#else
#error Redefined.
#endif

#ifndef CUBE
#define CUBE(n) ((n) * (n) * (n)) /* n^3 */
#else
#error Redefined.
#endif

#ifndef QUAD
#define QUAD(n) ((n) * (n) * (n) * (n)) /* n^4 */
#else
#error Redefined.
#endif

#ifndef POW2
#define POW2(i) (1 << (i))        /* 2^i */
#else
#error Redefined.
#endif

#ifndef POW4
#define POW4(i) (1 << ((i) << 1)) /* 4^i */
#else
#error Redefined.
#endif

#ifndef POW8
#define POW8(i) (1 << (3 * (i)))  /* 8^i */
#else
#error Redefined.
#endif

#ifndef ROUND
#define ROUND(x) ( (int)round(x) )
#else
#error Redefined.
#endif


#endif /* MATHMACROS_H */
