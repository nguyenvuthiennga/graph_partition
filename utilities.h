#ifndef __UTILITIES_H
#define __UTILITIES_H
#include <float.h>
#include <sys/time.h> 
#define DEF_MAX   DBL_MAX
#define DEF_MIN   (0.0 - DBL_MAX)


typedef struct timeval timeval_t;
void timeStart(timeval_t *tv); 
double timeEnd(timeval_t *tstart);

inline int getIndex(int row, int col, int n);

unsigned long int random_seed();

#endif  /* __UTILITIES_H */
