/* 
 * Jan 2014
 * Author: Vu Thien Nga Nguyen
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "utilities.h"

void timeStart(timeval_t *tv) {
  gettimeofday(tv, NULL);
}

double timeEnd(timeval_t *tstart) {
  timeval_t tend;
  gettimeofday(&tend, NULL);
  double x = ((double)tend.tv_sec + 1.0e-6*tend.tv_usec) -
    ((double)tstart->tv_sec + 1.0e-6*tstart->tv_usec);
  return x;
}


inline int getIndex(int row, int col, int n) {
  assert(row < n && col < n);
  return row * n + col;
}

unsigned long int random_seed()
{

  unsigned int seed;
  struct timeval tv; 
  FILE *devrandom;

  if ((devrandom = fopen("/dev/random","r")) == NULL) {
    gettimeofday(&tv,0);
    seed = tv.tv_sec + tv.tv_usec;
  } else {
    fread(&seed,sizeof(seed),1,devrandom);
    fclose(devrandom);
  }
  return(seed);
}

