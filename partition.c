/* 
 * Jan 2014
 * Author: Vu Thien Nga Nguyen
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <sys/time.h>
#include <assert.h>
#include <math.h>
#include "utilities.h"
#include "graph.h"

void printUsage() {
  printf("usage: -gf <graph_file> -tf <target_file> -m <method> -of <output_file>\n");
  printf("      <method>: 1 = CA, 2 = KL, default = CA\n");
}

int main(int argc, char **argv) {
  char *graph_fn = NULL;
  char *tar_fn = NULL;
  char *out_fn = NULL;
  int m = 1;
  int i;
  for (i = 0; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      printUsage();
      exit(1);
    }

    if (strcmp(argv[i], "-gf") == 0)
      graph_fn = argv[++i];
    else if (strcmp(argv[i], "-tf") == 0)
      tar_fn = argv[++i];
    else if (strcmp(argv[i], "-m") == 0)
      m = atoi(argv[++i]);
    else if (strcmp(argv[i], "-of") == 0)
      out_fn = argv[++i];
  }

  if (graph_fn == NULL || tar_fn == NULL || out_fn == NULL) {
    printUsage();
    exit(1);
  }


  const gsl_rng * r = gsl_rng_alloc (gsl_rng_env_setup()) ;

  //  gsl_ieee_env_setup ();
  gsl_rng_set (r, random_seed());  

  search_t method = CA;
  if (m == 2)
    method = KL;

  graph_t *graph = loadGraph(graph_fn);
  target_t *target = loadTarget(tar_fn);

  double val;
  timeval_t tv;
  double t;

  timeStart(&tv);
  partition(graph, target, method, &val, r);
  t = timeEnd(&tv);
  saveAlloc(graph, out_fn);

  printf("%s: best tp = %lf, et = %lf\n", method == CA ? "CA" : "KL", val, t);
  delGraph(&graph);
  delTarget(&target);

  return 0;
}







