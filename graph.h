/* 
 * Jan 2014
 * Author: Vu Thien Nga Nguyen
 *
 */


#ifndef __GRAPH_H
#define __GRAPH_H
#include <gsl/gsl_rng.h>
#include "utilities.h"

#define RND_INIT   0
#define BAL_INIT   1
#define EPS       1.1920928955078125e-07


typedef enum {KL, CA} search_t;

typedef struct task_t task_t;
struct task_t {
  char *name;       /* unique name */
  int index;
  double weight;    /* time unit contributes to process 1 external messages */
  int partition;    /* partition id that the task is assigned to */
  int is_locked;    /* used for moving during refinement */
};



typedef struct {
  int index;

  /* hardware properties */
  double capacity;    /* in this case: number of cores/workers */

  /* partition properties */
  double weight;       /* total weight (execution time) of tasks that assigned to the partition */
} partition_t;



/*
 * Graph Structure
 *  Task index is used as task id, and is also used as row/column index of tasks & streams
 *  order of tasks is not changed
 */
typedef struct {
  int num_task;
  task_t *tasks;  /* array of tasks */
  double *streams;    /* upper triangular matrix, streams[i][j] = weight of stream S_ij (connecting task i and task j) */
} graph_t;

/*
 * target structure
 * each machine in the target is represented as a partition
 * partition index is used as partition id, and is also used as row/colum index of cuts and channels
 */
typedef struct {
  /* hardware properties */
  int num_par;
  partition_t *pars;
  double *channels;    /* in this case, it is the bandwidth = amount of data/time unit */

  /* partition properties */
  double *cuts;  /* upper triangular matrix, cuts[i][j] = cut between partition i and partition j */
} target_t;


/* congestion point structure
 * if par1 == -1 --> congestion point lies in par2
 * else bottleneck lies in communication between par1 & par2 (par1 <= par2)
 */
typedef struct {
  int par1;
  int par2;
} congestion_t;


/* move structure */
typedef struct move_t move_t;
struct move_t{
  partition_t *src;
  partition_t *dst;
  task_t *node;
  double value; 
};



/* load task graph and target from file */
target_t *loadTarget(const char *file_name);
graph_t *loadGraph(const char *file_name);

void delGraph(graph_t **n);
void delTarget(target_t **target);
void saveAlloc(graph_t * graph, const char *file_name);
void writeAlloc(graph_t *graph, FILE *f);


/* utility functions */
void undoMove(void *m);
double tryMove(graph_t *graph, target_t *target, partition_t *src, partition_t *dst, task_t *node);
void applyMove(graph_t *graph, target_t *target, move_t *m);
void moveTask(graph_t *graph, target_t *target, partition_t *src, partition_t *dst, task_t *node);

double evaluatePar(target_t *target, congestion_t *bn);
inline void addTaskToPar(partition_t *par, task_t *n);
inline void removeTaskFromPar(partition_t *par, task_t *n);
void calParCut(graph_t *graph, target_t *target);
int isBetter(double, double);


/* Init partition */
void randomAlloc(const gsl_rng *r, graph_t *graph, target_t *target);

/* clear all allocations */
void resetPar(graph_t *graph, target_t *target);

/* refinement */
void searchBestKL(graph_t *graph, target_t *target, move_t *result);
int searchBestCA(graph_t *graph, target_t *target, move_t *result);
double KLPass(graph_t *graph, target_t *target);
double CAPass(graph_t *graph, target_t *target);
void refine(graph_t *graph, target_t *target, search_t method, double *best_val);


void partition(graph_t *graph, target_t *target, search_t method, double *best_val, const gsl_rng *r);

/* apply allocation
 * alloc: array of parition id of each task, alloc[i] = partittion id that task i is allocated to
 */
double applyAlloc(graph_t *graph, target_t *target, int *alloc);


/* compare two tasks */
int isTaskGreater(void *n1, void *n2);
int isTaskSmaller(void *n1, void *n2);


#endif    /* __GRAPH_H */

