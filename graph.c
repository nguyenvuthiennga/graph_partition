/* 
 * Jan 2014
 * Author: Vu Thien Nga Nguyen
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include "utilities.h"
#include "graph.h"

/****** LOAD graphWORK/TARGET ******/
target_t *loadTarget(const char *file_name) {
  FILE *file = fopen(file_name, "r");
  assert(file != NULL && "Failed to open file \n");

  target_t *g = (target_t *) malloc(sizeof(target_t));

  int num_par;
  //first number is num_par
  fscanf(file, "%d", &num_par);
  assert (num_par > 0 && "Invalid number of node\n");
  g->num_par = num_par;

  // first $num_par values are partition capacity
  g->pars = (partition_t *) malloc(sizeof(partition_t) * num_par);

  int i, j;
  for (i = 0; i < num_par; i++) {
    g->pars[i].index = i;
    fscanf(file, "%lf", &g->pars[i].capacity);
    g->pars[i].weight = 0.0;
  }

  // next is the matrix of channel bandwidth
  g->channels = (double *) malloc(sizeof(double) * num_par * num_par);
  g->cuts = (double *) malloc(sizeof(double) * num_par * num_par);
  double val;
  for (i = 0; i < num_par; i++)
    for (j = 0; j < num_par; j++) {
      fscanf(file, "%lf", &val);
      int index = getIndex(i, j, num_par);
      g->channels[index] = val;
      g->cuts[index] = 0.0;
    }

  fclose(file);
  return g;
}

void delTarget(target_t **tptr) {
  target_t *target = *tptr;
  int i;

  free(target->pars);
  free(target->channels);
  free(target->cuts);
  free(target);
  tptr = NULL;
}

void delGraph(graph_t **nptr) {
  graph_t *graph = *nptr;
  int i;
  for (i = 0; i < graph->num_task; i++)
    free(graph->tasks[i].name);

  free(graph->tasks);
  free(graph->streams);
  free(graph);
  nptr = NULL;
}


graph_t *loadGraph(const char *file_name) {
  FILE *file = fopen(file_name, "r");
  assert(file != NULL && "Failed to open file \n");

  graph_t *g = (graph_t *) malloc(sizeof(graph_t));

  int num_task;
  //first number is num_task
  fscanf(file, "%d", &num_task);
  assert (num_task > 0 && "Invalid number of node\n");
  g->num_task = num_task;

  // first $num_task values are node weights
  g->tasks = (task_t *) malloc(sizeof(task_t) * num_task);

  int i, j;
  char buf[1000];
  for (i = 0; i < num_task; i++) {
    fscanf(file, "%s", buf);
    g->tasks[i].name = (char *)malloc(strlen(buf) + 1);
    strcpy(g->tasks[i].name, buf);
    g->tasks[i].index = i;
    fscanf(file, "%lf", &g->tasks[i].weight);
    g->tasks[i].partition = -1;
  }

  // next is the matrix of streams
  g->streams = (double *) malloc(sizeof(double) * num_task * num_task);
  double val;
  for (i = 0; i < num_task; i++)
    for (j = 0; j < num_task; j++) {
      fscanf(file, "%lf", &val);
      int index = getIndex(i, j, num_task);
      g->streams[getIndex(i, j, num_task)] = val;
    }

  fclose(file);
  return g;
}

void saveAlloc(graph_t *graph, const char *file_name) {
  FILE *f = fopen(file_name, "w");
  writeAlloc(graph, f);
  fclose(f);
}

void writeAlloc(graph_t *graph, FILE *f) {
  assert(f != NULL);
  int i;
  //number of 
  fprintf(f, "%d\n", graph->num_task);
  for (i = 0; i < graph->num_task; i++)
    fprintf(f, "%s %d\n", graph->tasks[i].name, graph->tasks[i].partition);
}


/****** COMPARIRSION FUNCTION TO SUPPORT LIST SORT *****/
int isTaskGreater(void *n1, void *n2) {
  return (((task_t *) n1)->weight - ((task_t*) n2)->weight);
}

int isTaskSmaller(void *n1, void *n2) {
  return (((task_t *) n1)->weight - ((task_t *) n2)->weight);
}



/******************** PARTITION UTILITIES *********/

double evaluatePar(target_t *target, congestion_t *bn) {
  double min = DEF_MAX;
  double val;
  int i, j;
  bn->par1 = -1;
  for (i = 0; i < target->num_par; i++) {
    if (target->pars[i].weight == 0.0)
      continue;
    val = target->pars[i].capacity / target->pars[i].weight;
    if (val < min) {
      min = val;
      bn->par2 = i;
    }
  }

  for (i = 0; i < target->num_par - 1; i++) {
    for (j = i + 1; j < target->num_par ; j++) {
      if (target->cuts[getIndex(i, j, target->num_par)] == 0.0)
        continue;
      val = target->channels[getIndex(i, j, target->num_par)]/target->cuts[getIndex(i, j, target->num_par)];
      if (val < min) {
        min = val;
        bn->par1 = i;
        bn->par2 = j;
      }
    }
  }

  return min;
}

/*
 * add a node to a partition
 *  - cut needs to be updated separatedly
 */
inline void addTaskToPar(partition_t *par, task_t *n) {
  n->partition = par->index;
  par->weight += n->weight;
}

inline void removeTaskFromPar(partition_t *par, task_t *n) {
  assert(n->partition == par->index);
  par->weight -= n->weight;
  n->partition = -1;    // not neccessary but should
}


void calParCut(graph_t *graph, target_t *target) {
  int i, j;
  int pi, pj;
  double eweight;

  for (pi = 0; pi < target->num_par; pi++)
    for (pj = 0; pj < target->num_par; pj++) {
      target->cuts[getIndex(pi, pj, target->num_par)] = 0.0;
      target->cuts[getIndex(pj, pi, target->num_par)] = 0.0;
    }

  for (i = 0; i < graph->num_task - 1; i++)
    for (j = i + 1; j < graph->num_task; j++) {
      eweight = graph->streams[getIndex(i, j, graph->num_task)];
      pi = graph->tasks[i].partition;
      pj = graph->tasks[j].partition;
      target->cuts[getIndex(pi, pj, target->num_par)] += eweight;
      /* if pi = pj, update only once */
      if (pi != pj) target->cuts[getIndex(pj, pi, target->num_par)] += eweight;
    }
}

/********* PARTITION **********/
void randomAlloc(const gsl_rng *r, graph_t *graph, target_t *target) {
  int par_id;
  int i;
  int num_par = target->num_par;
  for (i = 0; i < graph->num_task; i++) {
    par_id = gsl_rng_uniform_int(r, num_par);
    addTaskToPar(&target->pars[par_id], &graph->tasks[i]);
  }
  calParCut(graph, target);
}


/*  
 * exhaustive search for the best move
 *  for each node, test if it is the best to move to another partition
 * 
 * */
void searchBestKL(graph_t *graph, target_t *target, move_t *result) {
  assert(result != NULL);
  int i, j;
  task_t *node;
  partition_t *src;
  partition_t *dst;
  congestion_t bn;
  double val;
  result->value = DEF_MIN;

  for (i = 0; i < graph->num_task; i++) {
    if (graph->tasks[i].is_locked)
      continue;

    node = &graph->tasks[i];
    src = &target->pars[node->partition];
    for (j = 0; j < target->num_par; j++) {
      if (j == src->index)
        continue;

      dst = &target->pars[j];
      val = tryMove(graph, target, src, dst, node);
      if (val > result->value) {    /* it is ok to use absolute compare here */
        result->value = val;
        result->src = src;
        result->dst = dst;
        result->node = node;
      }
    }
  }
}


void savePar(graph_t *graph, int *alloc) {
  int i;
  for (i = 0; i < graph->num_task; i++)
    alloc[i] = graph->tasks[i].partition;
}

/* 
 * CA search: based on the congestion point
 */
int searchBestCA(graph_t *graph, target_t *target, move_t *result) {
  double val = DEF_MIN;
  congestion_t bn;
  partition_t *src1, *dst, *src2;
  task_t *cur_node1, *cur_node2;
  int i, j, k;
  int flag = 0;

  result->value = DEF_MIN;

  evaluatePar(target, &bn);
  if (bn.par1 == -1) {      // congestion lies on a partition's weight
    assert(bn.par2 >= 0 && bn.par2 < target->num_par);
    src2 = &target->pars[bn.par2];

    for (i = 0; i < graph->num_task; i++) {
      cur_node2 = &graph->tasks[i];
      if (cur_node2->partition != src2->index || cur_node2->is_locked)
        continue;

      for (j = 0; j < target->num_par; j++) {
        if (j == src2->index) continue;
        dst = &target->pars[j];
        val = tryMove(graph, target, src2, dst, cur_node2);
        if (val > result->value) {  // it is ok to use absolute compare here
          result->value = val;
          result->node = cur_node2;
          result->src = src2;
          result->dst = dst;
          flag = 1;
        }
      }
    } 

  } else {    /* congestion lies on communication channel */
    //    printf("comm\n");
    src1 = &target->pars[bn.par1];
    src2 = &target->pars[bn.par2];

    for (i = 0; i < graph->num_task; i++) {
      cur_node1 = &graph->tasks[i];
      if (cur_node1->partition != src1->index)
        continue;

      for (j = 0; j < graph->num_task; j++) {
        if(graph->tasks[j].partition == src2->index &&   /* connecting to the dst partition */
            graph->streams[getIndex(i, j, graph->num_task)] > 0.0)    /* having edge */
        {
          cur_node2 = &graph->tasks[j];

          // only consider to move cur_node if it is not locked 
          if (cur_node1->is_locked == 0) {
            // try to move cur_node first
            for (k = 0; k < target->num_par; k++) {
              if (k == src1->index)
                continue;
              dst = &target->pars[k];
              val = tryMove(graph, target, src1, dst , cur_node1);
              if (val > result->value) {  // it is ok to use absolute compare here 
                result->value = val;
                result->node = cur_node1;
                result->src = src1;
                result->dst = dst;
                flag = 1;
              }
            }
          }


          // only consider to move tmp_node if it is not locked 
          if (cur_node2->is_locked == 0) {
            for (k = 0; k < target->num_par; k++) {
              if (k == src2->index)
                continue;
              dst = &target->pars[k];
              val = tryMove(graph, target, src2, dst, cur_node2);
              if (val > result->value) {
                result->value = val;
                result->node = cur_node2;
                result->src = src2;
                result->dst = dst;
                flag = 1;
              }
            }
          }
        }
      }
    }
  }
  return flag;
}



/* 
 * congestion pass:
 *
 */
double CAPass(graph_t *graph, target_t *target) {
  double best_val = DEF_MIN;
  move_t moves[graph->num_task];   // maximum moves in 1 pass
  int i, peak_index = -1;

  for (i = 0; i < graph->num_task; i++)
    graph->tasks[i].is_locked = 0;
  for (i = 0; i < graph->num_task; i++) {
    if (searchBestCA(graph, target, &moves[i]) == 0)
      break;

    applyMove(graph, target, &moves[i]);
    moves[i].node->is_locked = 1;

    if (moves[i].value > best_val) {    // ok to use absolute compare
      best_val = moves[i].value;
      peak_index = i;
    }
  }

  int count = i;
  if (peak_index < 0) {   // no move found
    assert(count == 0);
    return -1.0;
  }
  // undo bad moves
  for (i = peak_index + 1; i < count; i++) {
    moveTask(graph, target, moves[i].dst, moves[i].src, moves[i].node);
  }

  return best_val;
}

double tryMove(graph_t *graph, target_t *target, partition_t *src, partition_t *dst, task_t *node) {
  moveTask(graph, target, src, dst, node);

  congestion_t bn;
  double val = evaluatePar(target, &bn);

  //undo
  moveTask(graph, target, dst, src, node);

  return val;
}

/* moveTask:
 *    move one node from src to dst
 */
void moveTask(graph_t *graph, target_t *target, partition_t *src, partition_t *dst, task_t *node) {
  int i;
  int connect;
  double eweight;


  assert(node->partition == src->index);
  assert(src->index != dst->index);
  node->partition = dst->index;
  /* update weight */
  src->weight -= node->weight;
  dst->weight += node->weight;

  /* update cut */
  for (i = 0; i < graph->num_task; i++) {
    if (i == node->index)
      continue;

    eweight = graph->streams[getIndex(i, node->index, graph->num_task)];
    if (eweight > 0.0){
      connect = graph->tasks[i].partition;
      if (connect == src->index) {   // now should increase the cut between src & dst
        target->cuts[getIndex(src->index, dst->index, target->num_par)] +=eweight;
        target->cuts[getIndex(dst->index, src->index, target->num_par)] +=eweight;

        //not neccessary but nice to have for future
        target->cuts[getIndex(src->index, src->index, target->num_par)] -=eweight;

      } else if (connect == dst->index) { // now should decrease the cut between src & dst
        target->cuts[getIndex(src->index, dst->index, target->num_par)] -=eweight;
        target->cuts[getIndex(dst->index, src->index, target->num_par)] -=eweight;

        //not neccessary but nice to have for future
        target->cuts[getIndex(dst->index, dst->index, target->num_par)] +=eweight;

      } else {  // decrease the cut between src & connect, increase the cut between connect & dst
        target->cuts[getIndex(src->index, connect, target->num_par)] -= eweight;
        target->cuts[getIndex(connect, src->index, target->num_par)] -= eweight;

        target->cuts[getIndex(connect, dst->index, target->num_par)] += eweight;
        target->cuts[getIndex(dst->index, connect, target->num_par)] += eweight;
      }
    }
  }
}


void applyMove(graph_t *graph, target_t *target, move_t *m) {
  moveTask(graph, target, m->src, m->dst, m->node);
}

/* 
 * KL pass:
 *    - search for the best move, lock the moved node, continue untill all
 *     have locked/moved
 *    - find the best value, apply moves to reach that value
 *    - in KL, no need to physically move
 */
double KLPass(graph_t *graph, target_t *target) {
  move_t m[graph->num_task];
  double best_val = DEF_MIN;
  int peak_index = -1;
  int i;
  // unlocked all 
  for (i = 0; i < graph->num_task; i++)
    graph->tasks[i].is_locked = 0;

  int count = 0;
  for (i = 0; i < graph->num_task; i++) {   
    /* every bestMove, one node is moved and locked --> apply this num_task time is a KL pass */
    searchBestKL(graph, target, &m[i]);
    if (m[i].value < 0)   // no move found
      break;

    applyMove(graph, target, &m[i]);
    // locked the node that is already moved
    m[i].node->is_locked = 1;

    if (m[i].value > best_val) {    // ok to use absolute compare
      best_val = m[i].value;
      peak_index = i;
    }
    count++;
  }

  if (peak_index < 0)
    return -1.0;
  /* undo bad moves */
  for (i = peak_index + 1; i < count; i++)
    moveTask(graph, target, m[i].dst, m[i].src, m[i].node);

  return best_val;
}

/* 
 * compare 2 double values within the tolerance of EPS
 * new_val is better only if it is at least EPS greater than old_val 
 */ 
int isBetter(double new_val, double old_val) {
  if ((new_val - old_val)/old_val > EPS)
    return 1;

  return 0;
}


void refine(graph_t *graph, target_t *target, search_t method, double *best_val) {
  int best_alloc[graph->num_task];
  congestion_t bn;
  *best_val = evaluatePar(target, &bn);
  savePar(graph, best_alloc);

  timeval_t tv;
  double t = 0;
  double cur_val;
  double (*pass_func)(graph_t *graph, target_t *target);

  if (method == KL)
    pass_func = &KLPass;
  else  // CA
    pass_func = &CAPass;

  int stop = 0;
  do {
    cur_val = pass_func(graph, target);
    if (isBetter(cur_val, *best_val)) {     /* better pass */
      /* theory: cur_val > *best_val
         but isBetter is used to avoid infinite loop due to double comparision */

      *best_val = cur_val;
      //save the best
      savePar(graph, best_alloc);
    } else
      stop = 1;
  } while (!stop);

  *best_val = applyAlloc(graph, target, best_alloc);
}



void partition(graph_t *graph, target_t *target, search_t method, double *best_val, const gsl_rng *r) {

  randomAlloc(r, graph, target);
  refine(graph, target, method, best_val);

}

void resetPar(graph_t *graph, target_t *target) {
  int i;
  for (i = 0; i < target->num_par; i++)
    target->pars[i].weight = 0.0;

  for (i = 0; i < graph->num_task; i++)
    graph->tasks[i].partition = -1;
}

double applyAlloc(graph_t *graph, target_t *target, int *alloc) {
  int i;
  resetPar(graph, target);

  for (i = 0; i < graph->num_task; i++)
    addTaskToPar(&target->pars[alloc[i]], &graph->tasks[i]);

  calParCut(graph, target);
  congestion_t bn;
  return evaluatePar(target, &bn);
}







