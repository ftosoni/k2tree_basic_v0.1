#ifndef PAGERANK_UTILS_H
#define PAGERANK_UTILS_H

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <pthread.h>
#include <stdio.h>
#include <assert.h>
#include <pthread.h>
#include <sched.h>

#include "kTree.h"

extern int set_core(pthread_t *thread, int tid, const int ncores);

extern int my_pthread_create(pthread_t *thread, const pthread_attr_t *attr, void *(*start_routine)(void *), void *arg, int tid, const int ncores);

extern void swap(unsigned* a, unsigned* b);

extern void minHeapify(double v[], unsigned arr[], unsigned n, unsigned i);

extern void kLargest(double v[], unsigned arr[], unsigned n, unsigned k);

#endif // PAGERANK_UTILS_H