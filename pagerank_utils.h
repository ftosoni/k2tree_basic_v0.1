#include <stdio.h>
#include "kTree.h"
#include <assert.h>
//#include <time.h>

/*
void debug_double_vec(double* vec, size_t len, char* title){
    printf("%s\n", title);
    double sum = 0.0;
    for(size_t i=0; i<len; ++i){
        printf("%f\n", vec[i]);
        sum += vec[i];
    }
    printf("(sum: %f)\n", sum);
}

void debug_uint_vec(uint* vec, size_t len, char* title){
    printf("%s\n", title);
    for(size_t i=0; i<len; ++i){
        printf("%d\n", vec[i]);
    }
    printf("\n");
}

void debug_size_t_vec(size_t* vec, size_t len, char* title){
    printf("%s\n", title);
    for(size_t i=0; i<len; ++i){
        printf("%ld\n", vec[i]);
    }
    printf("\n");
}

void debug_bitmat(MREP *rep) {
    printf("nnodes: %d\n", rep->numberOfNodes);
    uint* arr = (uint*)calloc(rep->numberOfNodes, sizeof(uint));
    for(int n=0; n<rep->numberOfNodes; ++n) {
        if (n%(K*K) == 0)
            printf("\n");
        uint* adjl = compactAdjacencyList(rep, n);
        for(size_t j=0;j<rep->numberOfNodes;j++) {
            arr[j] = 0x0;
        }
        for(size_t j=0;j<adjl[0];j++){
            size_t e = adjl[j+1];
            arr[e] = 1;
        }
        for(size_t j=0;j<rep->numberOfNodes;j++){
            if(j%(K*K) == 0) {
                printf(" ");
            }
            printf("%d", arr[j]);
        }
        printf("\n");
    }
    free(arr);
}

void debug_mrep(MREP *rep) {
    const int lmax = 1 + rep->maxLevel;
    int l = 0;
    uint to_read = K*K, pc=0, bgn=0;
    while(l++ < lmax) {
        for (uint i=bgn; i < bgn + to_read; ++i) {
            if((i%(K*K))==0)
                printf(" ");
            uint bit = (0 != isBitSet(rep->btl, i));
            pc += bit;
            printf("%u", bit != 0);
        }
        bgn += to_read;
        to_read = pc * K * K;
        pc = 0;
        printf("\n");
    }
    assert(bitlen = end);
}*/

//const size_t NITERS = 8;
//const double ALPHA = 0.3;
//const unsigned TOPK = 3;

// heap based algorithm for finding the k largest ranks
// in heap order

// A utility function to swap two elements
static void swap(unsigned* a, unsigned* b) {
    unsigned t = *a;
    *a = *b;
    *b = t;
}

// Heapify a subtree rooted with node i which is an index in arr[].
// n is size of heap. the key associated to entry arr[i] is v[arr[i]]
static void minHeapify(double v[], unsigned arr[], unsigned n, unsigned i) {
    unsigned smallest = i;  // Initialize smallest as root
    unsigned left = 2*i + 1;
    unsigned right = 2*i + 2;

    // If left child is smaller than root
    if (left < n && v[arr[left]] < v[arr[smallest]])
        smallest = left;

    // If right child is smaller than smallest so far
    if (right < n && v[arr[right]] < v[arr[smallest]])
        smallest = right;

    // If smallest is not root
    if (smallest != i) {
        swap(&arr[i], &arr[smallest]);
        // Recursively heapify the affected sub-tree
        minHeapify(v, arr, n, smallest);
    }
}

// Function to find the k'th largest elements in an array
// v[0..n-1], arr[0..k-1] is the output array already allocated
static void kLargest(double v[], unsigned arr[], unsigned n, unsigned k) {
    assert(k<=n);
    assert(k>0);
    // init arr[] with the first k elements
    for(unsigned i=0;i<k;i++) arr[i] = i;
    // Build a min heap of first (k) elements in arr[]
    for (long int i = k / 2 - 1; i >= 0; i--)
        minHeapify(v, arr, k, i);
    // Iterate through the rest of the array elements
    for (unsigned i = k; i < n; i++) {
        // If current element is larger than root of the heap
        if (v[i] > v[arr[0]]) {
            // Replace root with current element
            arr[0] = i;
            // Heapify the root
            minHeapify(v, arr, k, 0);
        }
    }
}