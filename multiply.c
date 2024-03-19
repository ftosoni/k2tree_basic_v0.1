#include <stdio.h>
#include "kTree.h"
#include <assert.h>

void debug_double_vec(double* vec, size_t len){
    for(size_t i=0; i<len; ++i){
        printf("%f ", vec[i]);
    }
    printf("\n");
}

void debug_size_t_vec(size_t* vec, size_t len){
    for(size_t i=0; i<len; ++i){
        printf("%ld ", vec[i]);
    }
    printf("\n");
}

void debug_bitmat(MREP *rep) {
    printf("nnodes: %d\n", rep->numberOfNodes);
    uint* arr = (size_t*)calloc(rep->numberOfNodes, sizeof(uint));
    for(int n=0; n<rep->numberOfNodes; ++n) {
        if (n%(K*K) == 0)
            printf("\n");
        size_t curr = 0;
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
    const uint bitlen = lenght_in_bits(rep->btl);
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
}

void right_multiply_helper(
        MREP *rep,
        double* invec,
        double* outvec,
        size_t l, //current lvl
        size_t* ps, //per-level offset
        size_t row,
        size_t col,
        size_t dim //submatrix dimension at the current lvl
        ) {
    assert(pos < rep->btl_len());
    assert(row < rep->numberOfNodes);
    assert(col < rep->numberOfNodes);
    assert(dim > 0);
    size_t curr_row, curr_col;
    //base case: leaves
    if(l == rep->maxLevel) {
        assert(dim == 1);
        for(size_t i=0; i<K*K; ++i) {
            uint bit = (0 != isBitSet(rep->btl, ps[l]));
            curr_row = row + (i/K)*dim;
            curr_col = col + (i%K)*dim;

            //multiply
            outvec[curr_row] += bit * invec[curr_col];

            //debug
//            for(size_t i=0; i<l; ++i) {
//                printf("-");
//            }
//            printf("pos %ld @ lvl %ld     ", ps[l], l);
//            printf("%d @ (%ld, %ld)  +%f %f \n", bit, curr_row, curr_col, bit * invec[curr_col], outvec[curr_row]);

            //update
            ps[l] += 1;
        }
        return; //done
    }
    //inductive case: internal nodes
    for(size_t i=0; i<K*K; ++i) {
        uint bit = (0 != isBitSet(rep->btl, ps[l]));
        curr_row = row + (i/K)*dim;
        curr_col = col + (i%K)*dim;

//        printf("pos %ld @ lvl %ld\n", ps[l], l);

        //update
        ps[l] += 1;

        //recursion
        if (bit) {
            right_multiply_helper(rep, invec, outvec,
                                  l + 1, //current lvl
                                  ps, //per-level index of the first position
                                  curr_row, //row
                                  curr_col, //col,
                                  dim / K //submatrix dimension at the current lvl
            );
        }
    }
    return;
}

void right_multiply(MREP *rep, double* invec, double* outvec) {
    assert(len <= rep->numberOfNodes());

    //allocation, 0-initialised
    size_t* ps = (size_t*)calloc(1+rep->maxLevel, sizeof(size_t));
    if (ps == NULL) {
        printf("Errore nell'allocazione della memoria.\n");
        exit(2);
    }

    //init
    size_t to_read=K*K, pos=0;
    for(size_t l=0; l<rep->maxLevel; ++l) {
        ps[l+1] = ps[l] + to_read;
        size_t pc = 0;

        while (to_read--) {
            assert(pos < btl_len);
            uint bit = (0 != isBitSet(rep->btl, pos));
            pc += bit;

            //update
            pos++;
        }

        //update
        to_read = pc*K*K;
    }

//    printf("ps_bgn\n");
//    debug_size_t_vec(ps_bgn, 1+rep->maxLevel);
//    printf("rs_bgn\n");
//    debug_size_t_vec(rs_bgn, 1+rep->maxLevel);

    //logic
    size_t dim=1;
    for(size_t l=0; l<rep->maxLevel; ++l) {
        dim *= K;
    }
    assert(dim < rep->numberOfNodes);
    assert(K*rep->numberOfNodes <= dim);

    right_multiply_helper(rep, invec, outvec,
        0, //current lvl
        ps, //per-level offset
        0, //row
        0, //col,
        dim //submatrix dimension at the current lvl
    );

}

int main(int argc, char* argv[]) {

    if (argc != 2+1) {
        fprintf(stderr, "USAGE: %s <GRAPH> <invecpath>\n", argv[0]);
        exit(1);
    }

    MREP *rep = loadRepresentation(argv[1]);
//    debug_mrep(rep);
//    debug_bitmat(rep);
    double* invec = (size_t*)calloc(rep->numberOfNodes, sizeof(double));
    double* outvec = (size_t*)calloc(rep->numberOfNodes, sizeof(double));

    for(int i=0; i<14; ++i) {
        if(outvec[i] != 0x0){
            perror("Bad init.");
            exit(3);
        }
        invec[i] = i%10;
    }

    right_multiply(rep, invec, outvec);
    debug_double_vec(invec, 14);
    debug_double_vec(outvec, 14);

    destroyRepresentation(rep);
    printf("Done.");
    return 0;
}


