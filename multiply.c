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

void fill_double_vec(char * infilepath, double * tofill, const size_t len) {
    FILE *file;
    size_t fileSize, numDoubles;

    // Apertura del file in modalità binaria
    file = fopen(infilepath, "rb");

    if (file == NULL) {
        perror("Error while opening file");
        exit(1);
    }

    // Calcolo della dimensione del file
    fseek(file, 0, SEEK_END);
    fileSize = ftell(file);
    fseek(file, 0, SEEK_SET);

    // Calcolo del numero di double nel file
    numDoubles = fileSize / sizeof(double);
    assert(numDoubles == len);

    // Lettura dei dati dal file
    fread(tofill, sizeof(double), numDoubles, file);

    // Liberazione della memoria e chiusura del file
    fclose(file);

    return;
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
        printf("Memory allocation error.\n");
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

    if (argc != 3+1) {
        fprintf(stderr, "USAGE: %s <GRAPH> <invecpath> <outvecpath>\n", argv[0]);
        exit(-1);
    }

    MREP *rep = loadRepresentation(argv[1]);
//    debug_mrep(rep);
//    debug_bitmat(rep);

    double* invec = (double*)calloc(rep->numberOfNodes, sizeof(double));
    double* outvec = (double*)calloc(rep->numberOfNodes, sizeof(double));
    fill_double_vec(argv[2], invec, rep->numberOfNodes);
//    for(int i=0; i<rep->numberOfNodes; ++i) {
//        if(outvec[i] != 0x0){
//            perror("Bad initialisation for output vector.");
//            exit(3);
//        }
//        invec[i] = i%10;
//    }

    //invec
    FILE *in_invec = fopen(argv[2], "r");
    if(in_invec == NULL){
        perror("Cannot open input vector file");
        exit(2);
    }
    fseek(in_invec, 0, SEEK_END);
    size_t len_invec = ftell(in_invec);
    fseek(in_invec, 0, SEEK_SET);
    fread(&invec[0], 1, len_invec, in_invec);

    //multiplication
    right_multiply(rep, invec, outvec);
//    debug_double_vec(invec, rep->numberOfNodes);
//    debug_double_vec(outvec, rep->numberOfNodes);

    free(invec);

    //outfile
    FILE *out_outvec = fopen(argv[3], "wb");  // Open in binary format
    if(out_outvec == NULL) {
        perror("Unable to open output vector file");
        exit(3);
    }
    fwrite(&outvec[0], sizeof(double), rep->numberOfNodes, out_outvec);
    fclose(out_outvec);  // Close the file

    destroyRepresentation(rep);
    free(outvec);
    return 0;
}


