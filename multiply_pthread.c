#include <stdio.h>
#include "kTree.h"
#include <assert.h>
#include <pthread.h>

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
        perror(infilepath);
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

typedef struct {
    uint tid;
    MREP **reps;
    char* argv1;
    char* argv4;
} targs_read_t;

void* pthread_init(void* args) {
    //args
    targs_read_t* targs = (targs_read_t*)args;
    uint tid = targs->tid;
    char* argv1 = targs->argv1;
    char* argv4 = targs->argv4;
    MREP** reps = targs->reps;

    //construct the infilepath
    uint len_argv1 = strlen(argv1);
    char* infilepath = calloc(len_argv1 + 4 + 1, sizeof(char));
    uint pos = 0;
    strcpy(infilepath + pos, argv1); pos += len_argv1;
    strcpy(infilepath + pos, "."); pos += 1;
    strcpy(infilepath + pos, argv4); pos += 1;
    strcpy(infilepath + pos, "."); pos += 1;
    sprintf(infilepath + pos, "%d", tid); pos+= 1;

    printf("%d: %s\n", tid, infilepath);
    reps[tid] = loadRepresentation(infilepath);
    free(infilepath);
    return NULL;
}

typedef struct {
    uint tid;
    MREP **reps;
    double* invec;
    double* outvec;
} targs_mult_t;

void* pthread_mult(void* args) {
    targs_mult_t* targs = (targs_mult_t*)args;
    uint tid = targs->tid;
    MREP** reps = targs->reps;
    double* invec = targs->invec;
    double* outvec = targs->outvec;

    right_multiply(reps[tid], invec, outvec);
    return NULL;
}

int main(int argc, char* argv[]) {

    if (argc != 4+1) {
        fprintf(stderr, "USAGE: %s <GRAPH> <invecpath> <outvecpath> <par. degree>\n", argv[0]);
        exit(-1);
    }

    //args
    const uint NT = atoi(argv[4]);
    assert(NT < 10); //todo case with >9 threads

    //read
    MREP **reps = (MREP **) calloc(NT, sizeof(MREP *));
    pthread_t *threads = (pthread_t *) calloc(NT, sizeof(pthread_t));

    targs_read_t *targs_read = (targs_read_t *) calloc(NT, sizeof(targs_read_t));
    for (uint tid = 0; tid < NT; ++tid) {
        targs_read[tid].tid = tid;
        targs_read[tid].argv1 = argv[1];
        targs_read[tid].argv4 = argv[4];
        targs_read[tid].reps = reps;
        pthread_create(&threads[tid], NULL, pthread_init, &targs_read[tid]);
    }
    for (uint tid = 0; tid < NT; ++tid) {
        pthread_join(threads[tid], NULL);
    }
    free(targs_read);

    const uint nnodes = reps[0]->numberOfNodes;
//    debug_mrep(rep);
//    debug_bitmat(rep);

    double* invec  = (double*)calloc(nnodes, sizeof(double));
    double* outvec = (double*)calloc(nnodes, sizeof(double));
    fill_double_vec(argv[2], invec, nnodes);

    //invec
    FILE *in_invec = fopen(argv[2], "r");
    if(in_invec == NULL){
        perror(argv[2]);
        exit(2);
    }
    fseek(in_invec, 0, SEEK_END);
    size_t len_invec = ftell(in_invec);
    fseek(in_invec, 0, SEEK_SET);
    fread(&invec[0], 1, len_invec, in_invec);

    //multiplication
    targs_mult_t *targs_mult = (targs_mult_t *) calloc(NT, sizeof(targs_mult_t));
    for (uint tid = 0; tid < NT; ++tid) {
        targs_mult[tid].tid = tid;
        targs_mult[tid].outvec = outvec;
        targs_mult[tid].invec = invec;
        targs_mult[tid].reps = reps;
        pthread_create(&threads[tid], NULL, pthread_mult, &targs_mult[tid]);
    }
    for (uint tid = 0; tid < NT; ++tid) {
        pthread_join(threads[tid], NULL);
    }
    free(targs_mult);
//    debug_double_vec(invec, rep->numberOfNodes);
//    debug_double_vec(outvec, rep->numberOfNodes);

    free(invec);

    //outfile
    FILE *out_outvec = fopen(argv[3], "wb");  // Open in binary format
    if(out_outvec == NULL) {
        perror(argv[3]);
        exit(3);
    }
    fwrite(&outvec[0], sizeof(double), nnodes, out_outvec);
    fclose(out_outvec);  // Close the file

    for(uint tid=0; tid<NT; ++tid) {
        destroyRepresentation(reps[tid]);
    }
    free(reps);
    free(outvec);
    return 0;
}


