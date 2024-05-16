#include <pthread.h>
#include "pagerank_utils.h"

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
}

//end debug

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
    if(fread(tofill, sizeof(double), numDoubles, file) != numDoubles){
        perror("Error while opening file");
        exit(1);
    };

    // Liberazione della memoria e chiusura del file
    fclose(file);

    return;
}

void right_multiply_helper(
        const MREP *rep,
        const double* invec,
        double* outvec,
        const size_t l, //current lvl
        size_t* ps, //per-level offset
        const size_t row,
        const size_t col,
        const size_t dim //submatrix dimension at the current lvl
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

//            printf("%d @ (%d, %d)\n", bit, curr_row, curr_col);

            //multiply
            outvec[curr_row] += bit * invec[curr_col];

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

void compute_outdeg_helper(
        const MREP *rep,
        uint* outdeg,
        const size_t l, //current lvl
        size_t* ps, //per-level offset
        const size_t row,
        const size_t col,
        const size_t dim //submatrix dimension at the current lvl
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

//            printf("%d @ (%d, %d)\n", bit, curr_row, curr_col);

            //multiply
            outdeg[curr_col] += bit;

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
            compute_outdeg_helper(rep, outdeg,
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

struct fetch_data_t {
    char *argv1;
    char *argv4;
    MREP **reps;
    uint tid;
};

void* fetch_f(void *arg){
    struct fetch_data_t* data = (struct fetch_data_t*)arg;

    //construct the infilepath
    uint len_argv1 = strlen(data->argv1);
    char* infilepath = calloc(len_argv1 + 4 + 1, sizeof(char));
    uint pos = 0;
    strcpy(infilepath + pos, data->argv1); pos += len_argv1;
    strcpy(infilepath + pos, "."); pos += 1;
    strcpy(infilepath + pos, data->argv4); pos += 1;
    strcpy(infilepath + pos, "."); pos += 1;
    sprintf(infilepath + pos, "%d", data->tid); pos+= 1;

//    printf("%d: %s\n", data->tid, infilepath);
    data->reps[data->tid] = loadRepresentation(infilepath);
    free(infilepath);
    return NULL;
}

struct xsum_data_t {
    double* xsum_helper;
    double* outvec;
    size_t nnodes;
    uint NT;
    uint tid;
};

void* xsum_f(void* arg) {
    struct xsum_data_t* data = (struct xsum_data_t*)arg;
    data->xsum_helper[data->tid] = 0.0;
    size_t col_block_size = (data->nnodes + data->NT - 1) / data->NT;
    size_t col_bgn = data->tid * col_block_size;
    size_t col_end = col_bgn + col_block_size;
    if (col_end > data->nnodes) {
        col_end = data->nnodes;
    }
    for (size_t c = col_bgn; c < col_end; c++) {
        data->xsum_helper[data->tid] += data->outvec[c];
    }
    return NULL;
}

struct xnormalise_data_t {
    double xsum;
    double* outvec;
    size_t nnodes;
    uint NT;
    uint tid;
};

void* xnormalise_f(void* arg) {
    struct xnormalise_data_t* data = (struct xnormalise_data_t*)arg;
    size_t col_block_size = (data->nnodes + data->NT - 1) / data->NT;
    size_t col_bgn = data->tid * col_block_size;
    size_t col_end = col_bgn + col_block_size;
    if (col_end > data->nnodes) {
        col_end = data->nnodes;
    }
    for (size_t c = col_bgn; c < col_end; c++) {
        data->outvec[c] /= data->xsum;
    }
    return NULL;
}

struct compute_ps_data_t {
    MREP *rep;
    size_t* ps;
};

void* compute_ps_f(void* arg) {
    struct compute_ps_data_t *data = (struct compute_ps_data_t *) arg;
    //assuming ps is 0-init
    assert(len <= rep->numberOfNodes());

    //init
    size_t to_read = K * K, pos = 0;
    for (size_t l = 0; l < data->rep->maxLevel; ++l) {
        data->ps[l + 1] = data->ps[l] + to_read;
        size_t pc = 0;

        while (to_read--) {
            assert(pos < btl_len);
            uint bit = (0 != isBitSet(data->rep->btl, pos));
            pc += bit;

            //update
            pos++;
        }

        //update
        to_read = pc * K * K;
    }

    //    printf("ps_bgn\n");
    //    debug_size_t_vec(ps_bgn, 1+rep->maxLevel);
    //    printf("rs_bgn\n");
    //    debug_size_t_vec(rs_bgn, 1+rep->maxLevel);
    return NULL;
}

struct compute_outdeg_data_t {
    MREP *rep;
    uint* outdeg;
    size_t* ps;
    size_t* ps_orig;
    size_t dim;
};

void* compute_outdeg_f(void *arg) {
    struct compute_outdeg_data_t *data = (struct compute_outdeg_data_t *) arg;
    //copy ps
    const size_t ps_len = 1+data->rep->maxLevel;
    memcpy(data->ps, data->ps_orig, ps_len * sizeof(size_t));
    //business logic
    compute_outdeg_helper(data->rep, data->outdeg,
                          0, //current lvl
                          data->ps, //per-level offset
                          0, //row
                          0, //col,
                          data->dim //submatrix dimension at the current lvl
    );
    return NULL;
}

//loop from here

struct loop_data_t {
    double *contrib_dn_helper;
    size_t nnodes;
    uint *outdeg;
    double *invec;
    double *outvec;
    size_t ps_len;
    size_t* ps;
    size_t* ps_orig;
    MREP *rep;
    size_t dim;
    uint NT;
    uint tid;
};

void* contrib_dn_f (void *arg) {
    struct loop_data_t *data = (struct loop_data_t *) arg;

    //swap invec & outvec
    {
        double *tmp = data->outvec;
        data->outvec = data->invec;
        data->invec = tmp;
    }

    data->contrib_dn_helper[data->tid] = 0.0;
    const size_t x_block_size = (data->nnodes + data->NT - 1) / data->NT;
    const size_t x_bgn = data->tid * x_block_size;
    size_t x_end = x_bgn + x_block_size;
    if (x_end > data->nnodes) {
        x_end = data->nnodes;
    }
    for (size_t x = x_bgn; x < x_end; x++) {
        data->contrib_dn_helper[data->tid] += (data->outdeg[x] == 0) * data->invec[x]; //contribution of dangling nodes
        if (data->outdeg[x]) data->invec[x] /= data->outdeg[x]; //divide input vector by outdegree
        data->outvec[x] = 0.0; //0-init
    }
    data->contrib_dn_helper[data->tid] /= data->nnodes;
    //copy ps
    memcpy(data->ps, data->ps_orig, data->ps_len * sizeof(size_t));
    return NULL;
}

void* right_multiply_f(void *arg) {
    struct loop_data_t *data = (struct loop_data_t*) arg;
    right_multiply_helper(data->rep, data->invec, data->outvec,
                          0, //current lvl
                          data->ps, //per-level offset
                          0, //row
                          0, //col,
                          data->dim //submatrix dimension at the current lvl
    );
    return NULL;
}

void* finalise_f(void* arg) {
    struct loop_data_t* data = (struct loop_data_t*)arg;
    size_t row_block_size = (data->nnodes + data->NT - 1) / data->NT;
    size_t row_bgn = data->tid * row_block_size;
    size_t row_end = row_bgn + row_block_size;
    if (row_end > data->nnodes) {
        row_end = data->nnodes;
    }
    for (size_t r = row_bgn; r < row_end; r++) {
        data->outvec[r] += data->contrib_dn_helper[0]; //dangling nodes
        data->outvec[r] = (1 - ALPHA) * data->outvec[r] + ALPHA / data->rep->numberOfNodes; //teleporting
    }
    return NULL;
}

int main(int argc, char* argv[]) {

    if (argc != 3+1) {
        fprintf(stderr, "USAGE: %s <GRAPH> <count col file> <par. degree>\n", argv[0]);
        exit(-1);
    }

    //args
    const uint NT = atoi(argv[3]);
    assert(NT < 10); //todo case with >9 threads

    //data
    pthread_t *threads = (pthread_t *) calloc(NT, sizeof(pthread_t));
    MREP **reps = (MREP **) calloc(NT, sizeof(MREP *));
    assert(reps != NULL);

    //read matrices
    struct fetch_data_t *fetch_data = (struct fetch_data_t *) calloc(NT, sizeof(struct fetch_data_t));
    for (uint tid = 0; tid < NT; ++tid) {
        fetch_data[tid].argv1 = argv[1];
        fetch_data[tid].argv4 = argv[3];
        fetch_data[tid].reps = reps;
        fetch_data[tid].tid = tid;
        pthread_create(&threads[tid], NULL, fetch_f, &fetch_data[tid]);
    }
    for (uint tid = 0; tid < NT; ++tid) {
        pthread_join(threads[tid], NULL);
        assert(reps[0]->numberOfNodes == reps[tid]->numberOfNodes);
    }
    free(fetch_data);

//    debug_mrep(rep);
//    debug_bitmat(rep);

    u_int32_t nnodes = reps[0]->numberOfNodes;  // Set the value of nnodes

    double* invec = (double*)calloc(nnodes, sizeof(double));
    double* outvec = (double*)calloc(nnodes, sizeof(double));

    //initialise input vector
    {
        for(size_t c=0; c<nnodes; ++c) {
            outvec[c] = 1.0/nnodes;
        }
    }

//    //sum
//    double* xsum_helper = (double*)calloc(NT, sizeof(double));
//    struct xsum_data_t* xsum_data = (struct xsum_data_t*) calloc(NT, sizeof(struct xsum_data_t));
//    for (uint tid = 0; tid < NT; ++tid) {
//        xsum_data[tid].xsum_helper = xsum_helper;
//        xsum_data[tid].outvec = outvec;
//        xsum_data[tid].nnodes = nnodes;
//        xsum_data[tid].NT = NT;
//        xsum_data[tid].tid = tid;
//        pthread_create(&threads[tid], NULL, xsum_f, &xsum_data[tid]);
//    }
//    double xsum = 0.0;
//    for (uint tid = 0; tid < NT; ++tid) {
//        pthread_join(threads[tid], NULL);
//        xsum += xsum_helper[tid];
//    }
//    free(xsum_data);
//    free(xsum_helper);
//
//    //normalise
//    struct xnormalise_data_t* xnormalise_data = (struct xnormalise_data_t*) calloc(NT, sizeof(struct xnormalise_data_t));
//    for (uint tid = 0; tid < NT; ++tid) {
//        xnormalise_data[tid].xsum = xsum;
//        xnormalise_data[tid].outvec = outvec;
//        xnormalise_data[tid].nnodes = nnodes;
//        xnormalise_data[tid].NT = NT;
//        xnormalise_data[tid].tid = tid;
//        pthread_create(&threads[tid], NULL, xnormalise_f, &xnormalise_data[tid]);
//    }
//    for (uint tid = 0; tid < NT; ++tid) {
//        pthread_join(threads[tid], NULL);
//    }
//    free(xnormalise_data);

    //ps
    const size_t ps_len = 1+reps[0]->maxLevel;
    size_t* pss_orig = (size_t*)calloc(NT*ps_len, sizeof(size_t)); //0-init
    size_t* pss = (size_t*)calloc(NT*ps_len, sizeof(size_t)); //0-init
    if (pss_orig == NULL || pss == NULL) {
        printf("Memory allocation error.\n");
        exit(2);
    }
    struct compute_ps_data_t* compute_ps_data = (struct compute_ps_data_t*) calloc(NT, sizeof(struct compute_ps_data_t));
    for (uint tid = 0; tid < NT; ++tid) {
        compute_ps_data[tid].ps = pss_orig+(tid*ps_len);
        compute_ps_data[tid].rep = reps[tid];
        pthread_create(&threads[tid], NULL, compute_ps_f, &compute_ps_data[tid]);
    }
    for (uint tid = 0; tid < NT; ++tid) {
        pthread_join(threads[tid], NULL);
    }
    free(compute_ps_data);

    //matdim
    size_t matdim=1;
    for(size_t l=0; l<reps[0]->maxLevel; ++l) {
        matdim *= K;
    }

    //outdegree
    u_int32_t *outdeg = (u_int32_t *) calloc(reps[0]->numberOfNodes, sizeof(u_int32_t));
    {
        FILE *file;
        file = fopen(argv[2], "rb");
        if (file == NULL) {
            perror("Error while opening outdeg file");
            exit(1);
        }
        // Calcolo della dimensione del file
        fseek(file, 0, SEEK_END);
        size_t fileSize = ftell(file);
        fseek(file, 0, SEEK_SET);
        // Calcolo del numero di double nel file
        size_t filelen = fileSize / sizeof(u_int32_t);
        assert(filelen == nnodes);
        // Lettura dei dati dal file
        if (fread(outdeg, sizeof(u_int32_t), filelen, file) != filelen) {
            perror("Error while reading outdeg file");
            exit(1);
        };
        // Liberazione della memoria e chiusura del file
        fclose(file);
    }

//    uint *outdeg_helper = (u_int32_t *) calloc(NT*reps[0]->numberOfNodes, sizeof(uint));
//    struct compute_outdeg_data_t* compute_outdeg_data = (struct compute_outdeg_data_t*) calloc(NT, sizeof(struct compute_outdeg_data_t));
//    for (uint tid = 0; tid < NT; ++tid) {
//        compute_outdeg_data[tid].rep = reps[tid];
//        compute_outdeg_data[tid].outdeg = outdeg_helper + (tid*nnodes);
//        compute_outdeg_data[tid].ps = pss + (tid*ps_len);
//        compute_outdeg_data[tid].ps_orig = pss_orig + (tid*ps_len);
//        compute_outdeg_data[tid].dim = matdim;
//        pthread_create(&threads[tid], NULL, compute_outdeg_f, &compute_outdeg_data[tid]);
//    }
//    for (uint tid = 0; tid < NT; ++tid) {
//        pthread_join(threads[tid], NULL);
//        for(uint c=0; c<reps[0]->numberOfNodes; ++c) {
//            outdeg[c] += outdeg_helper[tid*nnodes+c];
//        }
//    }
//    free(outdeg_helper);
//    free(compute_outdeg_data);

//    debug_uint_vec(outdeg, rep->numberOfNodes, "OUTDEG");
//    debug_bitmat(rep);
//    debug_size_t_vec(ps_orig, ps_len, "PS");

    //loop's start
    double *contrib_dn_helper = (double *) calloc(NT, sizeof(double));
    struct loop_data_t* loop_data = (struct loop_data_t*) calloc(NT, sizeof(struct loop_data_t));
    for (uint tid = 0; tid < NT; ++tid) {
        loop_data[tid].contrib_dn_helper = contrib_dn_helper;
        loop_data[tid].nnodes = nnodes;
        loop_data[tid].outdeg = outdeg;
        loop_data[tid].invec = invec;
        loop_data[tid].outvec = outvec;
        loop_data[tid].ps_len = ps_len;
        loop_data[tid].ps = pss + (tid*ps_len);
        loop_data[tid].ps_orig = pss_orig + (tid*ps_len);
        loop_data[tid].dim = matdim;
        loop_data[tid].rep = reps[tid];
        loop_data[tid].NT = NT;
        loop_data[tid].tid = tid;
    }
    for(size_t iter=0; iter < NITERS; ++iter) {
//        debug_uint_vec(outdeg, nnodes, "outdeg");
//        debug_double_vec(loop_data[0].outvec, nnodes, "OUTVEC prima");

        //contrib dn
        for (uint tid = 0; tid < NT; ++tid) {
            pthread_create(&threads[tid], NULL, contrib_dn_f, &loop_data[tid]);
        }
        for (uint tid = 0; tid < NT; ++tid) {
            pthread_join(threads[tid], NULL);
            if(tid) contrib_dn_helper[0] += contrib_dn_helper[tid]; //reduce
        }

//        debug_double_vec(loop_data[0].outvec, nnodes, "OUTVEC dopo");

        //right multiplication
        for (uint tid = 0; tid < NT; ++tid) {
            pthread_create(&threads[tid], NULL, right_multiply_f, &loop_data[tid]);
        }
        for (uint tid = 0; tid < NT; ++tid) {
            pthread_join(threads[tid], NULL);
        }

//        debug_double_vec(outvec, rep->numberOfNodes, "PART OUTVEC");
//        debug_double_vec(invec, rep->numberOfNodes, "PART INVEC");

        //finalisation
        for (uint tid = 0; tid < NT; ++tid) {
            pthread_create(&threads[tid], NULL, finalise_f, &loop_data[tid]);
        }
        for (uint tid = 0; tid < NT; ++tid) {
            pthread_join(threads[tid], NULL);
        }

    }//end loop
    free(contrib_dn_helper);
    free(threads);
    for(uint tid=0; tid<NT; ++tid) destroyRepresentation(reps[tid]);
    free(reps);

    debug_double_vec(loop_data[0].outvec, nnodes, "OUTVEC");

    // retrieve topk nodes
    double *_outvec = &(loop_data[0].outvec[0]);
    unsigned topk = TOPK;
    if (topk>nnodes) topk=nnodes;
    unsigned *top = (unsigned *) calloc(topk, sizeof(*top));
    unsigned *aux = (unsigned *) calloc(topk, sizeof(*top));
    if(top==NULL || aux==NULL){
        perror("Cannot allocate topk/aux array");
        exit(-1);
    }
    kLargest(_outvec,aux,nnodes,topk);
    // get sorted nodes in top
    for(long int i=topk-1;i>=0;i--) {
        top[i] = aux[0];
        aux[0] = aux[i];
        minHeapify(_outvec,aux,i,0);
    }
    // report topk nodes id's only on stdout
    fprintf(stdout,"Top:");
    for(int i=0;i<topk;i++) fprintf(stdout," %d",top[i]);
    fprintf(stdout,"\n");
    //deallocate
    free(top);
    free(aux);

    free(invec);
    free(outvec);
    free(loop_data);
    return 0;
}


