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

void right_multiply(const MREP *rep, const double* invec, double* outvec, size_t *ps, const size_t dim) {
    right_multiply_helper(rep, invec, outvec,
                          0, //current lvl
                          ps, //per-level offset
                          0, //row
                          0, //col,
                          dim //submatrix dimension at the current lvl
    );
}

void compute_outdeg(const MREP *rep, uint* outdeg, size_t *ps, const size_t dim) {
    compute_outdeg_helper(rep, outdeg,
                          0, //current lvl
                          ps, //per-level offset
                          0, //row
                          0, //col,
                          dim //submatrix dimension at the current lvl
    );
}

void compute_ps(const MREP *rep, size_t *ps){
    //assuming ps is 0-init
    assert(len <= rep->numberOfNodes());

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
}

int main(int argc, char* argv[]) {

    if (argc != 2+1) {
        fprintf(stderr, "USAGE: %s <GRAPH> <count col file>\n", argv[0]);
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

    //initialise input vector
    {
        for(size_t c=0; c<rep->numberOfNodes; ++c) {
            outvec[c] = 1.0/rep->numberOfNodes;
        }
    }

    //ps
    const size_t ps_len = 1+rep->maxLevel;
    size_t* ps_orig = (size_t*)calloc(ps_len, sizeof(size_t)); //0-init
    size_t* ps = (size_t*)calloc(ps_len, sizeof(size_t)); //0-init
    if (ps_orig == NULL) {
        printf("Memory allocation error.\n");
        exit(2);
    }
    compute_ps(rep, ps_orig);

    //matdim
    size_t matdim=1;
    for(size_t l=0; l<rep->maxLevel; ++l) {
        matdim *= K;
    }

    //outdegree
    u_int32_t *outdeg = (u_int32_t *) calloc(rep->numberOfNodes, sizeof(u_int32_t));
    FILE *ccol_file  = fopen(argv[2],"rb");
    {
        size_t e = fread(outdeg,sizeof(u_int32_t),rep->numberOfNodes,ccol_file);
        assert(e == rep->numberOfNodes);
    }

//    debug_uint_vec(outdeg, rep->numberOfNodes, "OUTDEG");
//    debug_bitmat(rep);
//    debug_size_t_vec(ps_orig, ps_len, "PS");

    for(size_t iter=0; iter < NITERS; ++iter) {
        //swap invec & outvec
        {
            double *tmp = outvec;
            outvec = invec;
            invec = tmp;
        }

        double contrib_dn = 0.0;
        for (uint x = 0; x < rep->numberOfNodes; ++x) {
            contrib_dn += (outdeg[x] == 0) * invec[x]; //contribution of dangling nodes
            if(outdeg[x]) invec[x] /= outdeg[x]; //normalisation
            outvec[x] = 0.0; //0-init
        }
        contrib_dn /= rep->numberOfNodes;

        //page rank
        memcpy(ps, ps_orig, ps_len * sizeof(size_t));
        right_multiply(rep, invec, outvec, ps, matdim);

//        debug_double_vec(outvec, rep->numberOfNodes, "PART OUTVEC");
//        debug_double_vec(invec, rep->numberOfNodes, "PART INVEC");

        for (uint r = 0; r < rep->numberOfNodes; ++r) {
            outvec[r] += contrib_dn; //dangling nodes
            outvec[r] = (1 - ALPHA) * outvec[r] + ALPHA / rep->numberOfNodes; //teleporting
        }

    }

    debug_double_vec(outvec, rep->numberOfNodes, "OUTVEC");
    free(invec);

    // retrieve topk nodes
    unsigned topk = TOPK;
    if (topk>rep->numberOfNodes) topk=rep->numberOfNodes;
    unsigned *top = (unsigned *) calloc(topk, sizeof(*top));
    unsigned *aux = (unsigned *) calloc(topk, sizeof(*top));
    if(top==NULL || aux==NULL){
        perror("Cannot allocate topk/aux array");
        exit(-1);
    }
    kLargest(outvec,aux,rep->numberOfNodes,topk);
    // get sorted nodes in top
    for(long int i=topk-1;i>=0;i--) {
        top[i] = aux[0];
        aux[0] = aux[i];
        minHeapify(outvec,aux,i,0);
    }
    // report topk nodes id's only on stdout
    fprintf(stdout,"Top:");
    for(int i=0;i<topk;i++) fprintf(stdout," %d",top[i]);
    fprintf(stdout,"\n");
    //deallocate
    free(top);
    free(aux);

    destroyRepresentation(rep);
    free(outvec);
    return 0;
}


