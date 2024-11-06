#define _GNU_SOURCE

#include "pagerank_utils.h"
#include <sched.h>
#include <pthread.h>

inline void set_core(pthread_t tid, int core) {
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(core, &cpuset);
    pthread_setaffinity_np(tid, sizeof(cpu_set_t), &cpuset);
}

inline int my_pthread_create(pthread_t *thread, const pthread_attr_t *attr, void *(*start_routine)(void *), void *arg, int core) {
    int result = pthread_create(thread, attr, start_routine, arg);
    if (result == 0) {
        set_core(*thread, core);
    }
    return result;
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

struct fetch_data_t {
    char *argv1;
    int NT;
    MREP **reps;
    uint tid;
};

void* fetch_f(void *arg){
    struct fetch_data_t* data = (struct fetch_data_t*)arg;

    //construct the infilepath

    const size_t len = snprintf(NULL, 0, "%s.%d.%d", data->argv1, data->NT, data->tid);
    char * fpath = (char *) malloc (len + 1);

    const int res = snprintf(fpath, len+1, "%s.%d.%d", data->argv1, data->NT, data->tid);
    assert(res == len);

//    printf("%d: %s\n", data->tid, fpath);
    data->reps[data->tid] = loadRepresentation(fpath);
    free(fpath);
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

//loop from here

struct loop_data_t {
    double *contrib_dn_helper;
    size_t nnodes;
    uint *outdeg;
    double *invec;
    double *outvec;
    double dampf;
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
        data->outvec[r] = data->dampf * data->outvec[r] + (1-data->dampf) / data->rep->numberOfNodes; //teleporting
    }
    return NULL;
}

static void usage_and_exit(char *name)
{
    fprintf(stderr,"Usage:\n\t  %s [options] matrix col_count_file\n",name);
    fprintf(stderr,"\t\t-v             verbose\n");
    fprintf(stderr,"\t\t-b num         parallelism degree, def. 2\n");
    fprintf(stderr,"\t\t-m maxiter     maximum number of iteration, def. 100\n");
//    fprintf(stderr,"\t\t-e eps         stop if error<eps (default ignore error)\n");
    fprintf(stderr,"\t\t-d df          damping factor (default 0.9)\n");
    fprintf(stderr,"\t\t-k K           show top K nodes (default 3)n\n");
    exit(1);
}

int main(int argc, char* argv[]) {
    extern char *optarg;
    extern int optind, opterr, optopt;
    int verbose=0;
    int c;
//    time_t start_wc = time(NULL);

    // default values for command line parameters
    int maxiter=100,topk=3,NT=2;
    double dampf = 0.9;//, eps = -1;

    /* ------------- read options from command line ----------- */
    opterr = 0;
    while ((c=getopt(argc, argv, "b:m:d:k:v")) != -1) {
        switch (c)
        {
            case 'v':
                verbose++; break;
            case 'm':
                maxiter=atoi(optarg); break;
            case 'b':
                NT=atoi(optarg); break;
//            case 'e':
//                eps=atof(optarg); break;
            case 'd':
                dampf=atof(optarg); break;
            case 'k':
                topk=atoi(optarg); break;
            case '?':
                fprintf(stderr,"Unknown option: %c\n", optopt);
                exit(1);
        }
    }
    if(verbose>0) {
        fputs("==== Command line:\n",stderr);
        for(int i=0;i<argc;i++)
            fprintf(stderr," %s",argv[i]);
        fputs("\n",stderr);
    }
    // check command line
    if(maxiter<1 || topk<1) {
        fprintf(stderr,"Error! Options -m and -k must be at least one\n");
        usage_and_exit(argv[0]);
    }
    if(NT<2) {
        fprintf(stderr,"Error! Options -b must be at least two\n");
        usage_and_exit(argv[0]);
    }
    if(dampf<0 || dampf>1) {
        fprintf(stderr,"Error! Options -d must be in the range [0,1]\n");
        usage_and_exit(argv[0]);
    }
    // virtually get rid of options from the command line
    optind -=1;
    if (argc-optind != 3) usage_and_exit(argv[0]);
    argv += optind; argc -= optind;

    // ----------- business logic from  here
    //data
    pthread_t *threads = (pthread_t *) calloc(NT, sizeof(pthread_t));
    MREP **reps = (MREP **) calloc(NT, sizeof(MREP *));
    assert(reps != NULL);

    //read matrices
    struct fetch_data_t *fetch_data = (struct fetch_data_t *) calloc(NT, sizeof(struct fetch_data_t));
    {
        for (uint tid = 0; tid < NT; ++tid) {
            fetch_data[tid].argv1 = argv[1];
            fetch_data[tid].NT = NT;
            fetch_data[tid].reps = reps;
            fetch_data[tid].tid = tid;
            my_pthread_create(&threads[tid], NULL, fetch_f, &fetch_data[tid], tid);
        }
        for (uint tid = 0; tid < NT; ++tid) {
            pthread_join(threads[tid], NULL);
            assert(reps[0]->numberOfNodes == reps[tid]->numberOfNodes);
        }
    }
    free(fetch_data);
    if(verbose>0) {
        fprintf(stderr,"Number of nodes: %d\n", reps[0]->numberOfNodes);
        size_t nedges = 0;
        for (uint tid = 0; tid < NT; ++tid) nedges += reps[tid]->numberOfEdges;
        fprintf(stderr,"Number of arcs: %ld\n", nedges);
    }

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
        my_pthread_create(&threads[tid], NULL, compute_ps_f, &compute_ps_data[tid], tid);
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
        loop_data[tid].dampf = dampf;
        loop_data[tid].ps_len = ps_len;
        loop_data[tid].ps = pss + (tid*ps_len);
        loop_data[tid].ps_orig = pss_orig + (tid*ps_len);
        loop_data[tid].dim = matdim;
        loop_data[tid].rep = reps[tid];
        loop_data[tid].NT = NT;
        loop_data[tid].tid = tid;
    }
    for(size_t iter=0; iter < maxiter; ++iter) {
//        debug_uint_vec(outdeg, nnodes, "outdeg");
//        debug_double_vec(loop_data[0].outvec, nnodes, "OUTVEC prima");

        //contrib dn
        for (uint tid = 0; tid < NT; ++tid) {
            my_pthread_create(&threads[tid], NULL, contrib_dn_f, &loop_data[tid], tid);
        }
        for (uint tid = 0; tid < NT; ++tid) {
            pthread_join(threads[tid], NULL);
            if(tid) contrib_dn_helper[0] += contrib_dn_helper[tid]; //reduce
        }

//        debug_double_vec(loop_data[0].outvec, nnodes, "OUTVEC dopo");

        //right multiplication
        for (uint tid = 0; tid < NT; ++tid) {
            my_pthread_create(&threads[tid], NULL, right_multiply_f, &loop_data[tid], tid);
        }
        for (uint tid = 0; tid < NT; ++tid) {
            pthread_join(threads[tid], NULL);
        }

//        debug_double_vec(outvec, rep->numberOfNodes, "PART OUTVEC");
//        debug_double_vec(invec, rep->numberOfNodes, "PART INVEC");

        //finalisation
        for (uint tid = 0; tid < NT; ++tid) {
            my_pthread_create(&threads[tid], NULL, finalise_f, &loop_data[tid], tid);
        }
        for (uint tid = 0; tid < NT; ++tid) {
            pthread_join(threads[tid], NULL);
        }

    }//end loop
    free(contrib_dn_helper);
    free(threads);
    for(uint tid=0; tid<NT; ++tid) destroyRepresentation(reps[tid]);
    free(reps);

//    debug_double_vec(loop_data[0].outvec, nnodes, "OUTVEC");

    if(verbose>0) {
        double sum = 0.0;
        for(int i=0;i<nnodes;i++) sum += loop_data[0].outvec[i];
        fprintf(stderr,"Sum of ranks: %f (should be 1)\n",sum);
    }

    // retrieve topk nodes
    double *_outvec = &(loop_data[0].outvec[0]);
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
    // report topk nodes sorted by decreasing rank
    if (verbose>0) {
        fprintf(stderr, "Top %d ranks:\n",topk);
        for(int i=0;i<topk;i++) fprintf(stderr,"  %d %lf\n",top[i],loop_data[0].outvec[top[i]]);
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
