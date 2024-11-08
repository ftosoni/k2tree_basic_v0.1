#include "pagerank_utils.h"

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

void right_multiply(const MREP *rep, const double* invec, double* outvec, size_t *ps, const size_t dim) {
    right_multiply_helper(rep, invec, outvec,
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

static void usage_and_exit(char *name)
{
    fprintf(stderr,"Usage:\n\t  %s [options] matrix col_count_file\n",name);
    fprintf(stderr,"\t\t-v             verbose\n");
//    fprintf(stderr,"\t\t-b num         number of row blocks, def. 1\n");
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
    int maxiter=100,topk=3,nblocks=1;
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
//            case 'b':
//                nblocks=atoi(optarg); break;
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
    if(maxiter<1 || nblocks < 1 || topk<1) {
        fprintf(stderr,"Error! Options -b -m and -k must be at least one\n");
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

    if (argc != 2+1) {
        fprintf(stderr, "USAGE: %s <GRAPH> <count col file>\n", argv[0]);
        exit(-1);
    }

    //bind to core 0
    {
        const int ncores = sysconf(_SC_NPROCESSORS_ONLN); // Number of available cores; not really used
        pthread_t main_thread = pthread_self(); // Get the identifier of the calling threadft
        int tid = 0;
        set_core(&main_thread,tid,ncores);
    }

    // ----------- business logic from here
    MREP *rep = loadRepresentation(argv[1]);

    if(verbose>0) {
        fprintf(stderr,"Number of nodes: %d\n", rep->numberOfNodes);
        fprintf(stderr,"Number of arcs: %ld\n", rep->numberOfEdges);
    }
//    debug_mrep(rep);
//    debug_bitmat(rep);

    double* invec = (double*)calloc(rep->numberOfNodes, sizeof(double));
    double* outvec = (double*)calloc(rep->numberOfNodes, sizeof(double));

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
        const size_t e = fread(outdeg,sizeof(u_int32_t),rep->numberOfNodes,ccol_file);
        assert(e == rep->numberOfNodes);
    }

//    debug_uint_vec(outdeg, rep->numberOfNodes, "OUTDEG");
//    debug_bitmat(rep);
//    debug_size_t_vec(ps_orig, ps_len, "PS");

    for(size_t iter=0; iter < maxiter; ++iter) {
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
            outvec[r] = dampf * outvec[r] + (1-dampf) / rep->numberOfNodes; //teleporting
        }

    }

//    debug_double_vec(outvec, rep->numberOfNodes, "OUTVEC");
    free(invec);

    if(verbose>0) {
        double sum = 0.0;
        for(int i=0;i<rep->numberOfNodes;i++) sum += outvec[i];
        fprintf(stderr,"Sum of ranks: %f (should be 1)\n",sum);
    }

    // retrieve topk nodes
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
    // report topk nodes sorted by decreasing rank
    if (verbose>0) {
        fprintf(stderr, "Top %d ranks:\n",topk);
        for(int i=0;i<topk;i++) fprintf(stderr,"  %d %lf\n",top[i],outvec[top[i]]);
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


