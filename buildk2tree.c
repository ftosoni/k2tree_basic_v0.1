#include <stdio.h>
#include <math.h>
#include <string.h>
#include "kTree.h"

#include <assert.h>

int main(int argc, char* argv[]) {
    FILE *f;
    uint nodes;
    ulong edges;
    register ulong i;

    if (argc != 1+1) {
        fprintf(stderr, "USAGE: %s <GRAPH no ext>\n", argv[0]);
        return (-1);
    }

    //open file
    {

        const size_t len = snprintf(NULL, 0, "%s.kt-plain", argv[1]);
        char * kt_plain_fpath = (char *) malloc (len + 1);

        int res = snprintf(kt_plain_fpath, len+1, "%s.kt-plain", argv[1]);
        assert(res == len);

        f = fopen(kt_plain_fpath, "r");
        if (f == NULL) {
            printf("Failed to open file %s.\n", kt_plain_fpath);
            return 1;
        }
        free(kt_plain_fpath);
    }
    assert(f is not NULL);

    fread(&nodes, sizeof(uint), 1, f);

    uint max_level = floor(log(nodes) / log(K));
    if (floor(log(nodes) / log(K)) == (log(nodes) / log(K))) {
        max_level = max_level - 1;
    }
    fread(&edges, sizeof(ulong), 1, f);
    uint nodes_read = 0;

    uint *xedges = (uint *) malloc(sizeof(uint) * edges);
    uint *yedges = (uint *) malloc(sizeof(uint) * edges);
    uint cedg = 0;
    for (i = 0; i < nodes + edges; i++) {
        int k;
        fread(&k, sizeof(uint), 1, f);
        if (k < 0) {
            nodes_read++;
        } else {
            k--;
            assert(k < nodes);
            xedges[cedg] = nodes_read - 1;
            yedges[cedg] = k;
            cedg++;
        }
    }

    MREP *rep;

    rep = compactCreateKTree(xedges, yedges, nodes, edges, max_level);

/*
 	rep->info = (uint *)malloc(sizeof(uint)*MAX_INFO);
	rep->element = (uint *)malloc(sizeof(uint)*MAX_INFO);	
	rep->basex = (uint *)malloc(sizeof(uint)*MAX_INFO);
	rep->basey = (uint *)malloc(sizeof(uint)*MAX_INFO);
	rep->iniq = -1;
	rep->finq =-1;
		rep->div_level_table = (uint *)malloc(sizeof(uint)*rep->maxLevel);
	for(i=0;i<rep->maxLevel;i++)
		rep->div_level_table[i]=exp_pow(K,rep->maxLevel-i);
*/

    saveRepresentation(rep, argv[1]);
    destroyRepresentation(rep);

    free(xedges);
    free(yedges);
    fclose(f);
    return 0;
}


