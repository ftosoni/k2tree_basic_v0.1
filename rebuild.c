#include <stdio.h>
#include <math.h>
#include <string.h>
#include "kTree.h"



int main(int argc, char* argv[]){

	if(argc<2){
		fprintf(stderr,"USAGE: %s <GRAPH>\n",argv[0]);
		return(-1);
	}

	char *filename = (char *)malloc(sizeof(char)*(strlen(argv[1])+4));
	MREP * rep = loadRepresentation(argv[1]);
	
  strcpy(filename,argv[1]);
  strcat(filename,".rb");
	FILE *fr = fopen(filename,"w");
	 fwrite(&(rep->numberOfNodes),sizeof(uint),1,fr);
  fwrite(&(rep->numberOfEdges),sizeof(ulong),1,fr);

	uint node, i, j, edge;
	
	
	//Rebuilding using successors list
	
	uint * listady;
	for(i=0;i<rep->numberOfNodes;i++){
  	listady = compactAdjacencyList(rep, i);
  	node = -(i+1);
  	fwrite(&node,sizeof(uint),1,fr);
  	for(j=0;j<listady[0];j++){
  		edge = listady[j+1] + 1;
  		fwrite(&edge,sizeof(uint),1,fr);
  	}
  }
  
  
  //Rebuilding using range function
  /*
  uint ** listady2;
	for(i=0;i<rep->numberOfNodes;i++){
  	listady2 = compactRangeQuery(rep, i, i, 0,rep->numberOfNodes-1 );
  	node = -(i+1);
  	fwrite(&node,sizeof(uint),1,fr);
  	for(j=0;j<listady2[0][0];j++){
  		edge = listady2[1][j+1] + 1;
  		fwrite(&edge,sizeof(uint),1,fr);
  	}
  }
  */
  
  fclose(fr);
  
  destroyRepresentation(rep);
  free(filename);
  
  return 0;
}


