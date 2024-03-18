#include <stdio.h>
#include <math.h>
#include <string.h>
#include "kTree.h"



int main(int argc, char* argv[]){

	if(argc<2){
		fprintf(stderr,"USAGE: %s <GRAPH>\n",argv[0]);
		return(-1);
	}

	char *filename = (char *)malloc(sizeof(char)*(strlen(argv[1])+5));
	MREP * rep = loadRepresentation(argv[1]);
	
  strcpy(filename,argv[1]);
  strcat(filename,".crb");
	FILE *fr = fopen(filename,"w");
	fwrite(&(rep->numberOfNodes),sizeof(uint),1,fr);
  fwrite(&(rep->numberOfEdges),sizeof(ulong),1,fr);

	uint node, i, j, edge;
 
	for(i=0;i<rep->numberOfNodes;i++){
  	node = -(i+1);
  	fwrite(&node,sizeof(uint),1,fr);
		for(j=0;j<rep->numberOfNodes;j++){
			if(compact2CheckLinkQuery(rep, i, j))	{
		  	edge = j + 1;
		  	fwrite(&edge,sizeof(uint),1,fr);
		  }
  	}
  }
  
  fclose(fr);
 
  destroyRepresentation(rep);
  free(filename);
  
  return 0;
}


