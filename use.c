#include <stdio.h>
#include <math.h>
#include <string.h>
#include "kTree.h"



int main(int argc, char* argv[]){

	if(argc<7){
		fprintf(stderr,"USAGE: %s <GRAPH> <param> <p1> <p2> <p3> <p4>\n",argv[0]);
		return(-1);
	}

	//char *filename = (char *)malloc(sizeof(char)*20);
	
	MREP * rep = loadRepresentation(argv[1]);
	fprintf(stderr,"Adjacency list of node %d:\n",atoi(argv[2]));
	uint i;
  uint * listady = compact2AdjacencyList(rep, atoi(argv[2]));
  fprintf(stderr,"Number of neighbours: %d\n",listady[0]);
  for(i=0;i<listady[0];i++)
  	fprintf(stderr,"%d\t",listady[i+1]);
  fprintf(stderr,"\n");
 // numberNodesPerLevel(rep);

  fprintf(stderr,"--------------\n");
 	
 	listady = compactInverseList(rep, atoi(argv[2]));
  fprintf(stderr,"Number of inverse neighbours: %d\n",listady[0]);
  for(i=0;i<listady[0];i++)
  	fprintf(stderr,"%d\t",listady[i+1]);
  fprintf(stderr,"\n");  
  
  uint p1=atoi(argv[3]);
  uint p2=atoi(argv[4]);
  uint q1=atoi(argv[5]);
  uint q2=atoi(argv[6]);
  uint ** respuesta = compactRangeQuery(rep, p1, p2, q1, q2);
  
  fprintf(stderr,"--------------\n");
  fprintf(stderr,"Range: [%d,%d]-[%d,%d], total of links %d\n",p1,p2,q1,q2,respuesta[0][0]);
 for(i=0;i<respuesta[0][0];i++)
  	fprintf(stderr,"(%d,%d)\t",respuesta[0][i+1],respuesta[1][i+1]);  
	fprintf(stderr,"\n");
	
	
	fprintf(stderr,"Checking (%d,%d): %d\n",p1,q1,compact2CheckLinkQuery(rep,p1,q1));


	fprintf(stderr,"Checking [%d,%d]-[%d,%d]: %d\n",p1,p2,q1,q2,compactCheckRangeQuery(rep,p1,p2,q1,q2));
	
  
//  int j;
//  for(j=p1;j<=p2;j++){
//  listady = compact2AdjacencyList(rep, j);
//  fprintf(stderr,"Number of neighbours: %d\n",listady[0]);
//  for(i=0;i<listady[0];i++)
//  	fprintf(stderr,"%d\t",listady[i+1]);
//  fprintf(stderr,"\n");
//	}
  destroyRepresentation(rep);

  
  return 0;
}


