#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include <igraph.h>
#include <time.h>
#include "random.h"
#include "cluster.h"
#include "extras.h"

#define kmin 1		/* minimum cutoff */
#define kmax 100	/* maximum cutoff */
#define PTS  12		/* number of kappa's */

double kappa;		/* exponential cutoff */
double tau;			/* power law decay */
igraph_integer_t nodes;			/* number of nodes */


/*
 * Routine generating the degree sequence with distribution
 * 		(k^tau)*exp(-k/kappa),
 * with 'tau' and 'kappa' globally defined.
 * Routine checks that the sum of the degrees is even; if not, routine
 * is repeated
 */
void gendegs (igraph_vector_t *degs){
	
	int i;		/* cycle variables */
	int check;		/* 1 when random number accepted, 0 when not */
	double u[2];	/* auxiliary variables: 2 random reals in [0,1) */
	int sum=0;		/* sum of degrees; check it to be even at the end */
	
	for(i=0; i < igraph_vector_size(degs); i++){
		check = 0;
		while(check != 1){
			ranlxd(u,2);
			
			/* generate random integer with exponential decay */
			VECTOR(*degs)[i] = (int)(-kappa*log(1-u[0]));
			
			/* accept number with power law probability */
			if (u[1] < 1/pow(VECTOR(*degs)[i],tau))
				check = 1;
			
			/* printf("%d\n", check);*/
		};
		sum += VECTOR(*degs)[i];
	}
	
	if(sum%2 == 1)
		gendegs(degs);
	else return;
}



int main(int argc, char **argv){

/************* CHECK NUMBER OF ARGUMENTS ******************************/

	if(argc<4){
		printf("Set following parameters:\n \
	1) number of vertices;\n \
	2) power law decay;\n \
	3) number of samples.\n");
		exit(EXIT_FAILURE);
	}
	
	
/************* DECLARATION OF VARIABLES *******************************/
	
	FILE *out;			/* output file pointer */
	char *outstring;	/* output file name */
	
	int i, j, rep;		/* counters in cycles */
	int REPS;			/* number of iteration (sample size) */
	
	igraph_t graph;	/* pointer to the main graph object */
	
	igraph_vector_t degrees;	/* degree sequence of 'graph' */
	igraph_vector_t sizes;		/* vector of component sizes */
	
	cluster *gsize; 	/* jackknife clusters for largest comps sizes */
	double logstep;		/* kappa *= logstep */

	
/************* INITIALIZATION/ASSIGNMENT OF MAIN VARIABLES ************/
	
	/* main arguments */
	nodes	= atoi(argv[1]);
	tau		= atof(argv[2]);
	REPS	= atoi(argv[3]);
	
	/* output file name and pointer */
	outstring = malloc(100*sizeof(char));
	sprintf(outstring, "output/giant_nodes=%d_tau=%1.2lf", nodes, tau);
	out = fopen(outstring, "w");
	
	logstep = pow(10, (log10(kmax)-log10(kmin))/(double)(PTS-1));
	
	/* Jackknife clusters containing sample, average and std. dev */
	gsize = malloc(PTS*sizeof(cluster));
	for(i=0; i<PTS; i++)
		cluster_init(gsize+i,REPS);

	
/************* CALCULATION OF BIGGEST COMPONENT SIZE ******************/

	/* initialization of the RANLUX random number generator */
	srand(time(NULL));
	rlxd_init(1,rand());
	
	printf("tau\n%2.1lf\n", tau);		
	printf("\tkappa\n");
	kappa = kmin;
	for(i=0; i<PTS; i++){
		printf("\t%3.3lf\n", kappa);
		for(rep=0; rep<REPS; rep++){
			
			/* initialization igraph vectors */
			igraph_vector_init(&degrees, nodes);
			igraph_vector_init(&sizes, nodes);
			
			/* generation of degrees */
			gendegs(&degrees);
			
			/* creation of random graph with fixed degree sequence */
			igraph_degree_sequence_game(&graph, &degrees, 0,\
				IGRAPH_DEGSEQ_SIMPLE);
		
			/* size of largest component */
			igraph_clusters(&graph, NULL, &sizes, NULL, IGRAPH_WEAK);
			
			/* assign (relative) of largest connected component to a
			 * jackknife vector */
			gsize[i].Vec[rep] = igraph_vector_max(&sizes)/(double)nodes;
			
			/* destroy graph and free space */
			igraph_destroy(&graph);
			igraph_vector_destroy(&degrees);
			igraph_vector_destroy(&sizes);
			
		}
		/* print average maximum component size and std. dev */
		clusterJK(&gsize[i]);
		fprintf(out, "%lf\t%lf\t%lf\n"\
			, kappa, gsize[i].Mean, gsize[i].Sigma);
		
		kappa *= logstep;
	}
	
	fclose(out);
	exit(EXIT_SUCCESS);
}

