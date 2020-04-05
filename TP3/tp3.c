#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>//to estimate the runing time

#define EXIT_FAILURE 1


#define EXIT_FAILURE 1

#define NLINKS 100000000 //maximum number of edges for memory allocation, will increase if needed


typedef struct {
	unsigned long s;
	unsigned long t;
} edge;

//edge list structure:
typedef struct {
	unsigned long n;//number of nodes
	unsigned long e;//number of edges
	edge *edges;//list of edges
	unsigned long *cd;//cumulative degree cd[0]=0 length=n+1
	unsigned long *adj;//concatenated lists of neighbors of all nodes
} adjlist;

//compute the maximum of three unsigned long
inline unsigned long max3(unsigned long a,unsigned long b,unsigned long c){
	a=(a>b) ? a : b;
	return (a>c) ? a : c;
}

//reading the edgelist from file
adjlist* readedgelist(char* input){
	unsigned long e1=NLINKS;
	FILE *file=fopen(input,"r");

	adjlist *g=malloc(sizeof(adjlist));
	g->n=0;
	g->e=0;
	g->edges=malloc(e1*sizeof(edge));//allocate some RAM to store edges

	while (fscanf(file,"%lu %lu", &(g->edges[g->e].s), &(g->edges[g->e].t))==2) { // g.edges[e] = [s, t]
		g->n=max3(g->n,g->edges[g->e].s,g->edges[g->e].t);
		if (++(g->e)==e1) {//increase allocated RAM if needed // increment g->e by 1
			e1+=NLINKS;
			g->edges=realloc(g->edges,e1*sizeof(edge));
		}
	}
	fclose(file);

	g->n++;

	g->edges=realloc(g->edges,g->e*sizeof(edge));

	return g;
}

//building the adjacency matrix
void mkadjlist(adjlist* g){
	unsigned long i,u,v;
	unsigned long *d=calloc(g->n,sizeof(unsigned long)); // Degree array

	for (i=0;i<g->e;i++) { // build the degree array
		d[g->edges[i].s]++;
		d[g->edges[i].t]++;
	}

	g->cd=malloc((g->n+1)*sizeof(unsigned long)); // Cumulative degree array
	g->cd[0]=0;
	for (i=1;i<g->n+1;i++) { // build the cumulative degree array
		g->cd[i]=g->cd[i-1]+d[i-1];
		d[i-1]=0;
	}

	g->adj=malloc(2*g->e*sizeof(unsigned long)); // adjacency array (of length 2*m)

	for (i=0;i<g->e;i++) {
		u=g->edges[i].s;
		v=g->edges[i].t;
		g->adj[ g->cd[u] + d[u]++ ]=v; // g.adj[g.cd[u] + 1] = v then g.adj[g.cd[u] + 2] = next_v ... until the last neighbor of u
		g->adj[ g->cd[v] + d[v]++ ]=u;
	}

	free(d);
	//free(g->edges);
}


//freeing memory
void free_adjlist(adjlist *g){
	free(g->edges);
	free(g->cd);
	free(g->adj);
	free(g);
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////


// Exercise 1 — Simple bechmark


double random_uniform() { return ((double)rand() / (double)RAND_MAX); } // real number in [0,1]


void create_random_graph(char* input, double p, double q){ //
	FILE *file=fopen(input,"w");

	int cluster_1[100], cluster_2[100], cluster_3[100], cluster_4[100]; // The graph has 400 nodes partition into 4 clusters of 100 nodes.

	for (int i=0;i<100;i++) {
		cluster_1[i] = i;
		cluster_2[i] = 100 + i;
		cluster_3[i] = 200 + i;
		cluster_4[i] = 300 + i;
	}

	for (int i=0;i<99;i++) { // Each pair of nodes (i,j) in the same cluster is connected with a probability p
		for (int j=i+1;j<100;j++) {
			double u1, u2, u3, u4;
			u1 = random_uniform();
			u2 = random_uniform();
			u3 = random_uniform();
			u4 = random_uniform();
			if (u1 < p) {
				fprintf(file, "%d %d\n", cluster_1[i], cluster_1[j]);
			}
			if (u2 < p) {
				fprintf(file, "%d %d\n", cluster_2[i], cluster_2[j]);
			}
			if (u3 < p) {
				fprintf(file, "%d %d\n", cluster_3[i], cluster_3[j]);
			}
			if (u4 < p) {
				fprintf(file, "%d %d\n", cluster_4[i], cluster_4[j]);
			}
		}
	}

	for (int i=0;i<100;i++) { // Each pair of nodes in different clusters is connected with a probability q <= p
		for (int j=0;j<100;j++) {
			// cluster 1 with clusters 2, 3 and 4
			if (random_uniform() < q){fprintf(file, "%d %d\n", cluster_1[i], cluster_2[j]);}
			if (random_uniform() < q){fprintf(file, "%d %d\n", cluster_1[i], cluster_3[j]);}
			if (random_uniform() < q){fprintf(file, "%d %d\n", cluster_1[i], cluster_4[j]);}
			// cluster 2 with clusters 3 and 4
			if (random_uniform() < q){fprintf(file, "%d %d\n", cluster_2[i], cluster_3[j]);}
			if (random_uniform() < q){fprintf(file, "%d %d\n", cluster_2[i], cluster_4[j]);}
			// cluster 3 with cluster 4
			if (random_uniform() < q){fprintf(file, "%d %d\n", cluster_3[i], cluster_4[j]);}
		}
	}

	fclose(file);

}




// Exercise 2 — Label propagation

void print_array(unsigned long *arr, unsigned long n){ // print an array of length n
	for (unsigned long i=0;i<n;i++) {
		printf("%lu\n", arr[i]);
	}
	printf("\n");
}

void *fisher_yates_shuffle(unsigned long *t, unsigned long n) { // Fisher–Yates shuffle

	unsigned long j, temp;
	for (unsigned long i=n-1;i>0;i--){
		j = rand() % (i+1); // random integer such that 0 ≤ j ≤ i
		temp = t[i]; // exchange a[j] and a[i]
		t[i] = t[j];
		t[j] = temp;
	}
}


bool check_label_frequency(adjlist *g, unsigned long *part) { // Check if there exists a node with a label that does not have the highest frequency among its neighbours.

	for (unsigned long u=0;u<g->n;u++) {

		unsigned long degree = g->cd[u+1] - g->cd[u]; // degree of node u
		unsigned long labels[degree]; // unique labels of the neighbors
		unsigned long counts[degree]; // count of each labels
		unsigned long nb_of_labels = 1;
		unsigned long most_freq_label = part[u];
		unsigned long count_most_freq_label = 1;
		labels[0] = part[u]; // the first label is the one of the node u
		counts[0] = 1;

		for (int i=0;i<degree;i++) { // for each neighbor v of u in the graph g
			unsigned long v = g->adj[ g->cd[u] + i ];

			bool label_exists = 0;
			for (int j=0;j<nb_of_labels;j++) {
				if (part[v] == labels[j]) { // more than one node with this label
					counts[j] += 1;
					label_exists = 1;
					break;
				}
			}
			if (label_exists == 0) { // first time that this label appears among the neighbors
				nb_of_labels += 1;
				labels[nb_of_labels-1] = part[v];
				counts[nb_of_labels-1] = 1;
			}
		}

		for (int j=0;j<nb_of_labels;j++) { // determine the most frequent label
			if (count_most_freq_label <= counts[j]) {
				most_freq_label = labels[j];
				count_most_freq_label = counts[j];
			}
		}

		if (part[u] != most_freq_label) { // there exists a node with a label that does not have the highest frequency among its neighbours.
			return 0;
		}

	}

	return 1; // There doesn't exist a node with a label that does not have the highest frequency among its neighbours.

}


void label_propagation(adjlist* g, unsigned long *part) { // Label propagation algorithm

	unsigned long *order = malloc(g->n * sizeof(unsigned long));

	for (unsigned long i=0;i<g->n;i++) { // Step 1 : give a unique label to each node in the network
		part[i] = i;
		order[i] = i;
	}

	while(check_label_frequency(g, part)==0) {
		fisher_yates_shuffle(order, g->n); // Step 2: Arrange the nodes in the network in a random order
		for (unsigned long u=0;u<g->n;u++) { 
// Step 3: for each node in the network (in this random order) set its label to a label occurring with the highest frequency among its neighbours
			unsigned long new_u = order[u];

			unsigned long degree = g->cd[new_u+1] - g->cd[new_u]; // degree of node new_u
			unsigned long labels[degree]; // unique labels of the neighbors
			unsigned long counts[degree]; // count of each labels
			unsigned long nb_of_labels = 1;
			unsigned long most_freq_label = part[new_u];
			unsigned long count_most_freq_label = 1;
			labels[0] = part[new_u]; // the first label is the one of the node new_u
			counts[0] = 1;

			for (int i=0;i<degree;i++) { // for each neighbor v of new_u in the graph g
				unsigned long v = g->adj[ g->cd[new_u] + i ];

				bool label_exists = 0;
				for (int j=0;j<nb_of_labels;j++) {
					if (part[v] == labels[j]) { // more than one node with this label
						counts[j] += 1;
						label_exists = 1;
						break;
					}
				}
				if (label_exists == 0) { // first time that this label appears in the neighbourhood
					nb_of_labels += 1;
					labels[nb_of_labels-1] = part[v];
					counts[nb_of_labels-1] = 1;
				}
			}

			for (int j=0;j<nb_of_labels;j++) { // determine the most frequent label
				if (count_most_freq_label <= counts[j]) {
					most_freq_label = labels[j];
					count_most_freq_label = counts[j];
				}
			}

			part[new_u] = most_freq_label; // set its label to a label occurring with the highest frequency among its neighbours

		}
	}

	free(order);

}





int main(int *argc, char *argv[]){
	adjlist* g;
	time_t t1,t2,t3,t4,t5;

	t1=time(NULL);

	// Exercise 1 - Simple benchmark (run ./tp3 create_graphs)

	//////////////////////
	if (strcmp(argv[1],"create_graphs")==0) {
		printf("Creation of the 6 random graphs\n");
		create_random_graph("graph_0.8_0.4.txt", 0.8, 0.4); // ratio = 2
		create_random_graph("graph_0.8_0.2.txt", 0.8, 0.2); // ratio = 4
		create_random_graph("graph_0.8_0.1.txt", 0.8, 0.1); // ratio = 8
		create_random_graph("graph_0.8_0.08.txt", 0.8, 0.08); // ratio = 10
		create_random_graph("graph_0.8_0.05.txt", 0.8, 0.05); // ratio = 16
		create_random_graph("graph_0.8_0.025.txt", 0.8, 0.025); // ratio = 32
	}
	//////////////////////



	// Exercise 2 - Label propagation (run ./tp3 graph.txt comm_graph.txt e.g. with graph = graph_0.8_0.025)

	//////////////////////
	else {
		printf("Reading edgelist from file %s\n",argv[1]);
		g=readedgelist(argv[1]);

		printf("Number of nodes: %lu\n",g->n);
		printf("Number of edges: %lu\n",g->e);

		printf("Building the adjacency list\n");
		mkadjlist(g);
		printf("Running the label propagation algorithm\n");
		unsigned long *part = malloc(g->n * sizeof(unsigned long));
		t2 = time(NULL);
		label_propagation(g ,part);
		t3 = time(NULL);
	  	printf("- Time to compute communities = %ldh%ldm%lds\n", (t3-t2)/3600, ((t3-t2)%3600)/60, ((t3-t2)%60));


		printf("Prints result in file %s\n", argv[2]);
	  	FILE* out = fopen(argv[2], "w");
	  	for(unsigned long i = 0; i < g->n; i++){
	    	fprintf(out, "%lu %lu\n", i, part[i]);
	  	}
	  	fclose(out);
	  	t4 = time(NULL);
	 	printf("- Time to export communities = %ldh%ldm%lds\n", (t4-t3)/3600, ((t4-t3)%3600)/60, ((t4-t3)%60));
		free_adjlist(g);
	}
	//////////////////////


	t5=time(NULL);

	printf("- Overall time = %ldh%ldm%lds\n",(t5-t1)/3600,((t5-t1)%3600)/60,((t5-t1)%60));

	return 0;
}