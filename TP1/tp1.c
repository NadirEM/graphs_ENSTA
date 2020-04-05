#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>//to estimate the runing time

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

typedef struct node node;
struct node
{
	unsigned long value;
	node *next;
};

typedef struct queue queue;
struct queue
{
	node *first;
};

inline unsigned long max(unsigned long a, unsigned long b) {
  return a > b ? a : b;
}

void print_array(unsigned long *arr, unsigned long n){ // print an array of length n
	for (unsigned long i=0;i<n;i++) {
		printf("%lu ", arr[i]);
	}
	printf("\n");
}

queue *init_queue(unsigned long s){ // initialisation of a queue with a node s
	queue *q = malloc(sizeof(*q));
	node *new = malloc(sizeof(*new));
	new->value = s;
	new->next = NULL;
	q->first = new;
	return q;
}




void add_queue(queue *q, unsigned long s){ // add node s to the queue q
	node *new = malloc(sizeof(*new));
	if(q == NULL || new == NULL)
	{
		printf("The queue is NULL\n");
		exit(EXIT_FAILURE);
	}

	new->value = s;
	new->next = NULL;

	if (q->first != NULL)  // More than one node in the queue
	{
		node *current_node = q->first;
		while (current_node->next != NULL)
		{
			current_node = current_node->next;
		}
		current_node->next = new;
	}
	else // No nodes in the queue, the node s will be the first node
	{
		q->first = new;
	}
}

unsigned long pop_queue(queue *q){ // delete the first node of the queue and return its value 
	if(q == NULL) {
		printf("The queue is NULL\n");
		exit(EXIT_FAILURE);
	}

	unsigned long v = -1; // return -1 if the queue is not NULL but the first node is NULL

	if(q->first != NULL){
		node *first_node = q->first;
		v = first_node->value;
		q->first = first_node->next;
		free(first_node);
	}

	return v;
}

void print_queue(queue *q){ // print all the queue values
	if(q == NULL) {
		printf("The queue is NULL\n");
		exit(EXIT_FAILURE);
	}

	node *current_node = q->first;

	while (current_node != NULL) {
		printf("value : %lu\n", current_node->value);
		current_node = current_node->next;
	}

	printf("\n");
}

void free_queue(queue *q){ // freeing the queue memory
	if(q == NULL) {
		printf("The queue is NULL\n");
		exit(0);
	}
	unsigned long temp = 0;
	while(q->first != NULL){
		temp = pop_queue(q);
	}
	free(q);
}


// Breadth First Search algorithm (slide 25, course 1)
unsigned long BFS(adjlist* g, unsigned long s, bool *Mark, unsigned long *pointeur_nb_nodes_largest_cc){  // BFS from node s with a graph defined as an adjacency array
	queue *q = init_queue(s);
	
	Mark[s] = 1;
	unsigned long nb_of_nodes = 0; // nb of nodes in the connected component starting from the node s
	unsigned long count = 0;
	while (q->first != NULL) {  // while FIFO not empty
		nb_of_nodes += 1;
		unsigned long u = pop_queue(q);

		count +=1;
		if (count % 1000 == 0) {
			printf("node %lu / %lu\n", count, g->n);
			printf("degree : %lu \n", g->cd[u+1] - g->cd[u]);
		}

		// printf("%lu ", u); // output u
		unsigned long degree = g->cd[u+1] - g->cd[u]; // degree of node u

		for (int i=0;i<degree;i++) { // for each neighbor v of u in the graph g
			unsigned long v = g->adj[ g->cd[u] + i ];
			if (Mark[v] == 0) { // v is not marked
				add_queue(q, v); // add v to the queue
				Mark[v] = 1; // v is marked
			}
		}
	}

	*pointeur_nb_nodes_largest_cc = max(*pointeur_nb_nodes_largest_cc, nb_of_nodes);
	// free(Mark);	
	free_queue(q);

	return nb_of_nodes;

}


void connected_components(adjlist *g, bool print) { 
// outputs all connected componenents and their sizes (uncomment the print)
// and the fraction of nodes in the largest connected component
	bool *Mark = calloc(g->n,sizeof(bool)); // Mark[i] = 1 if node i has been visited, O o/w 

	int nb_of_cc = 1; // nb of connected components
	unsigned long nb_nodes_largest_cc = 0; 
	unsigned long nb_of_nodes = 0;

	for (unsigned long s=0;s<g->n;s++) { // loop over all nodes
		if (Mark[s] == 1) { // node already marked
			continue;
		}

		nb_of_nodes = BFS(g, s, Mark, &nb_nodes_largest_cc);
		if(print == 1) {
			printf("\n");
			printf("Connected componenent n° %d size : %lu", nb_of_cc, nb_of_nodes);
		}


		nb_of_cc += 1;
	}
	printf("\nThe number of nodes in the largest connected component is %lu\n", nb_nodes_largest_cc);
	printf("The fraction of nodes in the largest connected component is %.3f\n\n", (float)(nb_nodes_largest_cc)/(g->n));

	free(Mark);
}


typedef struct {
	unsigned long node;
	int distance;
} res_bfs;

res_bfs BFS_procedure(adjlist *g, unsigned long s) { 
// outputs the node with maximal distance from source s (in the cc) and the distance
	queue *q = init_queue(s);
	bool *Mark = calloc(g->n,sizeof(bool)); // Mark[i] = 1 if node i has been visited, O o/w 
	int *dist = calloc(g->n, sizeof(int)); // stores the distances from source s
	res_bfs res; // result
	res.node = s;
	res.distance = 0;
	
	Mark[s] = 1;

	while (q->first != NULL) {  // while FIFO not empty
		unsigned long u = pop_queue(q);
		// printf("%lu ", u); // output u
		unsigned long degree = g->cd[u+1] - g->cd[u]; // degree of node u

		for (int i=0;i<degree;i++) { // for each neighbor v of u in the graph g
			unsigned long v = g->adj[ g->cd[u] + i ];
			if (Mark[v] == 0) { // v is not marked
				dist[v] = dist[u] + 1;
				if (dist[v] > res.distance) {
					res.node = v;
					res.distance = dist[v];
				}
				add_queue(q, v); // add v to the queue
				Mark[v] = 1; // v is marked
			}
		}
	}

	free(Mark);	
	free(dist);
	free_queue(q);

	return res;
}


// Lower bound diameter algorithm: // source : slide 34 of https://www.irif.fr/~habib//Documents/Diametre.pdf
// (Randomized BFS procedure)

// Repeat alpha times:
// Randomly choose a vertex w in V
// u <- BFS(w)
// v <- BFS(u)

// Select the vertices u_O, v_O s.t. dist(u_O,v_O) is maximal

int lower_bound(adjlist *g, int alpha) { 
// computes a lower bound to the graph diameter using a randomizd BFS procedure
// alpha is the number of BFS procedures
	int bound = 0;

	for (int i=0;i<alpha;i++) {
		unsigned long w = rand() % (g->n); // randomly choose a node w
		// printf("w = %lu\n", w);
		unsigned long u = BFS_procedure(g, w).node; // node u with maximal distance from source w
		// printf("u = %lu\n", u);
		int bound_test = BFS_procedure(g, u).distance; // maximal distance from source u
		// printf("bound test = %d\n", bound_test);
		bound = max(bound, bound_test);
	}
	return bound;
}



static int
int_compare(const void *p1, const void *p2) {
 unsigned long left =	*(const	unsigned long *)p1;
 unsigned long right = *(const unsigned long *)p2;

 return ((left > right) - (left < right));
}

void *re_index(adjlist *g, unsigned long *new_label) { // sort the nodes in non-increasing order of degree and re-index the graph

	unsigned long *degrees = malloc(g->n * sizeof(unsigned long));
	unsigned long *degrees_sorted = malloc(g->n * sizeof(unsigned long));

	bool *mark = calloc(g->n,sizeof(bool));
	for (unsigned long u=0;u<g->n;u++) {
		degrees_sorted[u] = g->cd[u+1] - g->cd[u];
		degrees[u] = g->cd[u+1] - g->cd[u];
	}

	qsort(degrees_sorted, g->n, sizeof(unsigned long), int_compare); // sort the degree list

	unsigned long i = 0;
	unsigned long index;
	for (index=(g->n)-1;index>0;index--) { // The higher the index the lower the degree
		for (unsigned long j=0;j<g->n;j++) { // Loop over all nodes to find one with degree = degrees_sorted[i]
			if (mark[j] == 0 && (degrees[j] == degrees_sorted[i])) {
				mark[j] = 1;
				new_label[j] = index; // new label to the node j
				i+=1;
				break;
			}
		}
	}
	index=0; // case index = 0;
	for (unsigned long j=0;j<g->n;j++) { // Loop over all nodes to find one with degree = degrees_sorted[i]
		if (mark[j] == 0 && (degrees[j] == degrees_sorted[i])) {
			mark[j] = 1;
			new_label[j] = index; // new label to the node j
			i+=1;
			break;
		}
	}	
	

	free(degrees);
	free(degrees_sorted);
	free(mark);

}






unsigned long triangles(adjlist *g) { // outputs the number of triangles

	unsigned long nb_of_triangles = 0;
	unsigned long **tsl = malloc(g->n * sizeof(unsigned long *)); // truncated and sorted list of neighbors of each node
	unsigned long *len_tsl = malloc(g->n * sizeof(unsigned long));
	unsigned long *new_label = malloc(g->n * sizeof(unsigned long));

	re_index(g, new_label); // sort the nodes in non-increasing order of degree and re-index the graph

	for (unsigned long u=0;u<g->n;u++) { // for each node u
		unsigned long degree = g->cd[u+1] - g->cd[u]; // degree of u
		tsl[u] = (unsigned long *)malloc(degree*sizeof(unsigned long));
		int index = 0;
		for (int i=0;i<degree;i++) { // for each neighbor v of u in the graph g
			unsigned long v = g->adj[ g->cd[u] + i ];
			unsigned long degree_neighbor = g->cd[v+1] - g->cd[v];
			// printf("node %lu : ", v)
			if (new_label[u] < new_label[v]) { // we consider only the neighbors with smaller degree
				tsl[u][index] = new_label[v];
				index+=1;
			}
	
		}
		len_tsl[u] = index;

		qsort(tsl[u], index, sizeof(tsl[u][0]), int_compare);
	}

	for (int k=0;k<g->e;k++) { // for each edge (u,v) in G, we look for the intersection of tsl[u] and tsl[v]
		unsigned long u=g->edges[k].s;
		unsigned long v=g->edges[k].t;
		int i = 0;
		int j = 0;

		while (i < len_tsl[u] && j < len_tsl[v]) {
			if (tsl[u][i]==tsl[v][j]) { // node in the intersection of tsl[u] and tsl[v]
				nb_of_triangles += 1;
				i+=1;
				j+=1;
			}
			else if (tsl[u][i] < tsl[v][j]) {
				i+=1;
			}
			else {
				j+=1;
			}
		}
	}

	for (unsigned long u=0;u<g->n;u++) { 
		free(tsl[u]);
	}
	free(len_tsl);
	free(new_label);

	return nb_of_triangles;
}







int main(int argc,char** argv){
	adjlist* g;
	time_t t1,t2,t3,t4;

	t1=time(NULL);

	printf("Reading edgelist from file %s\n",argv[1]);
	g=readedgelist(argv[1]);

	printf("Number of nodes: %lu\n",g->n);
	printf("Number of edges: %lu\n",g->e);

	printf("Building the adjacency list\n");
	mkadjlist(g);




	// Exercise 5 - BFS

	char *ex = argv[2];

	if(strcmp(ex, "cc")==0) { // run "./tp1 graph.txt cc print" with print = O or 1 
		printf("Running the connected components algorithm\n");
		char *print = argv[3];
		if(strcmp(print, "1")==0) {
			printf("ok\n");
			connected_components(g, 1); // outputs the sizes of each connected components and the fraction of nodes in the largest connected component 	
		}
		else {
			connected_components(g, 0); // output only the fraction of nodes in the largest connected component 
		}
	}

	else if (strcmp(ex, "lower_bound")==0) { // run "./tp1 graph.txt lower_bound"
		t2=time(NULL);
		printf("Running the lower bound algorithm\n");
		printf("lower bound : %d with alpha = %d\n", lower_bound(g, 5), 5); // email (alpha = 5), amazon (5), lj (), orkut (), friendster ()
		t3=time(NULL);
		printf("- Running time lower bound algorithm = %ldh%ldm%lds\n",(t3-t2)/3600,((t3-t2)%3600)/60,((t3-t2)%60));
	}




	// Exercise 6 - Triangles

	else if (strcmp(ex, "triangle")==0) { // run "./tp1 graph.txt triangle"
		t2=time(NULL);
		printf("Running the triangle algorithm\n");
		printf("nb of triangles : %lu\n", triangles(g));
		t3=time(NULL);
		printf("- Running time triangle algorithm = %ldh%ldm%lds\n",(t3-t2)/3600,((t3-t2)%3600)/60,((t3-t2)%60));
	}





	free_adjlist(g);

	t4=time(NULL);

	printf("- Overall time = %ldh%ldm%lds\n",(t4-t1)/3600,((t4-t1)%3600)/60,((t4-t1)%60));
	

	return 0;
}
