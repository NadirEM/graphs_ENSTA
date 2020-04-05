#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#define NLINKS 100000000

unsigned long max3(unsigned long a,unsigned long b,unsigned long c){
	a=(a>b) ? a : b;
	return (a>c) ? a : c;
}


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

//building the adjacency matri

//freeing memory
void free_adjlist(adjlist *g){
	free(g->edges);
	free(g->cd);
	free(g->adj);
	free(g);
}

void decrease_degree(adjlist *g, unsigned long u, unsigned long *degree){
		unsigned long deg = g->cd[u+1] - g->cd[u]; // degree of node u

		for (int i=0;i<deg;i++) { // for each neighbor v of u in the graph g
			unsigned long v = g->adj[ g->cd[u] + i ];
			degree[v] -= 1;
}
}


int nb_nodes(char* file_name)
{
    FILE* file = fopen (file_name, "r");
    unsigned long i = 0;
    unsigned long j =0;
    unsigned long n = 0;
    unsigned long e = 0;
    while (fscanf(file,"%lu %lu", &i, &j)==2)
    {
        n=max3(n,i,j);
        e++;
    }
    fclose (file);
    n++;

    return n;
}

unsigned long* nb_degrees(char*input){
    unsigned long n = nb_nodes(input);
    unsigned long* nb_degrees = (unsigned long*)calloc(n, sizeof(unsigned long));
    FILE* file = fopen (input, "r");
    unsigned long s =0;
    unsigned long t=0;

    while (fscanf(file,"%lu %lu", &s, &t)==2){
            nb_degrees[s] ++;
            nb_degrees[t]++;

    }
    fclose(file);
    return(nb_degrees);

};


unsigned long argmin(unsigned long* liste, unsigned long len){

    unsigned long min =liste[0];
    unsigned long arg = 0;
    unsigned long i=0;

    for (i=1;i<len;i++){

        if (liste[i]<min)
    {  ;
        min = liste[i];
        arg = i;
    }

}
 return(arg);

}

void mkadjlist(adjlist* g){
	unsigned long i,u,v;
	unsigned long *d=calloc(g->n,sizeof(unsigned long)); // Degree array

	for (i=0;i<g->e;i++) { // build the degree array
		d[g->edges[i].s]++;

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

	}

	free(d);
	//free(g->edges);
}



unsigned long maximum(unsigned long* liste, unsigned long len){

    unsigned long max =liste[0];

    unsigned long i=0;

    for (i=1;i<len;i++){

        if (liste[i]>max)
    {
        max = liste[i];

    }

}
 return(max);

}

unsigned long* k_decomp(char* input){
    unsigned long n = nb_nodes(input);
    unsigned long* cores = (unsigned long*)calloc(n, sizeof(unsigned long));
    unsigned long* degrees = nb_degrees(input);

    unsigned long v;
    unsigned long s = 0;
    unsigned long t = 0;
    v = argmin(degrees, n);
    unsigned long c = 0;
    if (c<degrees[v]){
        c= degrees[v];
    }
    FILE* file;
    adjlist* g;
    g=readedgelist(input);
    mkadjlist(g);

    while (degrees[v]< n){

        cores[v] = c;

    if (degrees[v]>0){

        decrease_degree(g, v, degrees);

    }



    degrees[v] = 2*n+1;
    v = argmin(degrees, n);
    if (c<degrees[v]){
        c= degrees[v];
    }





    }
    return(cores);

}


int main(){

	time_t t1,t2;
	t1=time(NULL);
	char* file= "net.txt";
	unsigned long n = nb_nodes(file);
	unsigned long* degrees = nb_degrees(file);

	FILE* out;
	out  = fopen("c_degrees.csv", "w");
	int i = 0;
	for(i=0;i<n;i++){

        fprintf(out, "%d %d\n", i, degrees[i]);
    }

    fclose(out);



	unsigned long* cores = k_decomp(file);

    unsigned long max = maximum(cores,n);
    printf("graph core %lu", max);
    out= open("cores.csv", "w");
    i=0;
    for(i=0;i<n;i++){

        fprintf(out,"%d %d", i, cores[i]);
    }
    fclose(out);



   // unsigned long* agmins = argmins(file, k, pr);


    t2=time(NULL);
	printf("- Overall time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));

	return 0;
}



