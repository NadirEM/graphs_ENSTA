
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#define NLINKS 100000000

typedef struct {
	unsigned long s;
	unsigned long t;
} edge;

typedef struct {
	unsigned long n;
	unsigned long e;
	edge *edges;
} edgelist;

unsigned long max3(unsigned long a,unsigned long b,unsigned long c){
	a=(a>b) ? a : b;
	return (a>c) ? a : c;
}


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


unsigned long nb_nodes(char* file_name)
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

int* argmins(char* file, int k, double* vec){
    unsigned long int n = nb_nodes(file);
    int* argmins = malloc(k*sizeof(int));
    int i =0;
    int j = 0;
    argmins[0] = 0;
    for (i=1;i<k;i++){
        j =0;
        argmins[i] = i;
        while (j<i){

            if (vec[argmins[j]] > vec[i]){
                argmins[j] = i;

                argmins[i] = j;
            }
            j =j+1;
            }
            }

    int m = k-1;
    while(m<n+1){

                m++;
                i = m;




                int j = 0;
                int r;
                int stop = 0;
                if (vec[i]<vec[argmins[j]]){
                        stop =1;
                    }
                    if (j<k){
                        stop=1;
                    }


                while (!stop)
                {

                    j ++;
                    if (vec[i]<vec[argmins[j]]){
                        stop =1;
                    }
                    if (j<k){
                        stop=1;
                    }

                };

                if (vec[i]<vec[argmins[j]]) {
                        printf("modification %d ", i);
                        int a = argmins[j];
                        int b;
                        for (r = j; r<k-1;r++){

                            b = argmins[r+1];
                            argmins[r+1] = a;
                            a = b;};
                        argmins[j] = i;

                };

    };

    int g =0;
    for(g=0;g<k;g++){
        printf("args %d", argmins[g]);
    }

    return argmins;}

unsigned long int* nb_degrees(char*input){
    int n = nb_nodes(input);
    unsigned long int* nb_degrees = (int*)calloc(n, sizeof(int));
    FILE* file = fopen (input, "r");
    int s =0;
    int t=0;

    while (fscanf(file,"%lu %lu", &s, &t)==2){
            nb_degrees[s] ++;
            nb_degrees[t]++;
    }
    fclose(file);
    return(nb_degrees);

}


int* degree_out(char* input)
{
    unsigned long n = nb_nodes(input);
    int* nodes = malloc(n*sizeof(int));
    int k=0;
    for (k=0; k<n ; k++)
    {
        nodes[k]=0;
    }
    FILE* file = fopen (input, "r");
    unsigned long i = 0;
    unsigned long j =0 ;
    while (fscanf(file,"%lu %lu", &i, &j)==2)
    {
        nodes[i] = nodes[i] + 1;
    }
    fclose (file);
    return nodes ;
}

int* degree_in(char* input)
{
    unsigned long n = nb_nodes(input);
    int* nodes = malloc(n*sizeof(int));
    int k=0;
    for (k=0; k<n ; k++)
    {
        nodes[k]=0;
    }
    FILE* file = fopen (input, "r");
    unsigned long i = 0;
    unsigned long j =0 ;
    while (fscanf(file,"%lu %lu", &i, &j)==2)
    {
        nodes[j] = nodes[j] + 1;
    }
    fclose (file);
    return nodes ;
}


void PR_pw_it(char *input, int* n_iter, double*P_k, double* err, double alpha){



    unsigned long n = nb_nodes(input);
    int* degrees = degree_out(input);

    double* P = malloc(n*sizeof(double));
    //vectors initialisation

    int i=0;
    for (i=0;i<n;i++){
            P[i] = 0;
            P_k[i] = (double)(1/(double)(n));
   }

    unsigned long s = 0;
    unsigned long t =0;
    int k=0;
    double error =  0.0;
    int r = 0;
    int m =0;
    double sum = 0;
    while (k<n_iter){


        for (i=0;i<n;i++){
            P[i] = (double)(0);

                }
    FILE* file = fopen (input, "r");
    adjlist* g;
    g=readedgelist(input);
    mkadjlist(g);


    unsigned long  node =0;

        for (node=0;node<n;node++){
            int deg = degrees[node];


            if (deg == 0){

                    int ver = 0;
                    for(ver=0;ver<n; ver++){
                        P[ver] = P[ver] + P_k[node]/((double)n);
                    }

            }
            else{
            for (int i=0;i<deg;i++) { // for each neighbor v of u in the graph g
			unsigned long v = g->adj[ g->cd[node] + i ];

                 P[v] = P[v] + P_k[node]/((double)deg);
}

                s = 0;
                t =0;



            }
        }
        sum = 0.0;

        int l;
        for (l=0;l<n;l++){


                P[l] = (double)(1-alpha)*P[l] + alpha/((double)n);

                sum = sum + fabs(P[l]);

        }
        m =0;
        for (m=0;m<n;m++)
        {
                P[m] = P[m] + (1- sum)/((double)n);
        }
        r = 0;

        error = 0.0;

        for(r=0; r<n; r++)

        {  error += fabs((double)P[r] - (double)P_k[r]);

        }





        r =0;
        for (r=0;r<n;r++){

            P_k[r] = P[r]; // memorizong the new value of the vector P

        }

        err[k] =error;
        k += 1;


    }

}

void free_edgelist(edgelist *g){
	free(g->edges);
	free(g);
}


int main(){
	edgelist* g;
	time_t t1,t2;
	t1=time(NULL);
	char* file= "soc-Epinions.txt";
	//g=readedgelist(file);
	unsigned long int n = nb_nodes(file);

    int* degrees_in = degree_in(file);
    int* degrees_out = degree_out(file);
    int k = 0;
    int n_iters = 15;
    double alpha = 0.15;
    double* pr = malloc(n*sizeof(double));
    double* errors = malloc(n*sizeof(double));



        printf("run for alpha= %f ...", alpha);
        PR_pw_it(file, n_iters, pr, errors, alpha);
       // unsigned long int* degrees = nb_degrees(file);
        printf("exporting the result alpha = %f ", alpha);

        FILE* out = fopen("results_015.csv", "w");
        int a=0;
        for (a=0;a<n;a++){
          fprintf(out, "%d %lf \n", a , pr[a]);
        }
       fclose(out);
       out = fopen("errors_15.csv", "w");
        a =  0;
        for (a=0;a<n_iters;a++){
          fprintf(out, "%d %lf \n", a , errors[a]);
       }
        fclose(out);



	printf("- Overall time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));

	return 0;
}



