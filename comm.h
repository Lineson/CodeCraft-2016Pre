
/*
 *  global var
 */
int data[700][700]; // raw topo marked by edge id
int w[5000]; // weight of each edge
int g[700][700]; // graph topo to work on
int adj[700][10]; // adjacency list -- adj[i][0] : number of neighbors
int prt[700][20]; // parent list -- prt[i][0] : number of parents
bool spec_2[700]; // specified vertexes
int spec[700];
int max_vs; // numbers of specified vertexes
int max_v; // number of vertexes
int max_e; // number of edges
int s; // vertex s
int t; // vertex t
int src; // current source
int dst; // current destination
int dist_to_dst[700]; // shortest distance from vertex i to vertex t

// Dijkstra data type
bool visited[700]; // visit hash table
int dist[700]; // shortest path distance
int heap[700]; // shortest path heap -- start from heap[1]
int heap_size; // current heap size
int loc[700]; // location of each vertex in heap
int path[700]; // path -- path[0] : length of path
int path_each[610][610][200];
int pre[700]; // path precursor list
bool xx[700]; // non-visitable vertexes

// DFS data type
int f[700][700]; // f[i][j] = 0 : g[i][j] unvisited; f[i][j] = 1 : visited
int a[700][700]; // a[i][j] = 0 : g[i][j] not available; a[i][j] = 1 : available
int r[700]; // r[i] = sum of a[i][j]
int c[700]; // c[i] = sum of a[j][i]
int row[700]; // row[i] = sum of f[i][j]
int col[700]; // col[i] = sum of f[j][i]
int best; // current best cost
int vs_cnt; // number of specified vertexes visited
int stack[20000];
int top; // stack top;
int max_dep;


//test
int RepairValid = 1;
int last_kt = 0;


#define max_cost 99999
#define PRANDMAX 1000000000
static int aNum = 0;
static int bNum = 0;
static int arr[55];
int Rand();
void initRand (int seed)
{
    int i, ii;
    int last, next;

    seed %= PRANDMAX;
    if (seed < 0) seed += PRANDMAX;

    arr[0] = last = seed;
    next = 1;
    for (i = 1; i < 55; i++) {
        ii = (21 * i) % 55;
        arr[ii] = next;
        next = last - next;
        if (next < 0)
            next += PRANDMAX;
        last = arr[ii];
    }
    aNum = 0;
    bNum = 24;
    for (i = 0; i < 165; i++)
        last = Rand ();
}

int Rand (void)
{
    int tp;

    if (aNum-- == 0)
        aNum = 54;
    if (bNum-- == 0)
        bNum = 54;
    if (aNum > 54 || bNum > 54)
    {
        aNum = aNum%55;
        bNum = bNum%55;
    }
    tp = arr[aNum] - arr[bNum];

    if (tp < 0)
        tp += PRANDMAX;

    arr[aNum] = tp;

    return tp;
}

#define RREAL ((double)Rand()/PRANDMAX)
#define RANDOM Rand
#define unifRand(n) (rand()%n)

#define INF 100000  // ZP 3-22
#define PARA_MAX 100000000
/*Accept probability of replacing node in [s,t] with [t+1,n]*/
#define PARA_01 ((double)Rand()/PARA_MAX)  
/*Accept probability of exchange between [s,t]*/
#define PARA_02 ((double)Rand()/PARA_MAX)


#define MOD(i,n)    ((i) % (n) >= 0 ? (i) % (n) : (i) % (n) + (n))
#define DirctD(x,y) DirectDist[(x)][(y)]    //get the distence of x,y BY Matrix[x][y]  

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))
#define sqr(x)   ((x)*(x))