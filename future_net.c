/*README
  # About
  The main algorithm is SAO(Simulated Annealing Optimization) and ACS(Ant Colony System). 
  All the codes are written by our group members, and no open source libs are used. 
  
  
  # Compile 
  For the best performance, you can complile with the following shell script by gcc4.8.4:
  gcc -O2 -std=c99 future_net.c -o future_net -lm
  
  
  # Authors 
  We come from department of Computer Science and Technology, Xi'an Jiaotong University. 
  郑鹏    zeepean@gmail.com
  卞正达  kurisusnowdeng@gmail.com
  王瑞龙  wangruilong0212@gmail.com
*/
#include <string.h>    /* memcpy */
#include <math.h>      /* exp    */
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>

#include <time.h>
#include <sys/timeb.h>
#include <sys/time.h>

#define DEBUG 1
#ifdef DEBUG
    #define DBG printf
#else
    #define DBG(format, arg...) do { ; } while (0)
#endif



//--------------------------------ACO below----------------------------------//

#include "malloc.h"
#include "math.h"


int  Agent_number  =  30; //30
int    Opt_time    =  100;
double UN_REACH    =  10e9;
double Phe_keep    =  0.90;
double Phe_fac     =  1.0;     //信息素因子，信息素的重要程度
double Inspir_fac  =  2.5;    //启发因子
double Toatal_Phe  =  100.0;
double Initial_Phe =  0.9;

int slength ;
int a_feasible_solution[700];

/* 定义个体信息 */  

typedef struct agent 
{
    int path[700] ;
    bool id_used[700] ;
    int VS_counter;
    int last_pos;
    double path_len;
}Agent;

bool vs[700] ;              // specified vertexes
int  spec_vex[700];
bool total_v[700];
// int  max_vs ;                 // numbers of specified vertexes
// int  s ;                      // source
// int  t ;                      // destination
int  city_num;
double length_table[700][700];
// int adj[700][10];               // adjacency list -- adj[i][0] : number of neighbors
bool global_id_used[700];
// int data[700][700];
int  spec_To_Change[700];
int  spec_Chan_num;

//ACS variable
double Deci_p;
double Phe_Matrix[700][700];
Agent Agents[40];
Agent Bst_Agent;
char *result;

//--------------------------------ACO above----------------------------------//



//-------------------------------global var----------------------------------//

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

double AS_T_INIT              = 100;
double AS_FINAL_T             = 0.1;
double AS_COOLING             = 0.99; /* to lower down T (< 1) */

int    AS_no_change_time      = 200;
int    AS_TRIES_PER_T         = 400;
int    AS_IMPROVED_PATH_PER_T = 5;

long   AS_seed                = -314159L;
bool   AS_small_case          = false;

/*
 * Defs-------------------------------------------------------------------------------
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

#define max_cost 99999

//test
int RepairValid = 1;
int last_kt = 0;

//-------------------------global var----------------------------------------//







/*****************************************************************************/ 
/**********************************ACS****************************************/ 
//返回指定范围内的随机浮点数
double rnd(double dbLow, double dbUpper)
{
    double dbTemp = rand() / ((double)RAND_MAX + 1.0);
    return dbLow + dbTemp*(dbUpper - dbLow);
}

//返回浮点数四舍五入取整后的浮点数
double ROUND(double dbA)
{
    return (double)((int)(dbA + 0.5));
}


void ACS_read_topo(const char* topo) {

    memset(length_table,1,sizeof(length_table));
    memset(data, 0xff, sizeof(data));
    memset(adj, 0, sizeof(adj));

    for (int i = 0; i < 700; ++i)
        for (int j = 0; j < 700; ++j)
            length_table[i][j]=UN_REACH;

    memset(total_v,0,sizeof(total_v));
    FILE *input = fopen(topo, "r");
    int id;
    int x; //source
    int y; //destination
    int weight;
    city_num = 0;
    //max_e = 0;
    while (!feof(input)) {
        if (fscanf(input, "%d,%d,%d,%d", &id, &x, &y, &weight) > 0) {
            if (length_table[x][y]==UN_REACH)
            {
                adj[x][0]++;
                adj[x][adj[x][0]] = y;

            }
            if (length_table[x][y]>weight)
            {
                length_table[x][y] = (double)weight;
                data[x][y] = id;
            }
            //length_table[x][y] = (length_table[x][y]<weight)?length_table[x][y]:weight;
            if (total_v[x]==0)
            {
                city_num+=1;
            }
            if (total_v[y]==0)
            {
                city_num+=1;
            }
            total_v[x]=1;
            total_v[y]=1;
        }
    }
    fclose(input);
}


void ACS_read_demand(const char* demand) {
    memset(vs, 0, sizeof(vs));
    max_vs = 0;
    int x;
    FILE *input = fopen(demand, "r");
    fscanf(input, "%d,%d,", &s, &t);
    while (!feof(input)) {
        if (fscanf(input, "%d|", &x) > 0) {
            vs[x] = 1;
            spec_vex[max_vs] = x;
            max_vs++;
        }
    }

    fclose(input);
}

void Change_spec()
{
    int count = 0 , i = 0;
    int node;
    spec_Chan_num = 0;
    int temp_spec_vex[700];
    memcpy(temp_spec_vex , spec_vex , sizeof(spec_vex));

    while(1)
    {
        count = 0;
        for (i = 0; i < max_vs; ++i)
            if (adj[temp_spec_vex[i]][0] == 1)
            {
                count ++ ;
                spec_Chan_num += 1;
                spec_To_Change[spec_Chan_num*2-2] = temp_spec_vex[i];
                spec_To_Change[spec_Chan_num*2-1] = adj[temp_spec_vex[i]][1];
                node = adj[temp_spec_vex[i]][1];
                temp_spec_vex[i] = node;
            }
        if (count == 0)
           return;
    }
}


void Change_adj()
{
    int i , j , k;
    int node , node_next;
    bool tag;
    for (i = 0; i < spec_Chan_num*2; i+=2)
    {
        node = spec_To_Change[i];
        node_next = spec_To_Change[i+1];
        for(j = 0; j < city_num ; ++j)
        {
            if( j == node)
                continue;
            tag = 0;
            for (k = 0; k <= adj[j][0]; ++k)
                if (adj[j][k] == node_next)
                {
                    length_table[j][adj[j][k]] = UN_REACH;
                    tag = 1;
                }
            // if(tag)
            //     adj[j][0]--; 
        }
    }
} 

/*****************************************************************************/ 
/*                                                                           */ 
/*                      Ant Conlony System                                   */ 
/*                                                                           */  
/*****************************************************************************/ 
void Initial_Phe_Matrix()
/*初始化信息素矩阵*/
{
    for (int i = 0; i < city_num; ++i)
        for (int j = 0; j < city_num; ++j)
            Phe_Matrix[i][j] = Initial_Phe;
}

void Initial_Agents()
{
    for (int i = 0; i < Agent_number; ++i)
    {        
        Agents[i].path[0] = s;

        Agents[i].VS_counter = 0;

        Agents[i].last_pos = 0;

        memset(Agents[i].id_used,0,sizeof(Agents[i].id_used));

        Agents[i].id_used[s] = 1;

        for (int j = 0; j < city_num; ++j)
            if (adj[j][0] <= 0 && j!=t)
                Agents[i].id_used[j] = 1;

        Agents[i].path_len = 0.0;
    }
}
int counter = 0;


int Chose_Nex_Node(Agent *Ag)
{
    // printf("%d|", Ag.path[Ag.last_pos]);
    int i = Ag->path[Ag->last_pos] ;

    int next_city = -1;

    double temp_q = 0.0;

    int k ;

    double total_valu = 0.0; //所有与i相关的邻接点的信息素与启发式信息乘积

    double temp_tol ;

    double chose_pro[700];

    memset(chose_pro,0,sizeof(chose_pro));

    /*计算各邻接点的概率*/
    for (k = 1; k <= adj[i][0]; ++k)
    {

        if (!Ag->id_used[adj[i][k]] && length_table[i][adj[i][k]]<UN_REACH)
        {
            
            if (vs[adj[i][k]])
            {
                chose_pro[adj[i][k]] = 0.8;
                
            }
            else
            {
                chose_pro[adj[i][k]] = pow(Phe_Matrix[i][adj[i][k]], Phe_fac)*pow(1.0/length_table[i][adj[i][k]], Inspir_fac);
                
            }
            if (Ag->VS_counter != max_vs && adj[i][k] == t)
            
                chose_pro[adj[i][k]] = 0.0;

            total_valu += chose_pro[adj[i][k]] ; 
        }
        else

            chose_pro[adj[i][k]] = 0.0;
    }

   
    if (total_valu > 0.0)//所谓轮盘选择，抄了，需修改
    {
        temp_tol = rnd(0.0,total_valu);
    
        for (k = 0; k < city_num; ++k)
        {
            if (!Ag->id_used[k]&& length_table[i][k]<UN_REACH)
            {
                temp_tol -= chose_pro[k];
               
                if (temp_tol < 0.0)
                {
                    next_city = k;
                    break;
                }
            }
        }//end for
    }//end if(total_valu)

    if (next_city == -1)
    {
        for (k = 1; k <= adj[i][0]; ++k)
            if (!Ag->id_used[adj[i][k]] && length_table[i][adj[i][k]]<UN_REACH)
            {
                next_city = adj[i][k];
                break;
            }
    }

    return next_city;
}


bool add_city(Agent *Ag)
{
    int next_city = Chose_Nex_Node(Ag);

    if (next_city == -1)
    {
        return 0;
    }
    Ag->last_pos += 1;
    Ag->path[Ag->last_pos] = next_city;
    Ag->id_used[next_city] = 1;
    Ag->path_len += length_table[Ag->path[Ag->last_pos-1]][next_city];
    return 1;
}


void Naive_Solution(Agent *Ag)
{

    int i , j;

    bool tag;
    Ag->VS_counter = 0;
    while(Ag->path[Ag->last_pos] != t || Ag->VS_counter != max_vs)
    {
        tag = add_city(Ag);
        if(!tag)
        {   

            if (vs[Ag->path[Ag->last_pos]])
            {
                Ag->path_len = UN_REACH;
                counter += 1;
                
                return;
            }
            Ag->last_pos -= 1;
            Ag->path_len -= length_table[Ag->path[Ag->last_pos]][Ag->path[Ag->last_pos+1]];
            if (vs[Ag->path[Ag->last_pos]])
            {
                Ag->VS_counter--;
            }
            counter += 1;
            
            continue;

        }//end if(!tag)
        if (vs[Ag->path[Ag->last_pos]])
        {
            Ag->VS_counter += 1;
        }
        
    }// end while

    if (Ag->path[Ag->last_pos] == t &&  Ag->path_len<Bst_Agent.path_len && Ag->VS_counter == max_vs)
    {
        memcpy(&Bst_Agent,Ag,sizeof(Bst_Agent));
    }
 
}

void BestWay_Update_Phe()
/*   利用当前Agent所计算最短路，全局Update 信息素矩阵  */ 
{
    
    int i , j ;

    double incre_Phe[700][700];

    memset(incre_Phe,0,sizeof(incre_Phe));
    for (i = 0; i < city_num; ++i)
        for (j = 0; j < city_num; ++j)
            incre_Phe[i][j] = 0.0;

    if (Bst_Agent.path_len == UN_REACH)
    
        for (i = 0; i < Agent_number; ++i)

            for (j = 0; j <= Agents[i].last_pos-1; ++j)
        
                incre_Phe [Agents[i].path[j]] [Agents[i].path[j+1]] = incre_Phe [Agents[i].path[j]] [Agents[i].path[j+1]] + Toatal_Phe/Agents[i].path_len;
    else

        for (j = 0; j <= Bst_Agent.last_pos-1; ++j)
        
            incre_Phe [Bst_Agent.path[j]] [Bst_Agent.path[j+1]] = incre_Phe [Bst_Agent.path[j]] [Bst_Agent.path[j+1]] + Toatal_Phe/Bst_Agent.path_len;

    for (i = 0; i < city_num; ++i)

        for (j = 0; j < city_num; ++j) 
        
            Phe_Matrix[i][j] = Phe_keep*Phe_Matrix[i][j] + incre_Phe[i][j]; 
}


void ACS_FOR_SearchWay()
{

    Initial_Phe_Matrix();

    int Opt_times; 

    for (Opt_times = 0; Opt_times < Opt_time; ++Opt_times)
    {

        Initial_Agents();

        for (int i = 0; i < Agent_number; ++i)
        
            Naive_Solution(&Agents[i]);
        


        BestWay_Update_Phe();

        if (Opt_times == Opt_time-1)
        {
            slength = Bst_Agent.last_pos+1;
            int pos = slength;
            memset(Bst_Agent.id_used,0,sizeof(Bst_Agent.id_used));
            for (int i = 0; i <= Bst_Agent.last_pos; ++i)
            {
                Bst_Agent.id_used[Bst_Agent.path[i]] = 1;
            }
            memcpy(a_feasible_solution,Bst_Agent.path,sizeof(a_feasible_solution));
            for (int i = 0; i < city_num; ++i)
            {
                if (!Bst_Agent.id_used[i])
                {
                    a_feasible_solution[pos] = i;
                    pos += 1;
                }
            }
            return ;
        }
    }

}
/********************************ACS******************************************/ 
/*****************************************************************************/ 




void up(int x) {
    int i = loc[x];
    while (i > 1 && dist[heap[i]] < dist[heap[i / 2]]) {
        heap[i] = heap[i] ^ heap[i / 2];
        heap[i / 2] = heap[i] ^ heap[i / 2];
        heap[i] = heap[i] ^ heap[i / 2];
        loc[heap[i]] = i;
        loc[heap[i / 2]] = i / 2;
        i = i / 2;
    }
}

void down(int x) {
    int i = loc[x];
    int k;
    while ((i + i <= heap_size && dist[heap[i]] > dist[heap[i + i]]) ||
        (i + i + 1 <= heap_size && dist[heap[i]] > dist[heap[i + i + 1]])) {
        k = i + i;
        if (k + 1 <= heap_size && dist[heap[k]] > dist[heap[k + 1]]) {
            k++;
        }
        heap[i] = heap[i] ^ heap[k];
        heap[k] = heap[i] ^ heap[k];
        heap[i] = heap[i] ^ heap[k];
        loc[heap[i]] = i;
        loc[heap[k]] = k;
        i = k;
    }
}

void push(int x) {
    heap_size++;
    heap[heap_size] = x;
    loc[x] = heap_size;
    up(x);
}

int pop() {
    int ans = heap[1];
    heap[1] = heap[heap_size];
    loc[heap[1]] = 1;
    heap_size--;
    down(heap[1]);
    return ans;
}


void get_path(int v) {
    if (pre[v] >= 0 && v != src) {
        get_path(pre[v]);
    }
    else {
        return;
    }
    path_each[src][dst][0]++;
    path_each[src][dst][path_each[src][dst][0]] = v;
}

void find_path(int begin, int end) {
    // Dijkstra
    memset(dist, 0x7f, sizeof(dist));
    memset(pre, 0xff, sizeof(pre));
    memset(heap, 0, sizeof(heap));
    memset(visited, 0, sizeof(visited));
    memset(loc, 0, sizeof(loc));
    visited[begin] = true;
    int u, v;
    heap_size = 0;
    int len = adj[begin][0];
    for (int i = 1; i <= len; i++) {
        v = adj[begin][i];
        dist[v] = g[begin][v];
        pre[v] = begin;
        if (!xx[v]) {
            push(v);
        }
    }
    while (!visited[end] && heap_size > 0) {
        u = pop();
        visited[u] = true;
        if (xx[u])
            continue;
        len = adj[u][0];
        for (int j = 1; j <= len; j++) {
            v = adj[u][j];
            if (!visited[v]) {
                if (dist[u] + g[u][v] < dist[v]) {
                    dist[v] = dist[u] + g[u][v];
                    pre[v] = u;
                    if (loc[v] == 0) {
                        push(v);
                    }
                    else {
                        up(v);
                    }
                }
            }
        }
    }

    src = begin;
    dst = end;
    get_path(dst);        
    // for (int i = 0; i < max_v; i++) {
    //     if (i != begin && dist[i] > 0) {
    //         src = begin;
    //         dst = i;
    //         get_path(dst);
    //     }
    // }
}


typedef struct pNode
{
    int node_id;
    int in_degree; 
    int out_degree;
    int* in_node_ids;  /* size = in_degree */
    int* out_node_ids; /* size = out_degree */
} PointNode;

typedef struct tspstruct {
    int s_point;             /* always = 0 */
    int t_point;
    int n;                   /* the number of point = max_v */
    int bestlen;             /* record best cost ever found */
    int best_t_point;        /* t_point of border */
    
    int *iorder;
    int *border;             /* record best order ever found */
    int *iNodePlace;         /* record the place of nodes in iorder */

    int* DD;                 /* all direct_distence; */
    PointNode* node;  /* all nodes infomation */
    bool *isVisted;          /* if visited between s and t */
    bool *isVspec;           /* mark the specific node vetex */

    int* pad;
} TSP;


#define MOD(i,n)    ((i) % (n) >= 0 ? (i) % (n) : (i) % (n) + (n))
#define DirctD(x,y) DirectDist[(x)][(y)]    //get the distence of x,y BY Matrix[x][y]  

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))
#define sqr(x)   ((x)*(x))

int DistD(int* DirectDistance, int max, int x, int y){
    return DirectDistance[x*max + y];
}

//-----------------------init func below--------------------------------------//
    /* 01. init direct distace of any two point. 
     *     input:  max_v, data[700,], w[5000]
     *     out  :  DD is DirectDist[max_v][max_v]
    */
void CF_init_direct_distance(int* DD){
    int i,j;
    memset(DD, 0, sizeof(DD));
    for (i = 0; i < max_v; ++i)
    {
        for (j = 0; j < max_v; j++)
        {
            if (w[data[i][j]] > 0)
            {
                DD[max_v*i+j] = w[data[i][j]]; //DD[i][j]
                // printf("%5d\t", w[data[i][j]]);
            }
            else{
                DD[max_v*i+j] = INF;
                // printf("100000\t");
            }
        }
        // printf("\n");
    }
    // printf("    CF_init_direct_distance: max_v = %d, max_vs = %d:\n", max_v, max_vs);
    return ;
}

    /* 02. init all point nodes infomation 
     *     input : DD is DirectDist , max_v
     *     out   : node[max_v]
    */
void CF_init_point_node(int* DD, PointNode* node)
{
    int i,j;

    for ( i = 0; i < max_v; i++ )
    {
        int in  = 0;
        int out = 0;
        int in_degree_size = 32;
        
        node[i].node_id = i;
        node[i].in_node_ids  = (int *)malloc(in_degree_size * sizeof(int));
        node[i].out_node_ids = (int *)malloc(8 * sizeof(int));

        for( j = 0; j < max_v; j++ ){
            if ( DistD(DD, max_v, j, i) != INF )
            {
                if ( in < in_degree_size )
                {
                    node[i].in_node_ids[in] = j;
                }
                else{  /* realloc memory*/
                    // printf(" in_degree_size=%d\n", in_degree_size);
                    in_degree_size += 16;
                    int* old_ptr = node[i].in_node_ids;
                    int* new_ptr = (int *) realloc(old_ptr, in_degree_size);
                    if (!new_ptr)
                    {
                        DBG(" realloc memory error!");
                        return;
                    }
                    if (new_ptr != old_ptr)
                    {
                        free(node[i].in_node_ids);
                        node[i].in_node_ids = new_ptr;
                    }
                    node[i].in_node_ids[in] = j;
                }
                in++;                
            }
            if ( DistD(DD, max_v, i, j) != INF )
            {
                node[i].out_node_ids[out] = j;
                out++;
            }
        }
        node[i].in_degree = in;
        node[i].out_degree = out;
   }
   return ;
}

int CF_init_st_path(TSP *tsp){
    int i, node;
    int k = 0;

    bool *isVisted;
    bool *isVspec;
    int n, *iorder, *border;
    int *iNodePlace;

    isVisted        = tsp->isVisted;
    isVspec         = tsp->isVspec; 
    iorder          = tsp->iorder;
    border          = tsp->border;
    iNodePlace      = tsp->iNodePlace;

    int t_point = max_vs + 1;
    /* init */
    tsp->s_point = 0;
    tsp->t_point = t_point;
    tsp->best_t_point = t_point;

    /* add s */
    iorder[k] = s;
    isVisted[s] = 1;
    isVspec[s] = 1;
    iNodePlace[s] = k;
    k++;

    /* add V'*/
    for (i = 0; i < max_vs; ++i) {
        node = spec[i];
        isVisted[node] = 1;
        isVspec[node] = 1;
        iorder[k] = node;
        iNodePlace[node] = k;
        k++;
    }

    /* add t */
    isVisted[t] = 1;
    isVspec[t] = 1;
    iorder[k] = t;
    iNodePlace[t] = k;
    k++;
    /* add the point in V-V' */
    for (i = 0; i < max_v; ++i) {
        if (!tsp->isVisted[i]) {
            iorder[k] = i;
            iNodePlace[i] = k;
            k++;
        }
    }

    for (i = 0; i <= t_point; ++i) {
        border[i] = iorder[i];
    }

    return 0;
}

/* slen : the nodes numbers of s->t path 
 * slen = 228, 
 **/
int CF_init_from_feasible_solution(int *order, int slen, TSP *tsp)
{

    int  i, n, node;

    bool *isVisted;
    bool *isVspec;
    int  *iorder;
    int  *iNodePlace;

    isVspec         = tsp->isVspec; 
    isVisted        = tsp->isVisted;
    iorder          = tsp->iorder;
    iNodePlace      = tsp->iNodePlace;
    n               = tsp->n;

    // printf("## CF_init_from_feasible_solution ##\n");
    // printf(" sizeof order = %ld, solution=%ld\n", sizeof(order), sizeof(a_feasible_solution));

    /* init isVspec */    
    isVspec[s] = 1;                       /* add s */
    for (i = 0; i < max_vs; ++i){         /* add V'*/
        node = spec[i];
        isVspec[node] = 1; 
    }  
    isVspec[t] = 1;                       /* add t */

    /* init t_point */
    tsp->s_point = 0;
    tsp->t_point = slen - 1;

    /* init iorder, iNodePlace, isVisted */
    for (i = 0; i < slen; ++i)
    {
        node = order[i];
        iorder[i]           = node;
        iNodePlace[node]    = i;
        isVisted[node]      = 1;
    }
    for (i = slen; i < n; ++i)
    {
        node             = order[i];
        iorder[i]        = node;
        iNodePlace[node] = i;
        isVisted[node]   = 0;
    }

}

int CF_getpathLength (TSP *tsp)
{
    unsigned int i;
    int len = 0;

    int *DD = tsp->DD;
    int *iorder = tsp->iorder;
    int  n      = tsp->n;
    int  t      = tsp->t_point;
    for (i = 0; i < t; i++)
    {
        len += DD[iorder[i]*n + iorder[i+1]];
    }
    return len;
}

//-----------------------init func above--------------------------------------//



//-----------------------debug functions below---------------------------------//

void DBG_CF_print_border(TSP *tsp)
{
    printf("    border is:\n");
    int i;
    
    printf("    (s_point=%d, best_t_point=%d, n=%d)\n", 
            tsp->s_point, tsp-> best_t_point, tsp->n);

    int vNum = 1;
    printf("    border is: ");
    for (i = 0; i <= tsp->best_t_point; ++i)
    {
        printf(" %d", tsp->border[i]);
        if ( tsp->isVspec[ tsp->border[i] ] )
        {
            printf("|(%d)", vNum++);
        }
    }
    printf("\n");

    int len = 0;
    int n = tsp->n;
    int *DD = tsp->DD;
    int *border = tsp->border;
    for (i = 0; i < tsp->best_t_point; i++)
    {
        len += DD[border[i]*n + border[i+1]];
        printf(" (%d, %d)=%d", border[i], border[i+1], DD[border[i]*n + border[i+1]] );
    }
    printf("\n    total border len = %d\n", len);
}

//-----------------------debug functin above--------------------------------//



//---------------------- AS solution lib function below---------------------//
int isConnect(PointNode *a, PointNode *b, TSP *tsp){
    int i;
    int a_out;
    int b_in;
    int *array_a;
    int *array_b;
    bool record[700]; /* record all nodes can reached from a and a->out_nodes
                      * witch is b and b->in_nodes*/

    a_out   = a->out_degree;
    b_in    = b->in_degree;
    array_a = a->out_node_ids;
    array_b = b->in_node_ids;
    memset(&record[0], 0, 700*sizeof(bool));

    // DBG_CF_printTSP(tsp);
    // printf("    a->node_id=%d\n", a->node_id);//-----------------------------------------------//
    /* get record value */
    record[a->node_id] = 1;
    for (i = 0; i < a_out; ++i)
    {
        record[array_a[i]] = 1;
    }

    /* check for connection */
    if (record[b->node_id])
    {
        return b->node_id;
    }
    for (i = 0; i < b_in; ++i)
    {
        int tmp = array_b[i];

        if ( record[tmp] && !(tsp->isVisted[tmp]) ){

            // if (tmp > (tsp->n))
            // {
            //     printf(" tmp = %d, tsp->n=%d, a-id=%d, b-id=%d\n",tmp, tsp->n, a->node_id, b->node_id );
            //     printf("    b:\n");
            //     DBG_CF_print_one_Node(b, tsp);
            //     printf("    a:\n");
            //     DBG_CF_print_one_Node(a, tsp);
            //     printf("\n");
            // }

            return tmp; /* find the first one and then return */
        }
    }
    return -1; /* not found */
}

/*
 * swap (index, tmp_node)<--->(tmp_index, node)
 */
void swap_node(int index, int node, int *iorder, int *iNodePlace, 
                      int *t_point, bool *isVisted)
{
    int tmp_node = iorder[index];
    int tmp_index = iNodePlace[node];

    if (node > 600 || tmp_node > 600)
    {
        printf(" SWAP: ERRRRRRRRRRR NODE=%d, tmp_node=%d > 600\n", node, tmp_node);
    }
    iorder[index] = node;
    iNodePlace[node] = index;

    iorder[tmp_index] = tmp_node;
    iNodePlace[tmp_node] = tmp_index;

    if (index > *t_point)    
        isVisted[node] = 0;    
    else
        isVisted[node] = 1;

    if (tmp_index > *t_point)    
        isVisted[tmp_node] = 0;
    else
        isVisted[tmp_node] = 1;
}

void move_nodes_left(int *iorder, int *iNodePlace, bool *isVisted,
                            int dest, int src, int *t_point)
{
    int i;
    int node = 0;
    int out_num = src - dest;
    int len = *t_point - src + 1;
    int *outs = (int *) malloc(out_num * sizeof(int));
    for (i = dest; i < src; i++)
    {
        outs[i - dest] = iorder[i];
    }
    for ( i = 0; i < len; ++i)
    {
        // if (node > 600 )
        // {
        //     printf(" move left: error NODE=%d > 600\n", node);
        // }
        node = iorder[src+i];
        iorder[dest+i] = node;
        iNodePlace[node] = dest+i;
    }
    for (i = 0; i < out_num; ++i)
    {
        // if (node > 600 )
        // {
        //     printf(" move left: error NODE=%d > 600\n", node);
        // }
        node = outs[i];
        iorder[*t_point] = node;
        iNodePlace[node] = *t_point;
        isVisted[node]   = 0;
        // (*t_point)--;
        *t_point = *t_point - 1;
    }
    return ;
}

/* src : t_point so far */
void move_nodes_right(int *iorder, int *iNodePlace, bool *isVisted,
                             int j_point, int *t_point, int dest)
{
    int i;
    int node = 0;
    int right_num = dest - *t_point;
    int len = *t_point - j_point + 1;  /* +1 is very important*/
    int *tmp = (int *) malloc(right_num * sizeof(int));
    for (i = 0; i < right_num; i++)
    {
        tmp[i] = iorder[dest - i];

    }
    for ( i = 0; i < len; ++i)
    {
        node = iorder[*t_point - i];
        if (node > 600 )
        {
            printf(" move right: ERRRRRRRRRRR NODE=%d > 600\n", node);
        }
        iorder[dest - i] = node;
        iNodePlace[node] = dest - i;
    }
    for (i = 0; i < right_num; ++i)
    {
        if (node > 600 )
        {
            printf(" move right: ERRRRRRRRRRR NODE=%d > 600\n", node);
        }
        node = tmp[i];
        iorder[j_point+i] = node;
        iNodePlace[node]  = j_point+i;
        isVisted[node]    = 1;
    }
    // printf("=====++===t_point=%d, dest=%d\n", *t_point, dest );
    *t_point = dest;
    // printf("==========t_point=%d, dest=%d\n", *t_point, dest );
    return ;
}

/* TODO */
void freeTSP(TSP *tsp){
    return;
}


/*iorder[a]--iorder[a+1] is not connected*/
void CFR_getVsNeighbor(int a, int *ahead, int *atail, int *iorder, bool *isVspec, int t_point){
    int ah, at;
    ah = a;              /* a-head is the index of the first V'-node before a;*/
    at = a + 1;          /* a-tail is the index of the first V'-node behind a;*/
    while( (ah >=  0) && !(isVspec[ iorder[ah] ]) ) ah--;
    while( (at <= t_point) && !(isVspec[ iorder[at] ]) ) at++;
    *ahead = ah;
    *atail = at;
}

int CFR_isValidpath(int s_node, int t_node, bool *isVisted) {
    int i;
    int len = path_each[s_node][t_node][0];

    if (len <= 0 || len > 398) return 0;    
    for (i = 1; i < len; ++i) {
        if (isVisted[ path_each[s_node][t_node][i] ])
            return 0;
    }
    return 1;
} 


int CFR_getSectionLen(int *iorder, int *DD, int hh, int ht) {
    int i;
    int len = 0;
    for (i = hh; i < ht; ++i)
    {
        len += DD[ iorder[i]*max_v + iorder[i+1] ];
    }
    return len;
}


int CFR_try3opt(TSP *tsp, int *pathlen, int *pathchg, double T) {

    int i, j;

    int hh_place, ih_place, jh_place;
    int hh, ht, ih, it, jh, jt;
    int hh_node, ht_node, ih_node, it_node, jh_node, jt_node;
    int *hi_path, *ij_path, *jh_path; 
    int hi_len, ij_len, jh_len;

    int n;
    int *t_point;
    int *iorder;
    int *iNodePlace;
    bool *isVisted;
    bool *isVspec;

    n           = tsp->n;
    t_point     = &tsp->t_point;
    iorder      = tsp->iorder;
    iNodePlace  = tsp->iNodePlace;
    isVisted    = tsp->isVisted;
    isVspec     = tsp->isVspec;

    

    /*get hh_place, ih_place, jh_place*/    
    do {
        hh_place = unifRand(max_vs-1);
    } while(hh_place < 0);

    do {
        ih_place = unifRand(max_vs);
    } while(ih_place < 0 || ih_place <= hh_place );

    do {
        jh_place = unifRand(max_vs+1);
    } while(jh_place < 0 || jh_place <= ih_place);

    /*get hh, ht, ih, it, jh, jt and nodes ID*/
    int count_vs = 0;
    for (i = 0; i <= *t_point; ++i) 
    {
        if ( isVspec[iorder[i]] ) {

            // printf(" i=%d, count_vs=%d\n", i, count_vs);
            if (count_vs == hh_place) {
                hh = i;
            }
            if (count_vs == hh_place+1) {
                ht = i;
            }
            if (count_vs == ih_place) {
                ih = i;
            }
            if (count_vs == ih_place+1) {
                it = i;
            }
            if (count_vs == jh_place) {
                jh = i;
            }
            if (count_vs == jh_place+1) {
                jt = i;
            }

            count_vs = count_vs + 1;

        }
    }

    hh_node = iorder[hh] ;
    ht_node = iorder[ht] ;
    ih_node = iorder[ih] ;
    it_node = iorder[it] ;
    jh_node = iorder[jh] ;
    jt_node = iorder[jt] ;

    /*get *hi_path, *ij_path, *jh_path*/
    memcpy(&xx[0], isVisted, n*sizeof(bool));
    for(i = hh+1; i < ht; i++) {
        xx[ iorder[i] ] = 0;
    }
    for(i = ih+1; i < it; i++) {
        xx[ iorder[i] ] = 0;
    }
    for(i = jh+1; i < jt; i++) {
        xx[ iorder[i] ] = 0;
    }

    if (!CFR_isValidpath(hh_node, it_node, xx)) {
        path_each[hh_node][it_node][0] = 0;
        find_path(hh_node, it_node);
    }
    hi_path = &path_each[hh_node][it_node][0];
    hi_len  = hi_path[0];
    if (hi_len == 0) {
        return -1;
    }
    for (i = 1; i < hi_len; ++i) {
        xx[ hi_path[i] ] = 1;
    }

    if (!CFR_isValidpath(ih_node, jt_node, xx)) {
        path_each[ih_node][jt_node][0] = 0;
        find_path(ih_node, jt_node);
    }
    ij_path = &path_each[ih_node][jt_node][0];
    ij_len  = ij_path[0];
    if (ij_len == 0) {
        return -1;
    }
    for (i = 1; i < ij_len; ++i) {
        xx[ ij_path[i] ] = 1;
    }

    if (!CFR_isValidpath(jh_node, ht_node, xx)) {
        path_each[jh_node][ht_node][0] = 0;
        find_path(jh_node, ht_node);
    }
    jh_path = &path_each[jh_node][ht_node][0];
    jh_len  = jh_path[0];
    if (jh_len == 0) {
        return -1;
    }

    if (jh_len > 398 || hi_len > 398 || ij_len > 398)
    {
        return -1;
        // DBG_CF_print_pathLength(tsp);
        // printf("errrrrror len path\n");
    }

    /*get path len*/
    int c1 = 0;
    int c2 = 0;
    int *DD = tsp->DD;
    c1 = CFR_getSectionLen(iorder, DD, hh, ht) + 
         CFR_getSectionLen(iorder, DD, ih, it) + 
         CFR_getSectionLen(iorder, DD, jh, jt);
    c2 = CFR_getSectionLen(hi_path, DD, 1, hi_len) + DD[hh_node*n + hi_path[1]] +
         CFR_getSectionLen(jh_path, DD, 1, jh_len) + DD[jh_node*n + jh_path[1]] +
         CFR_getSectionLen(ij_path, DD, 1, ij_len) + DD[ih_node*n + ij_path[1]]; 

    // DBG_CF_print_pathLength(tsp);
    int energyChange = c2 - c1;

    if ( energyChange < 0 || PARA_02 < exp(-energyChange/T) ) 
    {
        int in, out;
        int move_len;
        in  = hi_len  + jh_len  + ij_len  - 3;
        out = ht - hh + it - ih + jt - jh - 3; 

        int path_left[700];
        int path_right[700];
        path_left[0]  = ih - ht + 1;  /*path lenght*/
        path_right[0] = jh - it + 1;  
        for (i = 1; i <= path_left[0]; ++i) 
        {
            path_left[i] = iorder[ ht+i-1 ];
        }
        for (i = 1; i <= path_right[0]; ++i)
        {
            path_right[i] = iorder[ it+i-1 ];
        }

        int tmp_len  =  CF_getpathLength(tsp);
        int flag;
        if (out < in) {
            move_len = in - out;
            move_nodes_right(iorder, iNodePlace, isVisted, jt, t_point, *t_point+move_len);
        }
        else if (out > in) {
            move_len = out - in;
            move_nodes_left(iorder, iNodePlace, isVisted, jt - move_len, jt, t_point);
        }

        flag = hh;
        for (i = 1; i < hi_len; i++) {
            swap_node(flag+i, hi_path[i], iorder, iNodePlace, t_point, isVisted);
        }

        flag = flag + hi_len - 1;
        for(i = 1; i <= path_right[0]; i++) {
            swap_node(flag+i, path_right[i], iorder, iNodePlace, t_point, isVisted);
        }

        flag = flag + path_right[0];
        if (flag < 0) {
            return 0;
            printf("errerror\n");
        }
        for(i = 1; i < jh_len; i++) {
            swap_node(flag+i, jh_path[i], iorder, iNodePlace, t_point, isVisted);
        }

        flag = flag + jh_len -1;
        for(i = 1; i <= path_left[0]; i++) { /*YOU HAVE TO BE VERY CAREFUL in writing, [0] not [i]*/
            swap_node(flag+i, path_left[i], iorder, iNodePlace, t_point, isVisted);
        }

        flag = flag + path_left[0];
        for(i = 1; i < ij_len; i++) {
            swap_node(flag +i, ij_path[i], iorder, iNodePlace, t_point, isVisted);
        }

        // DBG_CF_print_pathLength(tsp);
        *pathchg++;
        *pathlen = CF_getpathLength(tsp);
        // *pathlen += energyChange;
        // printf(" pathlen = %d\n", *pathlen);
        // printf(" CF_getpathLength = %d\n", CF_getpathLength(tsp));
        // if (*pathlen == 1469 && tsp->t_point == 143)
        // {
        //     DBG_CF_print_pathLength(tsp);
        //     printf("====\n");
        // }
    }
    return 0;
}
//---------------------- AS solution lib function above---------------------//



void annealing(TSP *tsp)
{

    int    i, j, pathchg, nochgtimes;
    int    *pathlen;
    int    energyChange;
    double T;

    pathlen = (int *) malloc(sizeof(int));
    *pathlen = CF_getpathLength(tsp);
    tsp->bestlen = *pathlen;
    nochgtimes = 0;

    for (T = AS_T_INIT; T > AS_FINAL_T; T *= AS_COOLING)  /* annealing schedule */
    {
        pathchg = 0;
        for (j = 0; j < AS_TRIES_PER_T; j++)
        {

            CFR_try3opt(tsp, pathlen, &pathchg, T);              

            if ( *pathlen < (tsp->bestlen) ) {
                tsp->bestlen = *pathlen;
                tsp->best_t_point = tsp->t_point;
                for (i=0; i <= tsp->t_point; i++) 
                    tsp->border[i] = tsp->iorder[i];
                if (*pathlen < INF)
                {
                    return;
                }
            }

            if (pathchg > AS_IMPROVED_PATH_PER_T) break; /* finish early */

        }   
        if (pathchg == 0) nochgtimes++;
        if (nochgtimes == AS_no_change_time) break;
    }
}

/* INIT TSP DATA from data[700][700]*/
int find_Craft_solution(int n, int start, int end, int *res)
{
    int   i;
    TSP   tsp; 

    /*init random function*/
    // long  seed = -314159L;
    initRand (AS_seed);
    srand(AS_seed);

    /*initialize tsp struct*/
    tsp.n = n;
    tsp.iorder = NULL;
    tsp.border = NULL;
    tsp.DD     = NULL;
    // tsp.node   = NULL;
    
    tsp.iNodePlace = NULL;
    tsp.isVisted   = NULL;
    tsp.isVspec    = NULL;

    if (!(tsp.iorder = (int*) malloc (tsp.n * sizeof(int)))         ||
        !(tsp.border = (int*) malloc (tsp.n * sizeof(int)))         ||
        !(tsp.iNodePlace = (int*) malloc (n * sizeof(int)))         ||
        !(tsp.DD = (int*) malloc (n * n *sizeof(int)))              /*||
        !(tsp.node = (PointNode *) malloc (n * sizeof(PointNode)))   */  ) {
            return -1;
            printf( "Memory allocation failed!");
        }

    tsp.isVisted = (bool *) malloc (n * sizeof(bool));
    tsp.isVspec = (bool *) malloc (n * sizeof(bool)); 
    memset(tsp.isVisted, 0, sizeof(n * sizeof(bool)));
    memset(tsp.isVspec, 0, sizeof(n * sizeof(bool)));

    CF_init_direct_distance(tsp.DD);
    // CF_init_point_node(tsp.DD, tsp.node);
    CF_init_st_path(&tsp);

    annealing(&tsp);

    *res = tsp.bestlen;
    path[0] = tsp.best_t_point+1;
    for (i = 1; i <= path[0]; i++) 
        path[i] = tsp.border[i-1];

    DBG("Final Path Length(*total_len): %d\n", *res);
    // DBG_CF_print_border(&tsp);
    return 0;
}

int find_Better_Craft_solution(int n, int start, int end, int *res)
{
    int   i;
    TSP   tsp; 

    /*init random function*/
    // long  seed = -314159L;
    initRand (AS_seed);
    srand(AS_seed);

    /*initialize tsp struct*/
    tsp.n = n;
    tsp.iorder = NULL;
    tsp.border = NULL;
    tsp.DD     = NULL;
    // tsp.node   = NULL;
    
    tsp.iNodePlace = NULL;
    tsp.isVisted   = NULL;
    tsp.isVspec    = NULL;

    if (!(tsp.iorder = (int*) malloc (tsp.n * sizeof(int)))         ||
        !(tsp.border = (int*) malloc (tsp.n * sizeof(int)))         ||
        !(tsp.iNodePlace = (int*) malloc (n * sizeof(int)))         ||
        !(tsp.DD = (int*) malloc (n * n *sizeof(int)))              /*||
        !(tsp.node = (PointNode *) malloc (n * sizeof(PointNode)))   */  ) {
            return -1;
            printf( "Memory allocation failed!");
        }

    tsp.isVisted = (bool *) malloc (n * sizeof(bool));
    tsp.isVspec = (bool *) malloc (n * sizeof(bool)); 
    memset(tsp.isVisted, 0, sizeof(n * sizeof(bool)));
    memset(tsp.isVspec, 0, sizeof(n * sizeof(bool)));

    CF_init_direct_distance(tsp.DD);
    // CF_init_point_node(tsp.DD, tsp.node);
    // CF_init_st_path(&tsp);
    CF_init_from_feasible_solution(a_feasible_solution, slength, &tsp);

    annealing(&tsp);

    *res = tsp.bestlen;
    path[0] = tsp.best_t_point+1;
    for (i = 1; i <= path[0]; i++) 
        path[i] = tsp.border[i-1];

    DBG("Final Path Length(*total_len): %d\n", *res);
    // DBG_CF_print_border(&tsp);
    return 0;
}


//---------------------------IO code below ----------------------------------------//
void read_topo_for9(char * topo) {
    memset(data, 0xff, sizeof(data));
    memset(w, 0, sizeof(w));
    memset(adj, 0, sizeof(adj));
    memset(prt, 0, sizeof(prt));
    memset(g, 0x7f, sizeof(g));
    memset(a, 0, sizeof(a));
    FILE *input = fopen(topo, "r");
    int id;
    int x;
    int y;
    int weight;
    int i, j;
    max_v = 0;
    max_e = 0;
    // read data[][]
    // g[][] = data[][]
    while (!feof(input)) {
        if (fscanf(input, "%d,%d,%d,%d", &id, &x, &y, &weight) > 0) {
            if (id + 1 > max_e) max_e = id + 1;
            if (x + 1 > max_v) max_v = x + 1;
            if (y + 1 > max_v) max_v = y + 1;
            if (x != t && y != s) {
                if (weight < g[x][y]) {
                    data[x][y] = id;
                    g[x][y] = weight;
                    a[x][y] = 1;
                }
            }
            w[id] = weight;
        }
    }
    // reduce topo
    int len;
    for (int i = 0; i < max_vs; i++) {
        int u = spec[i];
        int v;
        // if out edge # = 1
        len = 0;
        for (int j = 0; j < max_v; j++) {
            if (data[u][j] >= 0) {
                len++;
                v = j;
            }
        }
        if (len == 1) {
            for (int j = 0; j < max_v; j++) {
                if (data[j][v] >= 0 && j != u) {
                    data[j][v] = -1;
                    g[j][v] = 0x7fffffff;
                    a[j][v] = 0;
                }
            }
        }
        // if in edge # = 1
        len = 0;
        for (int j = 0; j < max_v; j++) {
            if (data[j][u] >= 0) {
                len++;
                v = j;
            }
        }
        if (len == 1) {
            for (int j = 0; j < max_v; j++) {
                if (data[v][j] >= 0 && j != u) {
                    data[v][j] = -1;
                    g[v][j] = 0x7fffffff;
                    a[v][j] = 0;
                }
            }
        }
    }
    // build adj & prt list
    for (int i = 0; i < max_v; i++) {
        for (int j = 0; j < max_v; j++) {
            if (i != j && data[i][j] >= 0) {
                adj[i][0]++;
                adj[i][adj[i][0]] = j;
                prt[j][0]++;
                prt[j][prt[j][0]] = i;
            }
        }
    }
    // heuristic : visit shorter neighbor first
    // sort adj list by distance
    for (i = 0; i < max_v; i++) {
        len = adj[i][0];
        for (j = 1; j < len; j++) {
            int k = j + 1;
            while (k > 1 && g[i][adj[i][k]] < g[i][adj[i][k - 1]]) {
                adj[i][k] = adj[i][k] ^ adj[i][k - 1];
                adj[i][k - 1] = adj[i][k] ^ adj[i][k - 1];
                adj[i][k] = adj[i][k] ^ adj[i][k - 1];
                k--;
            }
        }
    }
    fclose(input);
}

void read_topo(char * topo) {
    memset(data, 0xff, sizeof(data));
    memset(w, 0, sizeof(w));
    memset(adj, 0, sizeof(adj));
    memset(prt, 0, sizeof(prt));
    memset(g, 0x7f, sizeof(g));
    memset(a, 0, sizeof(a));
    FILE *input = fopen(topo, "r");
    int id;
    int x;
    int y;
    int weight;
    int i, j;
    max_v = 0;
    max_e = 0;
    // read data[][]
    // g[][] = data[][]
    while (!feof(input)) {
        if (fscanf(input, "%d,%d,%d,%d", &id, &x, &y, &weight) > 0) {
            if (id + 1 > max_e) max_e = id + 1;
            if (x + 1 > max_v) max_v = x + 1;
            if (y + 1 > max_v) max_v = y + 1;
            if (x != t && y != s) {
                if (data[x][y] < 0) {
                    adj[x][0]++;
                    adj[x][adj[x][0]] = y;
                    prt[y][0]++;
                    prt[y][prt[y][0]] = x;
                    a[x][y] = 1;
                }
                if (weight < g[x][y]) {
                    data[x][y] = id;
                    g[x][y] = weight;
                }
            }
            w[id] = weight;
        }
    }
    // heuristic : visit shorter neighbor first
    // sort adj list by distance
    int len;
    for (i = 0; i < max_v; i++) {
        len = adj[i][0];
        for (j = 1; j < len; j++) {
            int k = j + 1;
            while (k > 1 && g[i][adj[i][k]] < g[i][adj[i][k - 1]]) {
                adj[i][k] = adj[i][k] ^ adj[i][k - 1];
                adj[i][k - 1] = adj[i][k] ^ adj[i][k - 1];
                adj[i][k] = adj[i][k] ^ adj[i][k - 1];
                k--;
            }
        }
    }
    fclose(input);
}

void read_demand(char* demand) {
    memset(spec_2, 0, sizeof(spec_2));
    memset(spec, 0, sizeof(spec));
    max_vs = 0;
    int x;
    FILE *input = fopen(demand, "r");
    fscanf(input, "%d,%d,", &s, &t);
    while (!feof(input)) {
        if (fscanf(input, "%d|", &x) > 0) {
            spec_2[x] = true;
            spec[max_vs] = x;
            max_vs++;
        }
    }
    fclose(input);
}

//---------------------- DFS solution below ---------------------//

void find_path_to_t() {
    // Dijkstra
    memset(dist, 0x7f, sizeof(dist));
    memset(pre, 0xff, sizeof(pre));
    memset(heap, 0, sizeof(heap));
    memset(visited, 0, sizeof(visited));
    memset(loc, 0, sizeof(loc));
    visited[t] = true;
    dist[t] = 0;
    int u, v;
    heap_size = 0;
    int len = prt[t][0];
    for (int i = 1; i <= len; i++) {
        v = prt[t][i];
        dist[v] = g[v][t];
        push(v);
    }
    while (heap_size > 0) {
        u = pop();
        visited[u] = true;
        len = prt[u][0];
        for (int j = 1; j <= len; j++) {
            v = prt[u][j];
            if (!visited[v]) {
                if (dist[u] + g[v][u] < dist[v]) {
                    dist[v] = dist[u] + g[v][u];
                    if (loc[v] == 0) {
                        push(v);
                    }
                    else {
                        up(v);
                    }
                }
            }
        }
    }
}

void init_g() {
    find_path_to_t();
    memset(dist_to_dst, 0x7f, sizeof(dist_to_dst));
    for (int i = 0; i < max_v; i++) {
        dist_to_dst[i] = dist[i];
    }
    // heuristic : visit neighbor of shorter distance to t first
    // sort adj list by distance
    int len;
    for (int i = 0; i < max_v; i++) {
        len = adj[i][0];
        for (int j = 1; j < len; j++) {
            int k = j + 1;
            while (k > 1 && dist_to_dst[adj[i][k]] < dist_to_dst[adj[i][k - 1]]) {
                adj[i][k] = adj[i][k] ^ adj[i][k - 1];
                adj[i][k - 1] = adj[i][k] ^ adj[i][k - 1];
                adj[i][k] = adj[i][k] ^ adj[i][k - 1];
                k--;
            }
        }
    }
}

void get_result(int v) {
    if (pre[v] >= 0) {
        get_result(pre[v]);
    }
    path[0]++;
    path[path[0]] = v;
}

void push_stack(int x) {
    stack[top] = x;
    top++;
}

int pop_stack() {
    top--;
    return stack[top];
}

void dfs(int u, int cost) {
    if (u == t && vs_cnt == max_vs) {
        if (cost < best) {
            memset(path, 0xff, sizeof(path));
            path[0] = 0;
            get_result(t);
            best = cost;
        }
        return;
    }
    int len;
    int v;
    bool check;
    int len2;
    int x;
    int cnt_in;
    int cnt_out;
    len = adj[u][0];
    for (int i = 1; i <= len; i++) {
        v = adj[u][i];
        check = true;
        if (a[u][v]) {
            if (spec_2[v])
                vs_cnt++;
            f[u][v]++;
            row[u]++;
            col[v]++;
            push_stack(pre[v]);
            pre[v] = u;
            len2 = prt[v][0];
            cnt_in = 0;
            for (int j = 1; j <= len; j++) {
                x = adj[u][j];
                if (a[u][x]) {
                    a[u][x]--;
                    r[u]--;
                    c[x]--;
                    if (spec_2[x] && c[x] + col[x] == 0)
                        check = false;
                    cnt_in++;
                    push_stack(x);
                }
            }
            cnt_out = 0;
            for (int j = 1; j <= len2; j++) {
                x = prt[v][j];
                if (a[x][v]) {
                    a[x][v]--;
                    r[x]--;
                    c[v]--;
                    if (spec_2[x] && r[x] + row[x] == 0) {
                        check = false;
                    }
                    cnt_out++;
                    push_stack(x);
                }
            }
            if (check && cost + g[u][v] + dist_to_dst[v] < best) {
                // go next iteration
                dfs(v, cost + g[u][v]);
            }
            // track back
            for (int j = 0; j < cnt_out; j++) {
                x = pop_stack();
                c[v]++;
                r[x]++;
                a[x][v] = 1;
            }
            for (int j = 0; j < cnt_in; j++) {
                x = pop_stack();
                c[x]++;
                r[u]++;
                a[u][x] = 1;
            }
            pre[v] = pop_stack();
            col[v]--;
            row[u]--;
            f[u][v]--;
            if (spec_2[v])
                vs_cnt--;
        }
    }
}

int solve_dfs() {
    init_g();
    memset(r, 0, sizeof(r));
    memset(c, 0, sizeof(c));
    for (int i = 0; i < max_v; i++) {
        a[i][s] = 0;
        a[t][i] = 0;
        r[i] = adj[i][0];
        c[i] = prt[i][0];
    }
    memset(row, 0, sizeof(row));
    memset(col, 0, sizeof(col));
    best = max_cost;
    vs_cnt = 0;
    memset(stack, 0, sizeof(stack));
    memset(pre, 0xff, sizeof(pre));
    top = 0;
    dfs(s, 0);
    return best;
}

//---------------------- DFS solution above ---------------------//

//---------------------- ACS solution below ---------------------//

void solve_acs(int *res) {
    // to do

    *res = max_cost;
}

//---------------------- ACS solution above ---------------------//

void print_result(char* filename, int res) {
    //printf("case %d-%d : ", max_v, max_vs);
    FILE *output = fopen(filename, "w");
    if (res < 12001) {
        //fprintf(output, "%d\n", res);
        //printf("result path length = %d\n", best);
        int len = path[0];
        for (int i = 2; i < len; i++) {
            fprintf(output, "%d|", data[path[i - 1]][path[i]]);
        }
        fprintf(output, "%d", data[path[len - 1]][path[len]]);
    }
    else {
        //printf("NA\n");
        fprintf(output, "NA\n");
    }
    fclose(output);
}

//---------------------------main function ----------------------------------------//
int main(int argc, char *argv[])
{
    char *topo = argv[1];
    char *demand = argv[2];
    char *output = argv[3];

    int res;
    read_demand(demand);
    read_topo(topo);
    
    /*BEST for 1-5*/
    if (max_v < 100) {
        AS_no_change_time = 800;
        AS_TRIES_PER_T = 500;
        AS_small_case = true;
        res = solve_dfs();
    }
    // for 14
    else if (max_v > 500 && max_vs > 25) {
        result = argv[3];
        ACS_read_topo(topo);
        ACS_read_demand(demand);
        Change_spec();
        Change_adj();
        Bst_Agent.path_len = UN_REACH;
        Agent_number = 34; //35
        Opt_time    =  500;  //500
        Phe_keep    =  0.98;    //0.90
        Phe_fac     =  1.0;     //1.0 信息素因子，信息素的重要程度
        Inspir_fac  =  2.5;     //2.5 启发因子
        ACS_FOR_SearchWay();

        AS_no_change_time = 40; //32
        AS_TRIES_PER_T = 3000; //3000
        AS_seed = -314159L;
        AS_IMPROVED_PATH_PER_T = 50; //50
        find_Better_Craft_solution(max_v, s, t, &res);
        print_result(output, res);
        return 0;
    }
    // for 15
    else if (max_v > 500 ) {
        result = argv[3];
        ACS_read_topo(topo);
        ACS_read_demand(demand);
        Change_spec();
        Change_adj();
        Bst_Agent.path_len = UN_REACH;
        Opt_time    =  500;
        Phe_keep    =  0.98;    //0.89 /0.98
        Phe_fac     =  1.0;     //1.0 /1.0 
        Inspir_fac  =  2.5;     //0.6 /0.6 
        ACS_FOR_SearchWay();

        AS_no_change_time = 25; //25
        AS_TRIES_PER_T = 3000; //3000
        AS_seed = -314159L;
        AS_IMPROVED_PATH_PER_T = 50; //50
        find_Better_Craft_solution(max_v, s, t, &res);
        // find_Craft_solution(max_v, s, t, &res);
        print_result(output, res);
        return 0;
    }
    else {
        // return 0;
        /*BEST for 6*/
        if (max_v == 100) {
            AS_no_change_time = 70; //so far best 75
        }
        /*BEST for 7*/
        else if (max_v <= 150) {
            AS_no_change_time = 65; //so far best 68
            for (int i = 0; i < max_v; i++) {
                for (int j = 0; j < max_v; j++) {
                    g[i][j] = data[i][j];
                }
            }
        }
        /*BEST for 8*/
        else if (max_v == 200) {
            AS_no_change_time = 45; //45
            AS_TRIES_PER_T = 500;
        }
        /*BEST for 9*/
        else if (max_v == 250) {
            read_topo_for9(topo);
            AS_no_change_time = 400; //400
            AS_TRIES_PER_T = 410; //410
            // AS_seed = -314159L;
        }
        /*BEST for 10*/
        else if (max_v == 300) {
            AS_no_change_time = 400; //400
            AS_TRIES_PER_T = 410; //410
        }
        /*BEST for 11*/
        else if (max_v == 500 && max_vs > 30) {
            AS_no_change_time = 400; //400
            AS_TRIES_PER_T = 730; //730
            AS_seed = -314159L;
            AS_IMPROVED_PATH_PER_T = 50; //50
        }
        /*for 12 */
        else if (max_v == 500 && max_vs < 23) {
            AS_no_change_time = 600; //so far best 600
            AS_TRIES_PER_T = 700; // 700
            AS_seed = -314159L;
            AS_IMPROVED_PATH_PER_T = 50; //50
        }
        /*for 13*/
        else if (max_v == 500) {
            AS_no_change_time = 600; //so far best 600
            AS_TRIES_PER_T = 700; // 700
            AS_seed = -314159L;
            AS_IMPROVED_PATH_PER_T = 50; //50
        }
        find_Craft_solution(max_v, s, t, &res);
    }

    print_result(output, res);
    return 0;
}  

