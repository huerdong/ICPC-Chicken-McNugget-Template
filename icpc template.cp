/*
 ID: areke
 PROG: fence
 LANG: C++
 */

#include <iostream>
#include <string>
#include <math.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <vector>
#include <map>
#include <set>
#include <stack>
#include <limits>
#include <queue>
#include <list>
#include <cstring>
using namespace std;

long long c[100000];

long long find(vector<long long> &v, long long x) {
    return (v[x] == x) ? x : v[x] = find(v, v[x]);
}


void merge(vector<long long> &v, long long x, long long y) {
    
    long long fx = find(v, x);
    long long fy = find(v, y);
    if (fx == fy) c[fx]++;
    else {
        v[fy] = fx;
        c[fx] += c[fy];
    }
}


long long gcd(long long a, long long b) {
    return b == 0 ? a : gcd(b, a % b);
}


int t[100000];
void kmpTable(string w) {
    int pos = 2;
    int cnd = 0;
    t[0] = -1;
    t[1] = 0;
    while (pos < w.length()) {
        if (w[pos - 1] == w[cnd]) {
            t[pos] = cnd + 1;
            cnd++;
            pos++;
        } else if (cnd > 0) {
            cnd = t[cnd];
            t[pos] = 0;
        } else {
            t[pos] = 0;
            pos++;
        }
    }
}

int kmpSubstr(string s, string w) {
    int m = 0;
    int i = 0;
    while (m + i < s.length()) {
        if (w[i] == s[m + i]) {
            if (i == w.length() - 1) return m;
            i++;
        } else {
            if (t[i] > -1) {
                m = m + i - t[i];
                i = t[i];
            }
            else {
                m++;
                i = 0;
            }
        }
    }
    return -1;
}

vector<vector<long long> > v;

int circuit[2048] = {0};
int s = 0;

void deleteNode(int a, int b) {
    for (int i = 0; i < v[a].size(); i++) {
        if (v[a][i] == b) {
            v[a].erase(v[a].begin() + i);
            break;
        }
    }
}

// Find Eulerian Tour

void findCircuit(int n) {
    if (v[n].size() == 0) {
        circuit[s] = n;
        s++;
    }
    else {
        while (v[n].size() > 0) {
            int u = v[n][0];
            deleteNode(n, u);
            deleteNode(u, n);
            findCircuit(u);
        }
        circuit[s] = n;
        s++;
    }
}

bool comp(const long long &a, const long long & b) {
    return a < b;
}
int main() {
    ifstream fin("fence.in");
    ofstream fout("fence.out");
    long long n;
    fin >> n;
    
    v.resize(500);
    for (int i = 0; i < n; i++) {
        long long x, y;
        fin >> x >> y;
        x--;
        y--;
        v[x].push_back(y);
        v[y].push_back(x);
    }
    
    for (int i = 0; i < 500; i++) {
        sort(v[i].begin(), v[i].end(), comp);
    }
    int ind = 0;
    for (int i = 0; i < 500; i++) {
        if (v[i].size() % 2 == 1) {
            ind = i;
            break;
        }
    }
    
    findCircuit(ind);
    /*
    stack<int> st;
    
    st.push(ind);
    
    while (!st.empty()) {
        int node = st.top();
        st.pop();
        if (v[node].size() == 0) {
            circuit[s] = node;
            s++;
        }
        else {
            st.push(node);
            int u = v[node][0];
            deleteNode(node, u);
            deleteNode(u, node);
            st.push(u);
        }
    }*/
    
    cout << s << endl;
    for (int i = s-1; i >= 0; i--) {
       fout << circuit[i] + 1 << endl;
    }
}

// maxflow

long long source = 0;
long long sink = m - 1;
long long c[200][200];
long long totalflow = 0;

if (source == sink) {
    totalflow = inf;
}

while (true) {
    long long prevnode[200];
    long long flow[200];
    bool visited[200];
    for (int i = 0; i < m; i++) {
        prevnode[i] = nil;
        flow[i] = 0;
        visited[i] = false;
    }
    
    flow[source] = inf;
    long long maxflow, maxloc;
    while (true) {
        maxflow = 0;
        maxloc = nil;
        for (int i = 0; i < m; i++) {
            if (flow[i] > maxflow && !visited[i]) {
                maxflow = flow[i];
                maxloc = i;
            }
        }
        if (maxloc == nil || maxloc == sink) break;
        visited[maxloc] = true;
        for (int i = 0; i < m; i++) {
            if (i != maxloc && flow[i] < min(maxflow, c[maxloc][i])) {
                prevnode[i] = maxloc;
                flow[i] = min(maxflow, c[maxloc][i]);
            }
        }
    }
    if (maxloc == nil) break;
    long long pc = flow[sink];
    totalflow += pc;
    long long curnode = sink;
    while (curnode != source) {
        long long nextnode = prevnode[curnode];
        c[nextnode][curnode] -= pc;
        c[curnode][nextnode] += pc;
        curnode = nextnode;
    }
}
fout << totalflow << endl;


/**
 *   ///////////////////////
 *   // MIN COST MAX FLOW //
 *   ///////////////////////
 *
 *   Authors: Frank Chu, Igor Naverniouk
 **/

/*********************
 * Min cost max flow * (Edmonds-Karp relabelling + Dijkstra)
 *********************
 * Takes a directed graph where each edge has a capacity ('cap') and a
 * cost per unit of flow ('cost') and returns a maximum flow network
 * of minimal cost ('fcost') from s to t. USE THIS CODE FOR (MODERATELY)
 * DENSE GRAPHS; FOR VERY SPARSE GRAPHS, USE mcmf4.cpp.
 *
 * PARAMETERS:
 *      - cap (global): adjacency matrix where cap[u][v] is the capacity
 *          of the edge u->v. cap[u][v] is 0 for non-existent edges.
 *      - cost (global): a matrix where cost[u][v] is the cost per unit
 *          of flow along the edge u->v. If cap[u][v] == 0, cost[u][v] is
 *          ignored. ALL COSTS MUST BE NON-NEGATIVE!
 *      - n: the number of vertices ([0, n-1] are considered as vertices).
 *      - s: source vertex.
 *      - t: sink.
 * RETURNS:
 *      - the flow
 *      - the total cost through 'fcost'
 *      - fnet contains the flow network. Careful: both fnet[u][v] and
 *          fnet[v][u] could be positive. Take the difference.
 * COMPLEXITY:
 *      - Worst case: O(n^2*flow  <?  n^3*fcost)
 * FIELD TESTING:
 *      - Valladolid 10594: Data Flow
 * REFERENCE:
 *      Edmonds, J., Karp, R.  "Theoretical Improvements in Algorithmic
 *          Efficieincy for Network Flow Problems".
 *      This is a slight improvement of Frank Chu's implementation.
 **/

#include <string.h>
#include <limits.h>
using namespace std;

// the maximum number of vertices + 1
#define NN 1024

// adjacency matrix (fill this up)
int cap[NN][NN];

// cost per unit of flow matrix (fill this up)
int cost[NN][NN];

// flow network and adjacency list
int fnet[NN][NN], adj[NN][NN], deg[NN];

// Dijkstra's successor and depth
int par[NN], d[NN];        // par[source] = source;

// Labelling function
int pi[NN];

#define CLR(a, x) memset( a, x, sizeof( a ) )
#define Inf (INT_MAX/2)

// Dijkstra's using non-negative edge weights (cost + potential)
#define Pot(u,v) (d[u] + pi[u] - pi[v])
bool dijkstra( int n, int s, int t )
{
    for( int i = 0; i < n; i++ ) d[i] = Inf, par[i] = -1;
    d[s] = 0;
    par[s] = -n - 1;
    
    while( 1 )
    {
        // find u with smallest d[u]
        int u = -1, bestD = Inf;
        for( int i = 0; i < n; i++ ) if( par[i] < 0 && d[i] < bestD )
            bestD = d[u = i];
        if( bestD == Inf ) break;
        
        // relax edge (u,i) or (i,u) for all i;
        par[u] = -par[u] - 1;
        for( int i = 0; i < deg[u]; i++ )
        {
            // try undoing edge v->u
            int v = adj[u][i];
            if( par[v] >= 0 ) continue;
            if( fnet[v][u] && d[v] > Pot(u,v) - cost[v][u] )
                d[v] = Pot( u, v ) - cost[v][u], par[v] = -u-1;
            
            // try edge u->v
            if( fnet[u][v] < cap[u][v] && d[v] > Pot(u,v) + cost[u][v] )
                d[v] = Pot(u,v) + cost[u][v], par[v] = -u - 1;
        }
    }
    
    for( int i = 0; i < n; i++ ) if( pi[i] < Inf ) pi[i] += d[i];
    
    return par[t] >= 0;
}
#undef Pot

int mcmf3( int n, int s, int t, int &fcost )
{
    // build the adjacency list
    CLR( deg, 0 );
    for( int i = 0; i < n; i++ )
        for( int j = 0; j < n; j++ )
            if( cap[i][j] || cap[j][i] ) adj[i][deg[i]++] = j;
    
    CLR( fnet, 0 );
    CLR( pi, 0 );
    int flow = fcost = 0;
    
    // repeatedly, find a cheapest path from s to t
    while( dijkstra( n, s, t ) )
    {
        // get the bottleneck capacity
        int bot = INT_MAX;
        for( int v = t, u = par[v]; v != s; u = par[v = u] )
            bot <?= fnet[v][u] ? fnet[v][u] : ( cap[u][v] - fnet[u][v] );
        
        // update the flow network
        for( int v = t, u = par[v]; v != s; u = par[v = u] )
            if( fnet[v][u] ) { fnet[v][u] -= bot; fcost -= bot * cost[v][u]; }
            else { fnet[u][v] += bot; fcost += bot * cost[u][v]; }
        
        flow += bot;
    }
    
    return flow;
}

//----------------- EXAMPLE USAGE -----------------
#include <iostream>
#include <stdio.h>
using namespace std;

int main()
{
    int numV;
    //  while ( cin >> numV && numV ) {
    cin >> numV;
    memset( cap, 0, sizeof( cap ) );
    
    int m, a, b, c, cp;
    int s, t;
    cin >> m;
    cin >> s >> t;
    
    // fill up cap with existing capacities.
    // if the edge u->v has capacity 6, set cap[u][v] = 6.
    // for each cap[u][v] > 0, set cost[u][v] to  the
    // cost per unit of flow along the edge i->v
    for (int i=0; i<m; i++) {
        cin >> a >> b >> cp >> c;
        cost[a][b] = c; // cost[b][a] = c;
        cap[a][b] = cp; // cap[b][a] = cp;
    }
    
    int fcost;
    int flow = mcmf3( numV, s, t, fcost );
    cout << "flow: " << flow << endl;
    cout << "cost: " << fcost << endl;
    
    return 0;
}

// short O(n^2) longest nondecreasing subsequence
#include <stdio.h>
#define SIZE 200000
#define MAX(x,y) ((x)>(y)?(x):(y))

int     best[SIZE];        // best[] holds values of the optimal sub-sequence

int main (void) {
    FILE *in  = fopen ("input.txt", "r");
    FILE *out = fopen ("output.txt", "w");
    int     i, n, k, x, sol = -1;
    
    fscanf (in, "%d", &n);	// N = how many integers to read in
    for (i = 0; i < n; i++) {
        best[i] = -1;
        fscanf (in, "%d", &x);
        for (k = 0; best[k] > x; k++)
            ;
        best[k] = x;
        sol = MAX (sol, k + 1);
    }
    printf ("best is %d\n", sol);
    return 0;
}

// O(nlogn) lnds
#include <stdio.h>
#define SIZE 200000
#define MAX(x,y) ((x)>(y)?(x):(y))

int     best[SIZE];        // best[] holds values of the optimal sub-sequence

int
main (void) {
    FILE *in  = fopen ("input.txt", "r");
    int i, n, k, x, sol;
    int low, high;
    
    fscanf (in, "%d", &n);	// N = how many integers to read in
    // read in the first integer
    fscanf (in, "%d", &best[0]);
    sol = 1;
    for (i = 1; i < n; i++) {
        best[i] = -1;
        fscanf (in, "%d", &x);
        
        if(x >= best[0]) {
            k = 0;
            best[0] = x;
        }
        else {
            // use binary search instead
            low = 0;
            high = sol-1;
            for(;;) {
                k = (int) (low + high) / 2;
                // go lower in the array
                if(x > best[k] && x > best[k-1]) {
                    high = k - 1;
                    continue;
                }
                // go higher in the array
                if(x < best[k] && x < best[k+1]) {
                    low = k + 1;
                    continue;
                }
                // check if right spot
                if(x > best[k] && x < best[k-1])
                    best[k] = x;
                if(x < best[k] && x > best[k+1])
                    best[++k] = x;
                break;
            }
        }
        sol = MAX (sol, k + 1);
    }
    printf ("best is %d\n", sol);
    fclose(in);
    return 0;
}


//lcsubsequence
int lcs( char *X, char *Y, int m, int n )
{
    int L[m+1][n+1];
    int i, j;
    
    /* Following steps build L[m+1][n+1] in bottom up fashion. Note
     that L[i][j] contains length of LCS of X[0..i-1] and Y[0..j-1] */
    //long long result = 0
    for (i=0; i<=m; i++)
    {
        for (j=0; j<=n; j++)
        {
            if (i == 0 || j == 0)
                L[i][j] = 0;
            
            else if (X[i-1] == Y[j-1])
                L[i][j] = L[i-1][j-1] + 1;
            //for substring result = max(result, LCSuff[i][j]);
            
            else
                L[i][j] = max(L[i-1][j], L[i][j-1]);
                // for substring l[i][j] = 0;
        }
    }
    
    /* L[m][n] contains length of LCS for X[0..n-1] and Y[0..m-1] */
    return L[m][n];
}

//totient
void phi(int n) {
    //Sieve of eratosthenes
    vector<int> sieve;
    sieve.push_back(2);
    for (int i = 3; i < n; ++i) {
        bool flag = true;
        for (int j = 0; sieve[j] * sieve[j] <= i; ++j) {
	    if (i % sieve[j] == 0) {
                flag = false;
                break;
            }
	}
        if (flag) sieve.push_back(i);
    }

    //dp
    int phi[n+1];
    memset(phi, 0, n+1);
    phi[1] = 0
    for (int p: sieve) {
        phi[p] = p - 1;
    }
    for (int i = 4; i <= n+1; ++i) {
        int a, b;
        for (int j = 0; sieve[j] * sieve[j] <= i; ++j) {
	    if (i % sieve[j] == 0) {
                a = sieve[j];
		b = i / sieve[j];
                break;
            }
        }
        phi[i] = phi[a] * phi[b];
    }      
}
