#include <list>
#include <stack>
#include <iostream>
#include <stdlib.h>

using namespace std;

// void find_scc(double *array, int *n, int *label);

class Graph
{
    list<int> *adj;
    void fillOrder(int v, bool visited[], stack<int> &Stack);
public:
    int V;
    int* label;
    Graph(int V);
    void addEdge(int v, int w);
    void getSccs(int *, int *);
    Graph getTranspose();
    void DFS(int v, bool visited[], int count);
};

Graph::Graph(int V)
{
    this->V = V;
    adj = new list<int>[V];
    this->label = (int*) malloc(V*sizeof(int));
}

void Graph::DFS(int v, bool visited[], int count)
{
    visited[v] = true;
    label[v] = count;
    list<int>:: iterator i;
    for(i = adj[v].begin(); i!=adj[v].end(); ++i)
        if(!visited[*i])
            DFS(*i, visited, count);
}

Graph Graph::getTranspose()
{
    Graph g(V);
    for (int v = 0; v<V; v++)
    {
        list<int>::iterator i;
        for(i = adj[v].begin(); i!=adj[v].end(); ++i)
        {
            g.adj[*i].push_back(v);
        }
    }
    return g;
}

void Graph::addEdge(int v, int w)
{
    adj[v].push_back(w);
}

void Graph::fillOrder(int v, bool visited[], stack<int> &stack){
    visited[v] = true;
    list<int>::iterator i;
    for(i = adj[v].begin(); i!=adj[v].end(); ++i)
        if(!visited[*i])
            fillOrder(*i, visited, stack);
    stack.push(v);
}

void Graph::getSccs(int *sccarray, int* numscc)
{
    stack<int> Stack;
    bool *visited = new bool[V];
    int count = 1;
    for (int i = 0; i<V; i++)
        visited[i] = false;
    for (int i = 0; i<V; i++)
        if(visited[i] == false)
            fillOrder(i, visited, Stack);
    Graph gr = getTranspose();
    for (int i =0; i<V; i++)
        visited[i] = false;
    while(Stack.empty() == false)
    {
        int v = Stack.top();
        Stack.pop();
        if (visited[v] == false)
        {
            gr.DFS(v, visited, count);
            count++;
        }
    }
    for (int i = 0 ; i<V; i++)
        sccarray[i] = gr.label[i];
    *numscc = count-1;
}

void realign(int num, int n, int* label)
{
    int maxcount = 0;
    int maxIndex = 0;
    for(int i = 0; i < num; i++)
    {
        int count = 0;
        for(int j = 0; j<n; j++)
        {
            if(label[j] == i+1) count++;
            if(count > maxcount)
            {
                maxcount = count;
                maxIndex = i+1;
            }
        }
    }
    if (maxIndex != 1) // swap with the maxindex with 1
    {
        for(int j = 0; j<n; j++){
            if(label[j] == 1) 
                label[j] = maxIndex;
            else if(label[j] == maxIndex) 
                label[j] = 1;
        }
    }
}

extern "C" void find_scc_(double *array, int *n, int *label, int *nums_scc){
    // array is the tcm
    // n is the number of states
    // label is the index of each states
    // nums_scc: how many scc in the tcm
    int num = *n;
    Graph g(num);
    for(int i = 0; i<num; ++i)
    {
        for (int j = 0 ; j<num; ++j)
        {
            if(array[i*num+j] > 0 && i!=j)
                g.addEdge(j,i);
        }
    }
    g.getSccs(label, nums_scc);
    realign(*nums_scc, *n, label);
}

extern "C" void find_undirected_scc_(double *array, int *n, int *label, int *nums_scc)
{
    int num = *n;
    Graph g(num);
    for (int i = 0; i<num; i++)
        for(int j = 0; j<num; j++)
            if(array[i*num+j] > 0.0 && i != j)
            {
                g.addEdge(i,j);
                g.addEdge(j,i);
            }
    stack<int> stack;
    bool *visited = new bool[g.V];
    for (int i = 0; i<g.V; i++)
        visited[i] = false;
    int count = 1;
    for(int i =0; i<g.V; i++){
        if(! visited[i])
            g.DFS(i, visited, count);
        count++;
    *nums_scc = count-1;
    for (int i = 0 ; i<g.V; i++)
        label[i] = g.label[i];
    }
    realign(*nums_scc, *n, label);
}

// int main()
// {
//     // double *a=(double*) malloc(16*sizeof(double));
//     double array[16] = {0.0,2.0,0.0,1.0,
//                         0.0,1.0,1.0,0.0,
//                         0.0,1.0,2.0,1.0,
//                         0.0,1.0,2.0,1.0};
    
//     int n = 4;
//     int num;
//     int *label = (int*) malloc(4*sizeof(int));
//     find_scc_(array, &n, label, &num);
//     for(int i = 0; i<4; i++)
//         cout << label[i] << " ";

//     cout << endl;
//     find_undirected_scc_(array, &n, label, &num);
//     for(int i = 0; i<4; i++)
//         cout << label[i] << " ";

// //     // int i = 0;
// //     // Graph g(5);
// //     // g.addEdge(1,0);
// //     // g.addEdge(0,2);
// //     // g.addEdge(2,1);
// //     // g.addEdge(0,3);
// //     // g.addEdge(3,4);
// //     // int *label = (int*) malloc(5*sizeof(int));
// //     // g.getSccs(label);

// //     // for(i=0; i<5; i++)
// //     //     cout << label[i] <<endl;
// //     // return 0;
// }
