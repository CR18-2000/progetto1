#include <iostream>
#include <queue>
#include <iterator>
#include <list>
#include <fstream>
#include <algorithm>

using namespace std;

# define INF 0x3f3f3f3f        //infinito per il primo Dijkstra
# define INF2 999999           //infinito per settare il peso degli arichi che partono dal nodo che togliamo
typedef pair<int, int> iPair;
vector<vector<int> > paths;    //percorsi ottenuti con Dijkstra del grafo dato in input da ogni nodo a Powarts
//vector<vector<int> > paths2;   //percorsi ottenuti con Dijkstra dal grafo a cui viene tolta una citta valutata come quella da bloccare
vector<int> dist;              //vettore delle distanze minime di ogni nodo da Powarts nel grafo di partenza
vector<int> dist2;             //vettore delle distanze minime di ogni nodo da Powarts nel grafo a cui viene tolta la citta da bloccare
int N;                         //numero di nodi del grafo

struct Graph
{
	//int V;
	list< pair<int, int> > *adj;
};

//algoritmo di Dijkstra che calcola i percorsi minimi da ogni nodo a Powarts(src)
void  shortestPath(Graph g,int src)
{
	priority_queue< iPair, vector <iPair> , greater<iPair> > pq;
	//priority_queue<int> pq;
	pq.push(make_pair(0, src));     //inserimento del nodo di partenza (che ha distanza 0 da sè stesso)
	//pq.push(src);
	dist[src] = 0;

	vector<bool> f(N, false);       //nessun nodo è ancora stato visitato completamente
	while (!pq.empty())
	{
		int u = pq.top().second;
		//int u=pq.top();
		pq.pop();
		f[u] = true;
		list< pair<int, int> >::iterator i;
		for (i = g.adj[u].begin(); i != g.adj[u].end(); ++i)    //aggiunta dei nodi vicini a quello che stiamo visitando alla coda
		{
			int v = (*i).first;
			int weight = (*i).second;
			if (f[v] == false && dist[v] > dist[u] + weight)
			{
                paths[v].push_back(u);
				dist[v] = dist[u] + weight;
				pq.push(make_pair(dist[v], v));
				//pq.push(v);
			}
		}
	}
/*

	printf("Vertex Distance from Source\n");
	for (int i = 0; i < N; ++i)
		printf("%d \t\t %d\n", i, dist[i]);

*/
}

void  shortestPath1(Graph g,int src)
{
	priority_queue< iPair, vector <iPair> , greater<iPair> > pq;
	//priority_queue<int> pq;
	pq.push(make_pair(0, src));
	//pq.push(0);
	dist2[src] = 0;

	vector<bool> f(N, false);
	while (!pq.empty())
	{
		int u = pq.top().second;
		//int u = pq.top();
		pq.pop();
		f[u] = true;
		list< pair<int, int> >::iterator i;
		for (i = g.adj[u].begin(); i != g.adj[u].end(); ++i)
		{
			int v = (*i).first;
			int weight = (*i).second;
			if (f[v] == false && dist2[v] > dist2[u] + weight)
			{
                //paths2[v].push_back(u);
                dist2[v] = dist2[u] + weight;
				pq.push(make_pair(dist2[v], v));
				//pq.push(v);
			}
		}
	}
/*

	printf("Vertex Distance from Source\n");
	for (int i = 0; i < N; ++i)
		printf("%d \t\t %d\n", i, dist[i]);

*/
}





/*bool trova (int p,int i){

  for (int j = 0; j < paths[i].size(); j++) {
    if (paths[i][j]==p) {
      return true;
    }
  }

return false;

}*/
//9-11  15-16-17-18-19

int main()
{
    ifstream in ("input.txt");
    //ifstream in2 ("input19.txt");
    int M, P, t1, t2, t3;   //N = numero nodi, M = numero archi, P = nodo di partenza
    in >> N >> M >> P;

    paths.resize(N);
    //paths2.resize(N);
	Graph g;
	//Graph g2;
	dist.resize(N,INF);     //setto tutte le distanze tra i nodi a infinito
	//dist2.resize(N,INF);    //setto tutte le distanze tra i nodi a infinito
    //g.V = N;
	//g2.V=N;
	g.adj = new list<iPair> [N];
	//g2.adj = new list<iPair> [N];

  for (int i = 0; i < M; i++) {
    in >> t1 >> t2 >> t3;
    g.adj[t1].push_back(make_pair(t2, t3));
    g.adj[t2].push_back(make_pair(t1, t3));
  }
  in.close();
	//in.clear();
	shortestPath(g,P); // log(n)


//fix paths
  for (int i = 0; i < N; i++) { //O(n^3)
		int index_next_path;
		if (i!=P) {
            while(find(paths[i].begin(),paths[i].end(),P)==paths[i].end()){
                index_next_path=paths[i].back();
                for (int j = 0; j < paths[index_next_path].size(); j++) {
                  paths[i].push_back(paths[index_next_path][j]);
                }
            }
		}
  }

	// CONTO LE OCCORRENZE g
	vector<int> occ(N,0);
	for (int i = 0; i < N; i++) { // O(n^2)
		//cout << endl << "path "<< i << " : ";
		for (int j = 0; j < paths[i].size(); j++) {
			occ[paths[i][j]]++;
		}
	}

	//TROVO LA CITTA DA TOGLIERE
    /*int citta = -1;
    int max = -1;
	for (int i = 0; i < occ.size(); i++) { //O(n)
		if (occ[i]>max && i!=P){
			 citta = i;
			 max = occ[i];
		}
	}*/

	occ[P]=-1;
	int citta=0;
	vector<int> sol;
    sol.resize(N);
    sol.clear();
    int max = -1;
    vector<int> sol2(N);
    Graph g2;

	//bool fine = false;
	while (citta > -1) {
        //paths2.resize(N);
        dist2.resize(N,INF);    //setto tutte le distanze tra i nodi a infinito
        citta = -1;
        //int max = -1;
        for (int i = 0; i < occ.size(); i++) { //O(n)
            if (occ[i]>max){
                citta = i;
                max = occ[i];
            }
        }
        occ[citta]=-1;
        //cout<<citta<<endl;
        //ifstream in2 ("input.txt");
        sol2.clear();
        dist2.clear();
        //ifstream in2 ("input0.txt");
        in.open("input0.txt");
        in >> N >> M >> P;
        g2.adj = new list<iPair> [N];
        for (int i = 0; i < M; i++) {
            in >> t1 >> t2 >> t3;
            if (t1 == citta || t2 == citta) {
                g2.adj[t1].push_back(make_pair(t2, INF2));
                g2.adj[t2].push_back(make_pair(t1, INF2));
            } else {
                g2.adj[t1].push_back(make_pair(t2, t3));
                g2.adj[t2].push_back(make_pair(t1, t3));

            }
        }
        in.close();
        shortestPath1(g2,P);
        for (int i = 0; i < N; i++) {
            if (dist[i]!=dist2[i]){
                sol2.push_back(i);
            }
        }
        //cout <<"massimo : "<< max << endl;
        if (sol2.size()>sol.size()) {
            sol.clear();
            for (int i=0; i<sol2.size(); i++) {
                sol.push_back(sol2[i]);
                //sol[i]=sol2[i];
                //cout<<sol[i]<<endl;
            }
        }
        max = sol.size();
	}

    /*in2 >> N >> M >> P;
	for (int i = 0; i < M; i++) {
		in2 >> t1 >> t2 >> t3;
		if (t1 == citta || t2 == citta) {
			g2.adj[t1].push_back(make_pair(t2, INF2));
			g2.adj[t2].push_back(make_pair(t1, INF2));
		} else {
			g2.adj[t1].push_back(make_pair(t2, t3));
			g2.adj[t2].push_back(make_pair(t1, t3));

		}
	}

	shortestPath1(g2,P);*/

    /*vector<int> sol;
    sol.resize(N);
    sol.clear();*/

	/*for (int i = 0; i < N; i++) {
		if (dist[i]!=dist2[i]){
			sol.push_back(i);
		}
	}

    vector<int> sol2(N);
    for (int i=0; i<altriMax.size(); i++) {
        ifstream in2 ("input.txt");
        sol2.clear();
        in2 >> N >> M >> P;
        Graph g2;
        for (int i = 0; i < M; i++) {
            in2 >> t1 >> t2 >> t3;
            if (t1 == altriMax[i] || t2 == altriMax[i]) {
                g2.adj[t1].push_back(make_pair(t2, INF2));
                g2.adj[t2].push_back(make_pair(t1, INF2));
            } else {
                g2.adj[t1].push_back(make_pair(t2, t3));
                g2.adj[t2].push_back(make_pair(t1, t3));

            }
        }
        shortestPath1(g2,P);
        for (int j = 0; j < N; j++) {
            if (dist[j]!=dist2[j]){
                sol2.push_back(j);
            }
        }
        if (sol2.size()>sol.size()) {
            for (int j=0; j<sol2.size(); j++) {
                sol[j]=sol2[j];
            }
        }
    }*/

   ofstream out ("output.txt");
   out << sol.size() <<endl;
   cout << sol.size() << endl;
   for (int i=0; i<sol.size(); i++) {
       out << sol[i] << endl;
   }

	return 0;
}
