#include <iostream>
#include <climits>
#include <cmath>
#include <iomanip>
#include <set>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include <utility>
#include <vector>
#include "gurobi_c++.h"
using namespace std;

struct comp 
{
  bool operator()(pair<int, int> const & a, pair<int, int> const & b) const 
  {
    if(a.first == b.first)  
    return a.second < b.second;
      
    return a.first > b.first;  
  }
};

int heur_upper_bnd = 0;
int NUM_THREADS = 4;
double TIME_LIMIT = 1800.0;

int conf;
bool dfs_val;

map<pair<int, int>, int> weights;
map<pair<int, int>, int> orig_weights;

unordered_map<int, int> freq;
set<pair<int, int>, comp> order;
vector<vector<int>> adj;
vector<vector<int>> orig_adj;

pair<int, int> best_edge_gl;
unordered_map<int, unordered_set<int>> freqEdges;

vector<bool> vis_res;

map<pair<int, int>, bool> visited_best_pairs;
vector<int> node_comp_idx;
vector<vector<int>> components;
vector<vector<int>> initialSol;

int n, m;
int best_reduction_val; 

int get_hash(pair<int, int> p1)
{
    if(p1.first == p1.second)
    return 0;

    int row = (p1.first < p1.second) ? p1.first : p1.second;
    int col = (p1.first < p1.second) ? p1.second : p1.first;
    
    return n*(row-1) + (col-1);

}

pair<int, int> get_edge(int val1)
{
    int row = val1/n;
    int col = val1 - row*n;
    
    return {row+1, col+1};

}

void dfs_res(int node, vector<int> &comp)
{
  vis_res[node] = 1;
  comp.push_back(node);

  for(int child: adj[node])
  {
    if(!vis_res[child])
    dfs_res(child, comp);
  }

}

void dfs_improve(int node, vector<int> &comp, int idx)
{
  vis_res[node] = 1;
  comp.push_back(node);
  node_comp_idx[node] = idx;

  for(int child: adj[node])
  {
    if(!vis_res[child])
    dfs_improve(child, comp, idx);
  }

}

void dfs_verify(int node, vector<int> &comp)
{
  vis_res[node] = 1;
  comp.push_back(node);

  for(int child: orig_adj[node])
  {
    if(!vis_res[child])
    dfs_verify(child, comp);
  }

}

void FindBest()
{
    for(int i=1; i<=n; ++i)
    { for(int j: adj[i])
      { for(int k: adj[j])
        {

          if(k <= i)
          continue;
          
          if(weights.find({i, k}) == weights.end())
          { 
              
            int val1 = get_hash({i, j});
            int val2 = get_hash({j, k});  
            int val3 = get_hash({i, k});
              
            freq[val1]++;
            freqEdges[val1].insert(val2);
            freqEdges[val2].insert(val1);
            
            freq[val2]++;
            freqEdges[val2].insert(val3);
            freqEdges[val3].insert(val2);
              
            
            freq[val3]++;
            freqEdges[val1].insert(val3);
            freqEdges[val3].insert(val1); 
             
            conf++; 
          }
        }
      }	
    }
 
   int max_freq = 0;
   int best_hash; 
    
   if(conf == 0)
   best_edge_gl = {-1, -1};   
   else 
   { for(pair<int, int> p: freq)
     { 
       order.insert({p.second, p.first});  
       
       if(p.second > max_freq)
       {  best_hash = p.first;
          max_freq = p.second;
       }
     }
    
     best_edge_gl = get_edge(best_hash);
   }
    
}

void ModifyEdge(pair<int, int> edge)
{   
    int hash_val_edge = get_hash(edge);
    
    order.erase({freq[hash_val_edge], hash_val_edge});
    
    freq[hash_val_edge] = 0;
    
    if(weights.find(edge) != weights.end())
    { weights.erase(edge);
      
      int a = edge.first, b = edge.second;
      
      for(int idx=0; idx<adj[a].size(); ++idx)
      { if(adj[a][idx] == b)
        { adj[a].erase(adj[a].begin() + idx);
          break;
        }
      }
     
      for(int idx=0; idx<adj[b].size(); ++idx)
      { if(adj[b][idx] == a)
        { adj[b].erase(adj[b].begin() + idx);
          break;
        }
      }
    }
    else
    { weights[edge] = 1;
      int a = edge.first, b = edge.second;
      adj[a].push_back(b);
      adj[b].push_back(a);
    }
    
    for(int idx1: freqEdges[hash_val_edge])
    {
      int cur_freq = freq[idx1];
      order.erase({cur_freq, idx1});
      if(cur_freq > 1) 
      order.insert({cur_freq - 1, idx1});  
          
      freq[idx1]--;
      freqEdges[idx1].erase(hash_val_edge);
    }
    
    best_edge_gl = {-1, -1};
    
    if(!order.empty())
    {   pair<int, int> front = *(order.begin());
    
        if(front.first > 0)
        best_edge_gl = get_edge(front.second); 
    }   
    
}

void ReduceConflicts()
{
    int prev_conf = INT_MAX;

    while(1)
    {
        conf = 0;
        FindBest();

        if(conf >= prev_conf)
        {
           vis_res.assign(n+1, 0);

           for(int i=1; i<=n; ++i)
           {
               if(!vis_res[i])
               {
                 vector<int> comp;
                 dfs_res(i, comp);

                 int sum_deg = 0;
                 int sz = comp.size();

                 for(int nd: comp)
                 sum_deg += adj[nd].size();

                 if(sum_deg != (sz*(sz-1)))
                 { int rem = sum_deg/2;
                   int add = (sz*(sz-1) - sum_deg)/2;

                   if(rem <= add)
                   {
                     for(int j=0; j<comp.size(); ++j)
                     {
                       for(int k=(j+1); k<comp.size(); ++k)
                       {
                         int node1 = comp[j];
                         int node2 = comp[k];

                         if(node1 > node2)
                         swap(node1, node2);

                         if(weights.find({node1, node2}) != weights.end())
                         weights.erase({node1, node2});
                       }
                     }
                    }
                    else
                    {
                      for(int j=0; j<comp.size(); ++j)
                      { 
                        for(int k=(j+1); k<comp.size(); ++k)
                        {
                          int node1 = comp[j];
                          int node2 = comp[k];

                          if(node1 > node2)
                          swap(node1, node2);

                          if(weights.find({node1, node2}) == weights.end())
                          weights[{node1, node2}] = 1;
                        }
                      }
                     }
                   }
                 }
            }

           break;
        }
        else
        prev_conf = conf;
        
        if(best_edge_gl == make_pair(-1, -1))
        break;

        while(1)
        { if(best_edge_gl == make_pair(-1, -1))
          break;

          ModifyEdge(best_edge_gl);
        }

        freq.clear();
        freqEdges.clear();
        order.clear();
    }

}

void IndexComponents()
{
    vis_res.assign(n+1, 0);
    node_comp_idx.assign(n+1, -1);
    
    int comp_idx = -1;

    for(int i=1; i<=n; ++i)
    {
      if(!vis_res[i])
      {
        vector<int> comp;
        comp_idx++; 
        dfs_improve(i, comp, comp_idx);

        components.push_back(comp);
      }
    }

}

int CalculateNodeRemovalCost(int node, int node_idx)
{
    int removal_cost = 0;

    for(int comp_node: components[node_idx])
    {
      if(comp_node == node)
      continue;

      int nd1 = node;
      int nd2 = comp_node;

      if(nd1 > nd2)
      swap(nd1, nd2);  

      if(orig_weights.find({nd1, nd2}) == orig_weights.end())
      removal_cost++;
      else 
      removal_cost--;
    }

    return removal_cost;

}

int CalculateNodeAdditionCost(int node, int comp_idx)
{
    int addition_cost = 0;

    for(int comp_node: components[comp_idx])
    {
      int nd1 = node;
      int nd2 = comp_node;

      if(nd1 > nd2)
      swap(nd1, nd2);  

      if(orig_weights.find({nd1, nd2}) != orig_weights.end())
      addition_cost++;
      else 
      addition_cost--;
    }

    return addition_cost;
}

void RemoveNodeComp(int best_node)
{
    for(int nd: adj[best_node])
    {
      if(best_node < nd)
      weights.erase({best_node, nd});
      else
      weights.erase({nd, best_node}); 

      for(int x=0; x<adj[nd].size(); ++x)
      { if(adj[nd][x] == best_node)
        {adj[nd].erase(adj[nd].begin() + x); 
         break;
        }
      }
    }

    adj[best_node].clear();

}

void AddNodeComp(int best_node, int best_comp)
{
    for(int nd: components[best_comp])
    { 
      adj[best_node].push_back(nd);

      if(best_node < nd)
      weights[{best_node, nd}] = 1;
      else
      weights[{nd, best_node}] = 1;

      adj[nd].push_back(best_node);
    }  
}

bool BestFitImprove()
{
  if(!dfs_val)
  IndexComponents();

  best_reduction_val = 0;
  int best_node = -1;
  int best_comp = -1;

  for(int i=1; i<=n; ++i)
  { 
    int curr_idx = node_comp_idx[i];

    int removal_cost = CalculateNodeRemovalCost(i, curr_idx);

    for(int j=0; j<components.size(); ++j)
    {
        if(j == curr_idx)
        continue;

        int addition_cost = CalculateNodeAdditionCost(i, j);

        int total_cost = addition_cost + removal_cost;

        if(total_cost >= removal_cost && total_cost >= best_reduction_val)
        {
          if((visited_best_pairs.find({i, j}) == visited_best_pairs.end() || visited_best_pairs.find({i, curr_idx}) == visited_best_pairs.end()))
          {
             best_reduction_val = addition_cost + removal_cost;
             best_node = i;
             best_comp = j;
          }
        }
        else if(removal_cost > best_reduction_val)
        {
          if((visited_best_pairs.find({i, -1}) == visited_best_pairs.end() || visited_best_pairs.find({i, curr_idx}) == visited_best_pairs.end()))
          {
             best_reduction_val = removal_cost;
             best_node = i;
             best_comp = -1;
          }
        }
    }
  }

  if(best_node == -1)
  return 0;

  RemoveNodeComp(best_node);

  if(best_comp != -1)
  AddNodeComp(best_node, best_comp);
  
  visited_best_pairs[{best_node, best_comp}] = 1; 
  visited_best_pairs[{best_node, node_comp_idx[best_node]}] = 1;

  for(int nd_idx = 0; nd_idx < components[node_comp_idx[best_node]].size(); ++nd_idx)
  {
    if(components[node_comp_idx[best_node]][nd_idx] == best_node)
    {
      components[node_comp_idx[best_node]].erase(components[node_comp_idx[best_node]].begin() + nd_idx);
      break;
    }
  } 

  if(best_comp != -1)
  { components[best_comp].push_back(best_node);
    node_comp_idx[best_node] = best_comp;
  }
  else
  {
    components.push_back({best_node});
    node_comp_idx[best_node] = components.size() - 1;
  }
  
  return 1;

}

void ReadInput()
{
    string s1, s2; cin>>s1>>s2;
    cin>>n>>m;

    adj.resize(n+1);
    orig_adj.resize(n+1);
    
    for(int i=1; i<=m; ++i)
    {
        int a, b; cin>>a>>b;

        adj[a].push_back(b); orig_adj[a].push_back(b);
        adj[b].push_back(a); orig_adj[b].push_back(a);
        
        if(a > b)
        swap(a, b);
            
        weights[{a, b}] = 1;
        orig_weights[{a, b}] = 1;
    }

}

void CalculateFinalModifications()
{

    for(int i=1; i<=n; ++i)
    { for(int j=(i+1); j<=n; ++j) 
      {
        if((weights.find({i, j}) != weights.end()) != (orig_weights.find({i, j}) != orig_weights.end())) 
        heur_upper_bnd++; 

        if(weights.find({i, j}) != weights.end())
        initialSol.push_back({i, j, 1});
        else
        initialSol.push_back({i, j, 0});  
      }
    }

}

void AddNeighbour(pair<int, int> edge)
{
   int node = edge.first;
   int nghr = edge.second;
  
   orig_adj[node].push_back(nghr);
   orig_adj[nghr].push_back(node);
}

void RemoveNeighbour(pair<int, int> edge)
{
   int node = edge.first;
   int nghr = edge.second;

   for(int idx = 0; idx < orig_adj[node].size(); ++idx)
   {  
      if(orig_adj[node][idx] == nghr)
      { orig_adj[node].erase(orig_adj[node].begin() + idx);
        break;
      }
   }

   for(int idx = 0; idx < orig_adj[nghr].size(); ++idx)
   { 
      if(orig_adj[nghr][idx] == node)
      { orig_adj[nghr].erase(orig_adj[nghr].begin() + idx);
        break;
      }
   }
}

bool VerifySolution()
{
  vis_res.assign(n + 1, false);
  
  for(int i=1; i<=n; ++i)
  { 
     if(!vis_res[i])
     {
        vector<int> comp;
        dfs_verify(i, comp);
        
        int num_nodes = comp.size();
        int nghrs = 0;
 
        for(int node: comp)
        nghrs += orig_adj[node].size();

        if(nghrs != (num_nodes*(num_nodes - 1)))
        return false; 
     }
  }

  return true;
}
  
void MIPSolverGurobi()
{
  try 
  {

    GRBEnv env = GRBEnv(true);
    env.start();

    GRBModel ClusterEditing = GRBModel(env);

    map<pair<int, int>, GRBVar> possible_edges;
    GRBLinExpr objective = GRBLinExpr();

    for(vector<int> var: initialSol)
    {
      int a = var[0], b = var[1];
      double val = var[2];

      GRBVar x = ClusterEditing.addVar(0.0, 1.0, 0.0, GRB_BINARY);
      x.set(GRB_DoubleAttr_Start, val);
      possible_edges[{a, b}] = x;

      if(orig_weights.find({a, b}) != orig_weights.end())
      objective += (1 - x);
      else
      objective += x;  
    }

    ClusterEditing.setObjective(objective, GRB_MINIMIZE);

    for(int u = 1; u<=n; ++u)
    { for(int v = (u + 1); v<=n; ++v)
      { for(int w = (v + 1); w<=n; ++w)
          { 
            if(orig_weights.find({u, v}) != orig_weights.end() and orig_weights.find({v, w}) != orig_weights.end() and orig_weights.find({u, w}) == orig_weights.end())
            ClusterEditing.addConstr(possible_edges[{u, v}] + possible_edges[{v, w}] - possible_edges[{u, w}] <= 1);
            else if(orig_weights.find({u, v}) != orig_weights.end() and orig_weights.find({v, w}) == orig_weights.end() and orig_weights.find({u, w}) != orig_weights.end())
            ClusterEditing.addConstr(possible_edges[{u, v}] - possible_edges[{v, w}] + possible_edges[{u, w}] <= 1);
            else if(orig_weights.find({u, v}) == orig_weights.end() and orig_weights.find({v, w}) != orig_weights.end() and orig_weights.find({u, w}) != orig_weights.end())
            ClusterEditing.addConstr(-possible_edges[{u, v}] + possible_edges[{v, w}] + possible_edges[{u, w}]  <= 1);
          }
      }
    }

    // Optimize model
    ClusterEditing.set(GRB_DoubleParam_TimeLimit, TIME_LIMIT);
    ClusterEditing.set(GRB_IntParam_Threads, NUM_THREADS);
    ClusterEditing.set(GRB_IntParam_NodeMethod, 2);

    ClusterEditing.optimize();

    int status = ClusterEditing.get(GRB_IntAttr_Status);

    vector<string> status_codes = {"LOADED", "OPTIMAL", "INFEASIBLE", "INF_OR_UNBD", "UNBOUNDED", "CUTOFF", "ITERATION_LIMIT", "NODE_LIMIT", "TIME_LIMIT", "SOLUTION_LIMIT", "INTERRUPTED", "NUMERIC", "SUBOPTIMAL", "INPROGRESS", "USER_OBJ_LIMIT", "WORK_LIMIT"};

    if((status == GRB_INF_OR_UNBD) || (status == GRB_INFEASIBLE) || (status == GRB_UNBOUNDED))
    { cout << "\nThe Cluster Editing model cannot be solved because it is infeasible or unbounded" << endl;
      return;
    }

    if(status != GRB_OPTIMAL)
    { cout << "\nOptimization was stopped with status " << status_codes[status - 1] << endl;
      return;
    }
     
    int modif = 0;
    vector<pair<int, int>> modifications;  
      
    for(pair<pair<int, int>, GRBVar> edge: possible_edges)
    {
         int val = round(edge.second.get(GRB_DoubleAttr_X));
         pair<int, int> curr_edge = edge.first;

         if(val == 1 && orig_weights.find(curr_edge) == orig_weights.end())
         {  
             modifications.push_back(curr_edge);
             AddNeighbour(curr_edge); 
             modif++;
         }
         else if(val == 0 && orig_weights.find(curr_edge) != orig_weights.end())
         {
             modifications.push_back(curr_edge);
             RemoveNeighbour(curr_edge);        
             modif++;
         }
    }
   
    cout<<"\nVerifying Solution....\n";

    if(VerifySolution())
    {
      cout << "\nEdge Modifications: \n";
     
      for(pair<int, int> edge: modifications)
      cout<<edge.first<<" "<<edge.second<<"\n";
 
      cout << "\nTotal Number of Modifications: " << modif << endl;
    }
    else
    cout << "\nUnable to find an OPTIMAL SOLUTION" << endl;  

  } 

  catch(GRBException e) 
  {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  }

  catch(...) 
  {
    cout << "Exception during optimization" << endl;
  }

}

int main(int argc, char* argv[])
{
    time_t start, heur_end, end;

    if(argc == 2)
    {
       string parameter = argv[1];
       size_t len = parameter.size();

       size_t equal_to_idx = parameter.find_first_of("=");
       
       string flag = parameter.substr(0, equal_to_idx);
       
       int val = stoi(parameter.substr(equal_to_idx + 1, len - equal_to_idx));
  
       if(flag == "--time-limit") 
       TIME_LIMIT = double(val);
       else if(flag == "--num-threads")
       NUM_THREADS = val;
    }
    else if(argc == 3)
    {
       string parameter = argv[1];
       size_t len = parameter.size();

       size_t equal_to_idx = parameter.find_first_of("=");

       string flag = parameter.substr(0, equal_to_idx);

       int val = stoi(parameter.substr(equal_to_idx + 1, len - equal_to_idx));

       if(flag == "--time-limit")
       TIME_LIMIT = double(val);
       else if(flag == "--num-threads")
       NUM_THREADS = val;

       parameter = argv[2];
       len = parameter.size();

       equal_to_idx = parameter.find_first_of("=");

       flag = parameter.substr(0, equal_to_idx);

       val = stoi(parameter.substr(equal_to_idx + 1, len - equal_to_idx));

       if(flag == "--time-limit")
       TIME_LIMIT = double(val);
       else if(flag == "--num-threads")
       NUM_THREADS = val;       
    }            
    
    ReadInput();
   
    time(&start);
    ReduceConflicts(); 

    while(BestFitImprove())
    {
      if(!dfs_val)
      dfs_val = 1;
    }

    for(int i=1; i<=n; ++i)
    { 
      int curr_idx = node_comp_idx[i];

      int removal_cost = CalculateNodeRemovalCost(i, curr_idx);

      int total_cost = INT_MAX, best_comp = -1;

      for(int j=0; j<components.size(); ++j)
      {
          if(j == curr_idx)
          continue;

          int addition_cost = CalculateNodeAdditionCost(i, j);
          int tot_cost = addition_cost + removal_cost;

          if(tot_cost <= 0)
          {  if(abs(tot_cost) < total_cost)
             { total_cost = abs(tot_cost);
               best_comp = j;
             }
          } 
          else
          total_cost = INT_MIN;  
      }
    }
	    
    CalculateFinalModifications();
   
    time(&heur_end);
    double time_elapsed = double(heur_end - start);
    TIME_LIMIT -= time_elapsed;

    cout<< "\nHeuristic Upper Bound: "<<heur_upper_bnd<<"\n";
    cout << "Time Elapsed in obtaining Heuristic Upper Bound: "<<fixed<<setprecision(5)<<time_elapsed<<" seconds\n\n";

    MIPSolverGurobi();

    time(&end);
    double time_taken = double(end - start);
    cout<<"\nTotal Time Taken: "<<fixed<<setprecision(5)<<time_taken<<" seconds\n\n";
    
    return 0;
}
