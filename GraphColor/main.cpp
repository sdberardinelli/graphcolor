/*******************************************************************************
 * Filename      : main.cpp
 * Header File(s): 
 * Description   :
 * Authors(s)    : 
 * Date Created  : 
 * Date Modified :
 * Modifier(s)   :
 *******************************************************************************/
/************************************
 * Included Headers 
 ************************************/
#include <iostream>
#include <lpsolve/lp_lib.h>
#include <random>
#include <vector>
#include <string>
#include <algorithm>
#include <sstream>
#include <boost/graph/undirected_graph.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/lexical_cast.hpp>
#include <ctime>

/************************************
 * Namespaces 
 ************************************/
using namespace std;

/************************************
 * Local Types 
 ************************************/
//https://github.com/codebrainz/color-names/


const char* color_data[] = {
"aliceblue",      "antiquewhite","antiquewhite1","antiquewhite2","antiquewhite3",
"antiquewhite4",  "aquamarine","aquamarine1","aquamarine2","aquamarine3",
"aquamarine4",    "azure","azure1","azure2","azure3",
"azure4",         "beige","bisque","bisque1","bisque2",
"bisque3",        "bisque4","black","blanchedalmond","blue",
"blue1",          "blue2","blue3","blue4","blueviolet",
"brown",          "brown1","brown2","brown3","brown4",
"burlywood",      "burlywood1","burlywood2","burlywood3","burlywood4",
"cadetblue",      "cadetblue1","cadetblue2","cadetblue3","cadetblue4",
"chartreuse",     "chartreuse1","chartreuse2","chartreuse3","chartreuse4",
"chocolate",      "chocolate1","chocolate2","chocolate3","chocolate4",
"coral",          "coral1","coral2","coral3","coral4"
"cornflowerblue", "cornsilk","cornsilk1","cornsilk2","cornsilk3",
"cornsilk4",      "crimson","cyan","cyan1","cyan2",
"cyan3",          "cyan4","darkgoldenrod","darkgoldenrod1","darkgoldenrod2",
"darkgoldenrod3", "darkgoldenrod4","darkgreen","darkkhaki","darkolivegreen",
"darkolivegreen1","darkolivegreen2","darkolivegreen3","darkolivegreen4","darkorange",
"darkorange1","	darkorange2","darkorange3","darkorange4","darkorchid",
"darkorchid1","	darkorchid2","darkorchid3","darkorchid4","darksalmon",
"darkseagreen","darkseagreen1","darkseagreen2","darkseagreen3","darkseagreen4",
"darkslateblue","darkslategray","darkslategray1	darkslategray2	darkslategray3",
"darkslategray4","darkslategrey","darkturquoise","darkviolet","deeppink",
"deeppink1","deeppink2","deeppink3","deeppink4","deepskyblue",
"deepskyblue1","deepskyblue2","deepskyblue3","deepskyblue4","dimgray",
"dimgrey","dodgerblue","dodgerblue1","dodgerblue2","dodgerblue3",
"dodgerblue4","	firebrick","firebrick1","firebrick2","firebrick3",
"firebrick4","floralwhite","forestgreen","gainsboro","ghostwhite",
"gold","gold1","gold2","gold3","gold4",
"goldenrod","goldenrod1","goldenrod2","goldenrod3","goldenrod4"
};
/************************************
 * Local Variables 
 ************************************/
bool still_running;
int total_nodes = 3;

/************************************
 * Local Functions 
 ************************************/

vector< vector<int> > GenerateGraph ( int );
void GraphDisplay ( vector< vector<int> > );
vector<string> using_lpsolve ( vector< vector<int> > );
void GraphOutput ( vector< vector<int> >, vector<string> );

/*******************************************************************************
* Function     : 
* Description  : 
* Arguments    : 
* Returns      : 
* Remarks      : 
********************************************************************************/
int main ( int argc, char* argv[] ) 
{
    vector< vector<int> > g;
    vector<string> colors;
    still_running = true;
    do
    {
        g = GenerateGraph(total_nodes);
        

        colors = using_lpsolve(g);

        GraphOutput(g,colors);
        total_nodes++;
        for ( int i = 0; i < g.size(); i++ )
            g[i].clear();
        g.clear();
        colors.clear();
    }while(still_running);
}
/*******************************************************************************
* Function     : 
* Description  : 
* Arguments    : 
* Returns      : 
* Remarks      : 
********************************************************************************/
vector<string> using_lpsolve ( vector< vector<int> > g )
{  
    vector<string> colors;
    int     Ncol = g.size();
    int     Mcol = g.size();

    lprec  *lp;
    lp = make_lp(0, Ncol*Mcol);
    
    vector<std::vector<int> > x(Ncol,vector<int>(Mcol));
    for (int v = 0, l = 1; v < Ncol; ++v) 
    {
        for (int k = 0; k < Mcol; ++k, l++ )
        {
            x[v][k] = l;
        }    
    }    
    
    {
        int rowno[1];
        double col[1];
        int idx = 0;
        rowno[idx] = idx+1;
        col[idx] = Mcol;                
        idx++; 
        set_int(lp,idx,TRUE);
        add_columnex(lp, idx, col, rowno);
    }

    set_add_rowmode(lp, TRUE);      
    
    {   // x_v1+x_v2+...+x_vK = 1, v in V
        for ( int ii = 0; ii < Ncol; ii++ )
        {
            int colno[Mcol*Ncol+1];
            double row[Mcol*Ncol+1];
            int idx = 0;  
            colno[idx] = idx+1;
            row[idx] = 0;
            idx++;  
            for ( int i = 0; i < Ncol; i++ )
            {                   
                for ( int j = 0; j < Mcol; j++ )
                {
                    colno[idx] = idx+1;
                    if ( ii == i )
                    {
                        row[idx] = 1;
                    }
                    else
                    {
                        row[idx] = 0;
                    }
                    set_int(lp,idx,TRUE);
                    idx++;
                }
            }
            add_constraintex(lp, idx, row, colno, EQ, 1);
        }
    } 
    
   {   // y >= k*x_vk, v in V, k=1...K
        for ( int ii = 0; ii < Ncol; ii++ )
        {        
            for ( int jj = 0; jj < Mcol; jj++ )
            {
                int colno[Ncol*Mcol+1];
                double row[Ncol*Mcol+1];
                int idx = 0;  
                colno[idx] = idx+1;
                row[idx] = -1;
                idx++;
                for ( int i = 0; i < Ncol; i++ )
                {              
                    for ( int j = 0; j < Mcol; j++ )
                    {
                        colno[idx] = idx+1;
                        if ( ii == i && jj == j )
                            row[idx] = j+1;
                        else
                            row[idx] = 0;
                        set_int(lp,idx,TRUE);
                        idx++;
                    }
                }
                add_constraintex(lp, idx, row, colno, LE, 0);     
            }
        }
    }    
    
    {   // 1 >= x_u,k + x_v,k u,v in E, k=1...K    
        for ( int src = 0; src < Ncol; src++ )
        {
            const vector<int>& succs = g[src];
            for ( int s = 0; s < succs.size(); s++ )
            {            
                int dst = succs[s];
                if (src > dst)
                {
                    for (int k = 0; k < Mcol; k++)
                    {                    
                        int colno[Ncol*Mcol+1];
                        double row[Ncol*Mcol*+1];
                        int idx = 0;    
                        colno[idx] = idx+1;
                        row[idx] = 0;
                        idx++;           
                        for ( int i = 0; i < Ncol; i++ )
                        {                          
                            for (int j = 0; j < Mcol; j++)
                            {
                               colno[idx] = idx+1;
                               row[idx] = 0;
                               set_int(lp,idx,TRUE);
                               idx++;
                            }
                        }
                        
                        row[x[src][k]] = 1;
                        row[x[dst][k]] = 1;
                        add_constraintex(lp, idx, row, colno, LE, 1);                        
                    }
                }
            }
        }
    } 
 
    set_add_rowmode(lp, FALSE);
    
    {   // min: y
        int colno[1];
        double row[1];
        int idx = 0;
        colno[idx] = idx+1;
        row[idx] = 1;
        idx++;
        set_obj_fnex(lp, idx, row, colno);
    }      
    
    set_verbose(lp, IMPORTANT);

    set_minim(lp);
    
    cout << "running lp solve" << endl;
    time_t start = time(0);
    if ( solve(lp) == OPTIMAL )
    {
        cout << "done lp solve with optimal" << endl;
        double seconds_since_start = difftime( time(0), start);
        
        cout << get_objective(lp) << "," << seconds_since_start << "," << total_nodes << endl;
        
//        cout << "objective value " << get_objective(lp) << endl;
        {
            double row[Ncol*Mcol+1];
            get_variables(lp, row);
            int idx = 0, i = 0, j = 0;
            for ( i = 0; i < Ncol; i++ )
            {
                for ( j = 0; j < Mcol ; j++ )
                {
                    if ( row[idx] == 1 )
                    {
//                        cout << "node " << i << " gets color " << color_data[(j*5)%150] << endl;
                        string str(color_data[(j*5)%150]);
                        colors.push_back(str);
                    }
                    idx++;
                }
            }
            if ( row[idx] )
            {
//                cout << "node " << i << " gets color " << color_data[(j*5)%150] << endl;
                string str(color_data[(j*5)%150]);
                colors.push_back(str);
            }
        }
        
    }
    else
    {
        cout << "done lp solve without optimal" << endl;
    }
    delete_lp(lp); 
    
    return colors;
}
/*******************************************************************************
* Function     : 
* Description  : 
* Arguments    : 
* Returns      : 
* Remarks      : 
********************************************************************************/
vector< vector<int> > GenerateGraph ( int n )
{
    vector< vector<int> > graph;
    random_device rd;
    mt19937 rnd(rd());
    uniform_int_distribution<int> dist(0,n-1);
    
    for ( vector<int>::size_type i = 0; i < n; i++ )
    {
        vector<int> tmp;
        tmp.push_back(i);
        graph.push_back(tmp);
    }
    for ( vector<int>::size_type i = 0; i < graph.size(); i++ )
    {
        for ( vector<int>::size_type j = 0; j < n; j++ )
        {
            int pick = dist(rnd);
            if ( i == j || pick == j )
                continue;
            
            graph[i].push_back(pick);
            graph[pick].push_back(i);

        }  
    }
    
    for ( vector<int>::size_type i = 0; i < graph.size(); i++ )
    {
        for ( vector<int>::size_type k = 1; k < graph[i].size(); k++ )
        {
            if ( graph[i][k] == i )
            {
                graph[i].erase(graph[i].begin()+k);
                k--;
            }
        }
        
        sort( graph[i].begin()+1, graph[i].end() );
        for ( vector<int>::size_type k = 1; k < graph[i].size()-1; k++ )
        {
            if ( graph[i][k]==graph[i][k+1] )
            {
                graph[i].erase(graph[i].begin()+k+1);
                k--;
            }
        }        
    }
    return graph;
}
/*******************************************************************************
* Function     : 
* Description  : 
* Arguments    : 
* Returns      : 
* Remarks      : 
********************************************************************************/
void GraphDisplay ( vector< vector<int> > g )
{
    for (int i = 0; i < g.size(); ++i)
    {
        for (int j = 0; j < g[i].size(); ++j)
        {
            cout << g[i][j] << " ";
        }
        cout << endl;
    }
}
/*******************************************************************************
* Function     : 
* Description  : 
* Arguments    : 
* Returns      : 
* Remarks      : 
********************************************************************************/
void GraphOutput ( vector< vector<int> > graph, vector<string> colors )
{
    ofstream fout("graph.dot");

    typedef boost::undirected_graph<boost::no_property> Graph;
    
    Graph g;
        
    vector<Graph::vertex_descriptor> vertices;
    for ( vector<int>::size_type i = 0; i < graph.size();  i++ )
    {
        Graph::vertex_descriptor v = g.add_vertex();
        vertices.push_back(v);
    }
    
    bool add;
    for ( vector<int>::size_type i = 0; i < graph.size();  i++ )
    {
        for ( vector<int>::size_type j = 0; j < graph[i].size();  j++ )
        {
            if ( i == graph[i][j] )
                continue;
            
            add = true;
            
            if ( graph[i][j] < i )
            {
                for ( vector<int>::size_type k = 0; k < graph[graph[i][j]].size();  k++ )
                {
                    if ( graph[graph[i][j]][k] == i )
                    {
                        add = false;
                    }
                }
            }
            
            if ( add )
            {
                boost::add_edge(vertices[i],vertices[graph[i][j]],g);
            }
        }
    }     

    boost::write_graphviz(fout,g);
    fout.close();
    ifstream dot1("graph.dot");
    ofstream dot2("graphc.dot");
    string line,str;
    stringstream ss;
    getline(dot1, line);
    dot2 << line << endl;
    for ( int i = 0; i < graph.size(); i++ )
    {
        getline(dot1, line);
        unsigned int node = boost::lexical_cast<unsigned int>(line.at(0));
        ss << node;
        dot2 << ss.str() << "[color = " << colors[node] << "];" << endl;
        ss.str("");
    }
    while (getline(dot1, line))
    {
        dot2 << line << endl;
    } 
    dot1.close();
    dot2.close();
    
    system("dot graphc.dot -Tpng -o image.png");
}