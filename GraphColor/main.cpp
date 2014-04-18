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
#include <glpk.h>
#include <boost/graph/undirected_graph.hpp>
#include <boost/graph/graphviz.hpp>

/************************************
 * Namespaces 
 ************************************/
using namespace std;

/************************************
 * Local Types 
 ************************************/

/************************************
 * Local Variables 
 ************************************/


/************************************
 * Local Functions 
 ************************************/
vector< vector<int> > GenerateGraph ( int );
void GraphDisplay ( vector< vector<int> > );
void using_glpk ( vector< vector<int> > );
void using_lpsolve ( vector< vector<int> > );
void GraphOutput ( vector< vector<int> > );

/*******************************************************************************
* Function     : 
* Description  : 
* Arguments    : 
* Returns      : 
* Remarks      : 
********************************************************************************/
int main ( int argc, char* argv[] ) 
{
    int     Ncol;
    
    if ( argc < 2 )
    {
        cout << "enter number of nodes" << endl;
        exit(0);
    }
    istringstream ss(argv[1]);
    if (!(ss >> Ncol))
    {
        cerr << "Invalid number " << argv[1] << '\n';    
        exit(0);
    }
    
    vector< vector<int> > g = GenerateGraph(Ncol);    
    
    cout << "-----------------------glpk" << endl;
    using_glpk(g);
    cout << "-----------------------lpsolve" << endl;
    using_lpsolve(g);

    GraphOutput(g);       
}
/*******************************************************************************
* Function     : 
* Description  : 
* Arguments    : 
* Returns      : 
* Remarks      : 
********************************************************************************/
void using_glpk ( vector< vector<int> > g )
{
    glp_prob* prob = glp_create_prob();
    glp_set_obj_dir(prob, GLP_MIN);     
    int num_vertices = g.size();
    int max_colors = g.size();
//    for (int i = 0; i < num_vertices; ++i)
//        max_colors = std::max(int(g[i].size()) + 1, max_colors);
    
    int y = glp_add_cols(prob, 1);
    glp_set_col_bnds(prob, y, GLP_DB, 1, max_colors); // DB = Double Bound
    glp_set_obj_coef(prob, y, 1.);
    glp_set_col_kind(prob, y, GLP_IV); // IV = Integer Variable
    
    std::vector<std::vector<int> > x(num_vertices, std::vector<int>(max_colors));
    for (int v = 0; v < num_vertices; ++v) 
    {
        for (int k = 0; k < max_colors; ++k)
        {
            x[v][k] = glp_add_cols(prob, 1);
            glp_set_col_kind(prob, x[v][k], GLP_BV); // BV = Binary Variable
        }    
    }
    
    std::vector<int> rows(1, 0);
    std::vector<int> cols(1, 0);
    std::vector<double> coeffs(1, 0.);
    // One vertex must have exactly one color:
    // for each vertex v, sum(x(v, k)) == 1
    for (int v = 0; v < num_vertices; ++v)
    {
        int row_idx = glp_add_rows(prob, 1);
        glp_set_row_bnds(prob, row_idx, GLP_FX, 1, 1); // FX: FiXed bound
        for (int k = 0; k < max_colors; ++k)
        {
            rows.push_back(row_idx); 
            coeffs.push_back(1);
            cols.push_back(x[v][k]);
        }
    }      

    // We ensure we use y colors max:
    // for each vertex v and for each color c,                
    //    y >= (k + 1) * x(v, k)
    for (int v = 0; v < num_vertices; ++v)                    
    {
        for (int k = 0; k < max_colors; ++k)
        {
            int row_idx = glp_add_rows(prob, 1);
            glp_set_row_bnds(prob, row_idx, GLP_LO, 0, -1); // LO = LOwer bound

            rows.push_back(row_idx);
            coeffs.push_back(1);
            cols.push_back(y);                      

            rows.push_back(row_idx);
            coeffs.push_back(- k - 1);
            cols.push_back(x[v][k]);
        }
    }  
    
    // Adjacent vertices cannot have the same color:        
    // for each edge (src, dst) and for each color k,                         
    //    x(src, k) + x(dst, k) <= 1     
    
    for (int src = 0; src < num_vertices; ++src)
    {
        const std::vector<int>& succs = g[src];
        for (int s = 0; s < succs.size(); ++s)
        {
            int dst = succs[s];
            // Ensure we don't add both (u, v) and (v, u)                                    
            if (src > dst)
            {
                for (int k = 0; k < max_colors; ++k)
                {
                    int row_idx = glp_add_rows(prob, 1);
                    glp_set_row_bnds(prob, row_idx, GLP_UP, -1, 1); // UP = UPper bound

                    rows.push_back(row_idx);
                    coeffs.push_back(1);
                    cols.push_back(x[src][k]);

                    rows.push_back(row_idx);
                    coeffs.push_back(1);
                    cols.push_back(x[dst][k]);
                    
                    //cout << x[src][k] << " " << x[dst][k] << endl;
                    //getchar();
                }
            }
        }
    }
    
    glp_load_matrix(prob, rows.size() - 1, &rows[0], &cols[0], &coeffs[0]);
    glp_iocp parm;
    glp_init_iocp(&parm);
    parm.presolve = GLP_ON;
    glp_intopt(prob, &parm);
        
    double solution = glp_mip_obj_val(prob);
    std::cout << "Colors: " << solution << std::endl;
    for (int i = 0; i < num_vertices; ++i)
    {
        std::cout << i << ": ";
        for (int j = 0; j < max_colors; ++j)
            std::cout << glp_mip_col_val(prob, x[i][j]) << " ";
        std::cout << std::endl;
    }
        
}
/*******************************************************************************
* Function     : 
* Description  : 
* Arguments    : 
* Returns      : 
* Remarks      : 
********************************************************************************/
void using_lpsolve ( vector< vector<int> > g )
{  
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
        add_columnex(lp, idx, col, rowno);
    }

    set_add_rowmode(lp, TRUE);      
    
    {   // x_v1+x_v2+...+x_vK = 1, v in V
        int colno[Ncol*Mcol*Mcol];
        double row[Ncol*Mcol*Mcol];
        int idx = 0;     
        for ( int i = 0; i < Ncol; i++ )
        {   
            colno[idx] = idx+1;
            row[idx] = 0;
            idx++;               
            for ( int j = 0; j < Mcol; j++ )
            {
                colno[idx] = idx+1;
                row[idx] = 1;
                idx++;
            }
        }    
        
        add_constraintex(lp, idx, row, colno, EQ, 1);
    } 
    
   {   // y >= k*x_vk, v in V, k=1...K
        int colno[Ncol*Mcol*Mcol];
        double row[Ncol*Mcol*Mcol];
        int idx = 0;  
        for ( int i = 0; i < Ncol; i++ )
        {    
            colno[idx] = idx+1;
            row[idx] = -1;
            idx++;              
            for ( int j = 0; j < Mcol; j++ )
            {
                colno[idx] = idx+1;
                row[idx] = -j-1;
                idx++;
            }
        }
        add_constraintex(lp, idx, row, colno, LE, 0);          
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
                        int colno[Ncol*Mcol*Mcol];
                        double row[Ncol*Mcol*Mcol];
                        int idx = 0;
                                
                        for ( int i = 0; i < Ncol; i++ )
                        {   
                            colno[idx] = idx+1;
                            row[idx] = 0;
                            idx++;                              
                            for (int j = 0; j < Mcol; j++)
                            {
                               colno[idx] = idx+1;
                               row[idx] = 0;
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
    
    write_LP(lp, stdout);
    
    //set_verbose(lp, IMPORTANT);

    set_minim(lp);
    
    if ( solve(lp) == OPTIMAL )
    {
        cout << "objective value " << get_objective(lp) << endl;
        
        {
            double row[Ncol*Mcol];
            get_variables(lp, row);
            for(int j = 0; j < Ncol; j++)
              printf("%s: %f\n", get_col_name(lp, j + 1), row[j]);
        }
        
    }
    delete_lp(lp);
 
    //GraphDisplay(g);    
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
//    random_device rd;
//    mt19937 rnd(rd());
//    uniform_int_distribution<int> dist(0,n-1);
//    
//    for ( vector<int>::size_type i = 0; i < n; i++ )
//    {
//        vector<int> tmp;
//        tmp.push_back(i);
//        for ( vector<int>::size_type j = 0; j < n; j++ )
//        {
//            if ( j == i )
//                continue;
//            
//            int pick = dist(rnd);
//            for ( vector<int>::size_type k = 0; k < tmp.size(); k++ )
//            {
//                if ( pick == tmp[k] )
//                {
//                    break;
//                }
//                else
//                {
//                    tmp.push_back(pick);
//                }
//            }
//        }
//        
//        sort( tmp.begin()+1, tmp.end() );
//        for ( vector<int>::size_type k = 1; k < tmp.size()-1; k++ )
//        {
//            if ( tmp[k]==tmp[k+1] )
//            {
//                tmp.erase(tmp.begin()+k+1);
//                k--;
//            }
//        }
//        
//        graph.push_back(tmp);
//    }
    
    
    vector<int> tmp;
    tmp.push_back(0);
    tmp.push_back(1);
    tmp.push_back(2);
    tmp.push_back(4);
    tmp.push_back(5);
    graph.push_back(tmp);
    vector<int> tmp1;
    tmp1.push_back(1);
    tmp1.push_back(0);
    tmp1.push_back(2);
    tmp1.push_back(3);
    tmp1.push_back(5);    
    graph.push_back(tmp1);
    vector<int> tmp2;
    tmp2.push_back(2);
    tmp2.push_back(0);
    tmp2.push_back(1);
    tmp2.push_back(3);
    tmp2.push_back(4);    
    graph.push_back(tmp2);
    vector<int> tmp3;
    tmp3.push_back(3);
    tmp3.push_back(1);
    tmp3.push_back(2);
    tmp3.push_back(4);
    tmp3.push_back(5);    
    graph.push_back(tmp3);
    vector<int> tmp4;
    tmp4.push_back(4);
    tmp4.push_back(0);
    tmp4.push_back(2);
    tmp4.push_back(3);
    tmp4.push_back(5);    
    graph.push_back(tmp4);
    vector<int> tmp5;
    tmp5.push_back(5);
    tmp5.push_back(0);
    tmp5.push_back(1);
    tmp5.push_back(3);
    tmp5.push_back(4);    
    graph.push_back(tmp5); 
    
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
void GraphOutput ( vector< vector<int> > graph )
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
    system("dot graph.dot -Tpng -o image.png");
}