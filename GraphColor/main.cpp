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
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <lpsolve/lp_lib.h>
#include <random>
#include <vector>
#include <boost/graph/undirected_graph.hpp>
#include <boost/graph/graphviz.hpp>
#include <string>
#include <fstream>

#include <glpk.h>

         
/************************************
 * Namespaces 
 ************************************/
using namespace std;
using namespace cv;

/************************************
 * Local Types 
 ************************************/
typedef boost::undirected_graph<boost::no_property> Graph;

/************************************
 * Local Variables 
 ************************************/
bool display;


/************************************
 * Local Functions 
 ************************************/
vector< vector<int> > GenerateGraph ( int );
void GraphOutput ( vector< vector<int> > );
void GraphDisplay ( void );

/*******************************************************************************
* Function     : 
* Description  : 
* Arguments    : 
* Returns      : 
* Remarks      : 
********************************************************************************/
//int main ( int argc, char* argv[] ) 
//{  
//    int     Ncol;
//    int     Mcol;
//    
//    if ( argc < 2 )
//    {
//        cout << "enter number of nodes" << endl;
//        exit(0);
//    }
//    istringstream ss(argv[1]);
//    if (!(ss >> Ncol))
//    {
//        cerr << "Invalid number " << argv[1] << '\n';    
//        exit(0);
//    }
//    
//    vector< vector<int> > g = GenerateGraph(Ncol);
//    
//    for (int i = 0; i < Ncol; ++i)
//        Mcol = max(int(g[i].size()) + 1, Mcol); 
//    
//    cout << Mcol << endl;
//
//    lprec  *lp;
//    lp = make_lp(0, Ncol*Mcol);
//    set_minim(lp);
//    
//    set_add_rowmode(lp, TRUE);
//    
//    {   // x_v1+x_v2+...+x_vK = 1, v in V
//        double row[Ncol*Ncol];
//        for ( int i = 0; i < Ncol; i++ )
//        {    
//            for ( int j = 0; j < Mcol; j++ )
//            {
//                row[i*Ncol+j] = 1;
//            }
//        }
//        
//        add_constraint(lp,row,EQ,1);
//    } 
//    
//   {   // y >= k*x_vk, v in V, k=1...K
//        double row[Ncol*Mcol];
//        for ( int i = 0; i < Ncol; i++ )
//        {    
//            for ( int j = 0; j < Mcol; j++ )
//            {
//                row[i*Ncol+j] = -j-1;
//            }
//        }
//        add_constraint(lp,row,LE,Mcol);  
//    }    
//    
//    {   // 1 >= x_u,k + x_v,k u,v in E, k=1...K
//        double row[Ncol*Mcol];
//        for ( int i = 0; i < Ncol; i++ )
//        {    
//            const std::vector<int>& succs = g[i];
//            for ( int j = 0; j < succs.size(); j++ )
//            {
//                int dst = succs[j];
//                // Ensure we don't add both (u, v) and (v, u)                                    
//                if (i > dst)
//                {
//                    for (int k = 0; k < Mcol; ++k)  
//                    {
//                        row[i*Ncol+j] = g[dst][k];
//                    }
//                }
//            }
//        }
//        add_constraint(lp,row,LE,1);
//    } 
//    
//    set_add_rowmode(lp, FALSE);
//    
//    {   // min: y
//        double row[Ncol+1];
//        row[0]= Mcol;         
//        set_obj_fn(lp,row);   
//        
//        double column[Ncol+1];
//        column[0]= Mcol;    
//        add_column(lp, column);
//    } 
//       
//    write_LP(lp, stdout);
//    
//    set_verbose(lp, IMPORTANT);
//
//    if ( solve(lp) == OPTIMAL )
//    {
//        cout << "objective value " << get_objective(lp) << endl;
//        
//        {
//            double row[Ncol];
//            get_variables(lp, row);
//            for ( int j = 0; j < Ncol; j++ )
//                cout << row[j] << endl;
//        }
//        
//    }
//    delete_lp(lp);
//
//    GraphOutput(g);
//    GraphDisplay();    
//     
//    return 0;
//}
int main ( int argc, char* argv[] ) 
{
    int     Ncol;
    int     Mcol;
    
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
    
    glp_prob* prob = glp_create_prob();
    glp_set_obj_dir(prob, GLP_MIN);     
    int num_vertices = g.size();
    int max_colors = num_vertices;
    for (int i = 0; i < num_vertices; ++i)
        max_colors = std::max(int(g[i].size()) + 1, max_colors);
    
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
    
    for (int i = 0; i < g.size(); ++i)
    {
        for (int j = 0; j < g[i].size(); ++j)
        {
            cout << g[i][j] << " ";
        }
        cout << endl;
    }
    GraphOutput(g);
    GraphDisplay();     
    
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
        for ( vector<int>::size_type j = 0; j < n; j++ )
        {
            if ( j == i )
                continue;
            
            int pick = dist(rnd);
            for ( vector<int>::size_type k = 0; k < tmp.size(); k++ )
            {
                if ( pick == tmp[k] )
                {
                    break;
                }
                else
                {
                    tmp.push_back(pick);
                }
            }
        }
        
        sort( tmp.begin()+1, tmp.end() );
        for ( vector<int>::size_type k = 1; k < tmp.size()-1; k++ )
        {
            if ( tmp[k]==tmp[k+1] )
            {
                tmp.erase(tmp.begin()+k+1);
                k--;
            }
        }
        
        graph.push_back(tmp);
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
void GraphOutput ( vector< vector<int> > graph )
{
    ofstream fout("graph.dot");

    Graph g;
        
    vector<Graph::vertex_descriptor> vertices;
    for ( vector<int>::size_type i = 0; i < graph.size();  i++ )
    {
        Graph::vertex_descriptor v = g.add_vertex();
        vertices.push_back(v);
        
        g[i][j]
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
/*******************************************************************************
* Function     : 
* Description  : 
* Arguments    : 
* Returns      : 
* Remarks      : 
********************************************************************************/
void GraphDisplay ( void )
{
    Mat image = imread("image.png", CV_LOAD_IMAGE_COLOR);
    namedWindow( "graph", WINDOW_AUTOSIZE );
    imshow( "graph", image );
    waitKey();
    destroyAllWindows();
}