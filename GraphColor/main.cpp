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

    lprec  *lp;
    lp = make_lp(0, Ncol*Ncol);
    set_minim(lp);
    
    set_add_rowmode(lp, TRUE);
    
    {   // x_v1+x_v2+...+x_vK = 1, v in V
        double row[Ncol*Ncol];
        for ( int i = 0; i < Ncol*Ncol; i++ )
        {    
            row[i] = 1;
        }
        
        add_constraint(lp,row,EQ,1);
    }
    
    {   // y >= k*x_vk, v in V, k=1...K
        double row[Ncol*Ncol];
        for ( int i = 0; i < Ncol*Ncol; i++ )
        {    
            row[i] = 1;
        }
        add_constraint(lp,row,LE,Ncol+1);
    }
    
    {   // 1 >= x_u,k + x_v,k u,v in E, k=1...K
        //double row[] = {0,10,12,8,18};
        //add_constraint(lp,row,LE,1);
    } 
    
    set_add_rowmode(lp, FALSE);
    
    {   // min: y
        double row[] = {1};
        set_obj_fn(lp,row);

        double column[1];
        column[0]= Ncol+1;        
        add_column(lp, column);        
    } 
    
    write_LP(lp, stdout);

    delete_lp(lp);

    GraphOutput(g);
    GraphDisplay();    
    
    return 0;
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
    vector< vector<int> > graph(n,vector<int>(n,0));
    random_device rd;
    mt19937 rnd(rd());
    uniform_int_distribution<int> dist(0,n-1);
    
    for ( vector<int>::size_type i = 0; i < n;  i++ )
    {
        int pick = dist(rnd);        
        graph[i][pick] = 1;
        graph[pick][i] = 1; 
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
    }
    
    
    for ( vector<int>::size_type i = 0; i < graph.size();  i++ )
    {
        for ( vector<int>::size_type j = i+1; j < graph[i].size();  j++ )
        {            
            if ( graph[i][j] != 0 )
            {              
                boost::add_edge(vertices[i],vertices[j],g);
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