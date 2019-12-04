// standard includes
#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
// includes for defining the Voronoi diagram adaptor
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Delaunay_triangulation_2.h>
//#include <CGAL/Voronoi_diagram_2.h>
//#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
//#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>


#include <CGAL/Segment_Delaunay_graph_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_2.h>
#include <CGAL/Segment_Delaunay_graph_adaptation_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_adaptation_policies_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Polygon_2.h>

#include <unordered_map>
#include <unordered_set>


typedef CGAL::Exact_predicates_inexact_constructions_kernel                  K;
typedef CGAL::Segment_Delaunay_graph_traits_2<K> Gt;
typedef CGAL::Segment_Delaunay_graph_2<Gt> SDG2;

typedef CGAL::Segment_Delaunay_graph_adaptation_traits_2<SDG2> AT;

typedef CGAL::Segment_Delaunay_graph_degeneracy_removal_policy_2<SDG2> AP;

typedef CGAL::Voronoi_diagram_2<SDG2, AT, AP> VD;

//typedef CGAL::Delaunay_triangulation_2<K>                                    DT;
//typedef CGAL::Delaunay_triangulation_adaptation_traits_2<DT>                 AT;
//typedef CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<DT> AP;
//typedef CGAL::Voronoi_diagram_2<DT,AT,AP>                                    VD;
// typedef for the result type of the point location
typedef AT::Site_2                    Site_2;
typedef AT::Point_2                   Point_2;
typedef VD::Locate_result             Locate_result;
typedef VD::Vertex_handle             Vertex_handle;
typedef VD::Face_handle               Face_handle;
typedef VD::Halfedge_handle           Halfedge_handle;
typedef VD::Ccb_halfedge_circulator   Ccb_halfedge_circulator;
typedef VD::Bounded_halfedges_iterator BHE_Iter;
typedef VD::Halfedge Halfedge;
typedef VD::Vertex Vertex;
typedef CGAL::Polygon_2<K> Polygon_2;


#define EPSILON 0.01
/*
void print_endpoint(Halfedge_handle e, bool is_src) {
  std::cout << "\t";
  if ( is_src ) {
    if ( e->has_source() )  std::cout << e->source()->point() << std::endl;
    else  std::cout << "point at infinity" << std::endl;
  } else {
    if ( e->has_target() )  std::cout << e->target()->point() << std::endl;
    else  std::cout << "point at infinity" << std::endl;
  }
}
*/
double radius_function(std::vector<Point_2> &points, Point_2 v);
double distance_to_segment(Point_2 a, Point_2 b, Point_2 x);

int main()
{

  CGAL::set_pretty_mode(std::cout);

  VD vd;
  Site_2 site;
  std::cout << "Using the vd interface" << std::endl;
  
  std::vector<Point_2> points = {
    Point_2(3.0, 3.0),
    Point_2(3.0, -3.0),
    Point_2(-3.0, -3.0),
    Point_2(-3.0, 3.0),
    Point_2(3.0, 3.0)
  };

  double eps = 0.1;


  
  // define sites
  
  for (std::size_t i = 0; i<points.size()-1; i++) {
    std::cout << points[i] << std::endl;
    site = Site_2::construct_site_2(points[i], points[i+1]);
    vd.insert(site);
  }

  assert( vd.is_valid() );

  std::cout << std::endl;

  std::cout << "Iterating over edges" << std::endl;
  std::cout << "number of halfedges: " << vd.number_of_halfedges() << std::endl;
  //std::cout << "example unbounded_halfedge: " << vd.bounded_halfedge() << std::endl;
  
  // associate vertex 

  


  // Assign each vertex an id

  std::unordered_map<int, std::vector<int>> connect;
  std::unordered_map<int, Point_2> id_vtx;

  BHE_Iter edge_iter = vd.bounded_halfedges_begin();
  int id_count = 0;
  bool found;
  for (;edge_iter != vd.bounded_halfedges_end(); edge_iter++) {
    Halfedge halfedge = *edge_iter;
    Vertex v1 = *(halfedge.source());
    int v1_id;
    found = false;
    for (auto j = id_vtx.begin(); j != id_vtx.end(); j++) {
      if (CGAL::squared_distance(v1.point(), j->second) < EPSILON) {
        v1_id = j->first;
        found = true;
      }
    }
    if (!found) {
      id_vtx[id_count] = v1.point();
      v1_id = id_count;
      id_count++;
    }


    Vertex v2 = *(halfedge.target());
    int v2_id;
    found = false;
    for (auto j = id_vtx.begin(); j != id_vtx.end(); j++) {
      if (CGAL::squared_distance(v2.point(), j->second) < EPSILON) {
        v2_id = j->first;
        found = true;
      }
    }
    if (!found) {
      id_vtx[id_count] = v2.point();
      v2_id = id_count;
      id_count++;
    }

    connect[v1_id].push_back(v2_id);

    //std::cout << halfedge << std::endl;
  }


  
  std::cout << "id points" << std::endl;
  for (auto i = id_vtx.begin(); i != id_vtx.end(); i++) {
    std::cout << i->first << ": " << i->second << std::endl;
  }

  std::cout << "connections" << std::endl;
  for (auto i = connect.begin(); i != connect.end(); i++) {
    std::cout << i->first << ": ";
    for (auto j = i->second.begin(); j != i->second.end(); j++) {
      std::cout << *j << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl << std::endl;
  std::cout << "Making Square" << std::endl;

  std::ofstream myfile;
  myfile.open ("square.txt");

  int id_p1, id_p2;
  for ( auto it = connect.begin(); it != connect.end(); ++it) {
    id_p1 = it->first;
    for (auto k = (it->second).begin(); k != (it->second).end(); k++) {
      id_p2 = *k;
      if (id_p1 < id_p2) {
        continue;
      }
      Point_2 p1 = id_vtx[id_p1];
      Point_2 p2 = id_vtx[id_p2];
      std::cout << p1 << " " << p2 << std::endl;
      double dist = sqrt(CGAL::squared_distance(p1, p2));
      std::cout << dist << std::endl;
      int num_pts = (int) (dist/EPSILON);
      double frac;
      for (int t = 1; t < num_pts; t++) {\
        frac = t*EPSILON/dist;
        Point_2 q = Point_2(p1.x()*frac + p2.x()*(1-frac), p1.y()*frac + p2.y()*(1-frac));
        myfile << q.x() << " " << q.y() << " " << radius_function(points, q) << std::endl;
      }
      myfile << p1.x() << " " << p1.y() << " " << radius_function(points, p1) << std::endl;
      myfile << p2.x() << " " << p2.y() << " " << radius_function(points, p2) << std::endl;
    }
  }
  myfile.close();
  /*
  BHE_Iter edge_iter = vd.bounded_halfedges_begin();
  for (;edge_iter != vd.bounded_halfedges_end(); edge_iter++) {
    Halfedge halfedge = *edge_iter;
    Vertex v1 = *(halfedge.source());
    Vertex v2 = *(halfedge.target());
    Point_2 p1 = v1.point();
    Point_2 p2 = v2.point();
    double dist = sqrt(CGAL::squared_distance(p1, p2));
    int num_pts = dist/eps;

    for (int k = 0; k < num_pts; k++)

    //std::cout << halfedge << std::endl;
  }
  */


}

double radius_function(std::vector<Point_2> &points, Point_2 v) {
  double min_sq_dist = -1, sq_dist;
  Point_2 p1, p2;
  for (std::size_t i = 0; i<points.size()-1; i++) {
    p1 = points[i];
    p2 = points[i+1];
    sq_dist = distance_to_segment(p1, p2, v);
    if (min_sq_dist < 0 || sq_dist < min_sq_dist) {
      min_sq_dist = sq_dist;
    }
  }
  return sqrt(min_sq_dist);
}

// return smallest distance to line segment
double distance_to_segment(Point_2 a, Point_2 b, Point_2 p) {
  double vx, vy, ux, uy, t, pa, pb;

  vx = b.x() - a.x();
  vy = b.y() - a.y();
  ux = a.x() - p.x();
  uy = a.y() - p.y();

  t = -(vx*ux + vy*uy)/(vx*vx + vy*vy);
  if (0 <= t && t <= 1) {
    Point_2 min = Point_2((1-t)*a.x() + t*b.x(), (1-t)*a.y() + t*b.y());
    return CGAL::squared_distance(p, min);
  }
  pa = CGAL::squared_distance(a, p);
  pb = CGAL::squared_distance(b, p);
  if (pa < pb) {
    return pa;
  }
  else {
    return pb;
  }
}