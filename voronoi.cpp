// standard includes
#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <string>
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


#define EPSILON 0.03
#define THRESHOLD 0.0001
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

int main(int argc, char*argv[])
{

  CGAL::set_pretty_mode(std::cout);


  std::cout << "Using the vd interface" << std::endl;
  std::vector<Point_2> points;
  

  if (argc > 1) {
    std::string line;
    std::ifstream pointfile;
    pointfile.open(argv[1]);
    if (pointfile.is_open()) {
      std::cout << "Reading from file " << argv[1] << std::endl;

      double x, y;
      while (getline(pointfile, line)) {
        std::stringstream ss(line);
        ss >> x;
        ss >> y;
        points.push_back(Point_2(x, y));
      }

      pointfile.close();
      std::cout << "Finished reading " << argv[1] << std::endl << std::endl;
    }
  }
  else {
    std::cout << "No polygon file to process." << std::endl;
    return 0;
  }


  double eps = 0.01;


  VD vd;
  Site_2 site;
  int n = points.size();
  // define sites
  
  std::cout << "Create Voronoi diagram" << std::endl;
  for (std::size_t i = 0; i< n-1; i++) {
    std::cout << points[i] << std::endl;
    site = Site_2::construct_site_2(points[i], points[i+1]);
    vd.insert(site);
  }
  // close the polygon - connect the last and first sites
  std::cout << points[n-1] << std::endl;
  site = Site_2::construct_site_2(points[0], points[n-1]);
  vd.insert(site);

  assert( vd.is_valid() );
  std:: cout << "Voronoi Diagram is valid" << std::endl << std::endl;




  std::cout << "Iterating over edges" << std::endl;
  std::cout << "number of halfedges: " << vd.number_of_halfedges() << std::endl;
  //std::cout << "example unbounded_halfedge: " << vd.bounded_halfedge() << std::endl;

  std::unordered_map<int, std::vector<int>> connect;
  std::unordered_map<int, Point_2> id_vtx;

  BHE_Iter edge_iter = vd.bounded_halfedges_begin();
  int id_count = 1;
  bool found;
  for (;edge_iter != vd.bounded_halfedges_end(); edge_iter++) {
    Halfedge halfedge = *edge_iter;
    Vertex v1 = *(halfedge.source());
    int v1_id;
    found = false;
    for (auto j = id_vtx.begin(); j != id_vtx.end(); j++) {
      if (CGAL::squared_distance(v1.point(), j->second) < THRESHOLD) {
        v1_id = j->first;
        found = true;
        break; // caution about this Anis
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
      if (CGAL::squared_distance(v2.point(), j->second) < THRESHOLD) {
        v2_id = j->first;
        found = true;
        break; // caution abut this Anis
      }
    }
    if (!found) {(the subjects for which you plan to apply for transfer credit) 
      id_vtx[id_count] = v2.point();
      v2_id = id_count;
      id_count++;
    }

    connect[v1_id].push_back(v2_id);

    //std::cout << halfedge << std::endl;
  }

  // Calculate reflex edges
  /*
  std::unordered_set<int> reflex_vertices;
  for (std::size_t i = 0; i< n; i++) {
    v = Point_2(); // really a vector
    w = Point_2(); // really a vector
    std::cout << points[i] << std::endl;
    site = Site_2::construct_site_2(points[(i+n)%n], points[(i+1+n)]);
    vd.insert(site);

  }
  */


  // Create medial axis vertex dump
  std::string medial_axis_nodes_file = std::string(argv[1]) + "_med_axis_nodes";
  std::string medial_axis_edges_file = std::string(argv[1]) + "_med_axis_edges";
  std::ofstream med_axis_nodes;
  std::ofstream med_axis_edges;
  
  
  std::cout << "Writing medial axis nodes points" << std::endl;
  med_axis_nodes.open(medial_axis_nodes_file);
  if (med_axis_nodes.is_open()) {
    for (auto i = id_vtx.begin(); i != id_vtx.end(); i++) {
      std::cout << i->first << ": " << i->second << std::endl;
      med_axis_nodes << i->first << " " << i->second.x() << " " << i->second.y() << std::endl;
    }
    med_axis_nodes.close();
  }
  else {
    std::cout << "Error dumping medial axis points to a file" << std::endl;
  }
  std::cout << std::endl;

  // Create medial axis edge dump
  std::cout << "Writing medial axis edges" << std::endl;
  med_axis_edges.open(medial_axis_edges_file);
  if (med_axis_edges.is_open()) {
    for (auto i = connect.begin(); i != connect.end(); i++) {
      std::cout << i->first << ": ";
      for (auto j = i->second.begin(); j != i->second.end(); j++) {
        std::cout << *j << " ";
				if (i->first < *j) {
					med_axis_edges << i->first << " " << *j << std::endl;
				}
      }
      std::cout << std::endl;
    }
    med_axis_edges.close();
  }

  std::cout << std::endl;
  std::cout << "Calculating segmented medial axis" << std::endl;
  std::string sampled_axis_nodes_file = std::string(argv[1]) + "_sampled_axis_nodes";
  std::string sampled_axis_edges_file = std::string(argv[1]) + "_sampled_axis_edges";
  std::ofstream sampled_axis_nodes;
  std::ofstream sampled_axis_edges;

  sampled_axis_nodes.open (sampled_axis_nodes_file);
  sampled_axis_edges.open(sampled_axis_edges_file);
  if (!sampled_axis_nodes.is_open() || !sampled_axis_edges.is_open()) {
    std::cout << "Error creating sampled medial axis files" << std::endl;
    return 0;
  }

  std::unordered_set<int> visited_nodes;
  int id_p1, id_p2;
  for ( auto it = connect.begin(); it != connect.end(); ++it) {
    id_p1 = it->first;
    for (auto k = (it->second).begin(); k != (it->second).end(); k++) {
      id_p2 = *k;
      if (id_p1 < id_p2) {
        // We are iterating over halfedges so edges appear twice
        continue;
      }
      Point_2 p1 = id_vtx[id_p1];
      Point_2 p2 = id_vtx[id_p2];

      std::cout << p1 << " " << p2 << std::endl;
      double dist = sqrt(CGAL::squared_distance(p1, p2));
      //std::cout << dist << std::endl;

      int num_pts = (int) (dist/EPSILON);
      double frac;
      int prev_point = id_p1;
      std::cout << "prev_point" << prev_point << std::endl;
      for (int t = 1; t < num_pts; t++) {
        frac = t*EPSILON/dist;
        Point_2 q = Point_2(p2.x()*frac + p1.x()*(1-frac), p2.y()*frac + p1.y()*(1-frac));
        sampled_axis_nodes << id_count << " " << q.x() << " " << q.y() << " " << radius_function(points, q) << std::endl;
        sampled_axis_edges << prev_point << " " << id_count << std::endl;
        prev_point = id_count;
        // for the last point, connect to the second edge point
        if (t == num_pts-1) {
          sampled_axis_edges << id_count << " " << id_p2 << std::endl;
        }

        // id_count is the count from earlier part of code that counts medial axis points
        id_count++;
      }
      if (visited_nodes.find(id_p1) == visited_nodes.end()) {
        sampled_axis_nodes << id_p1 << " " << p1.x() << " " << p1.y() << " " << radius_function(points, p1) << std::endl;
				visited_nodes.insert(id_p1);
      }
      if (visited_nodes.find(id_p2) == visited_nodes.end()) {
        sampled_axis_nodes << id_p2 << " " << p2.x() << " " << p2.y() << " " << radius_function(points, p2) << std::endl;
				visited_nodes.insert(id_p2);
      }
    }
  }
  sampled_axis_nodes.close();
  sampled_axis_edges.close();
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
	p1 = points[0];
	p2 = points[points.size()-1];
	sq_dist = distance_to_segment(p1, p2, v);
	if (min_sq_dist < 0 || sq_dist < min_sq_dist) {
		min_sq_dist = sq_dist;
	}
  return sqrt(min_sq_dist);
}

double distance_to_segment(Point_2 a, Point_2 b, Point_2 p) {
  double l2 = CGAL::squared_distance(a, b);
  if (l2 == 0.0) {
    return CGAL::squared_distance(a, p);
  }
  double dot = (p.x() - a.x())*(b.x() - a.x()) + (p.y() - a.y())*(b.y() - a.y());
  double t = std::max(0.0, std::min(1.0, dot / l2));
  Point_2 proj = Point_2(a.x() + t*(b.x() - a.x()), a.y() + t*(b.y() - a.y()));
  return CGAL::squared_distance(p, proj);
}

// return smallest distance to line segment
double distance_to_segment_2(Point_2 a, Point_2 b, Point_2 p) {


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
