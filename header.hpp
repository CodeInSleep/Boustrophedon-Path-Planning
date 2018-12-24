#ifndef _HEADER_H
#define _HEADER_H

#include <boost/geometry/geometry.hpp>

using namespace boost::geometry;
using namespace std;

typedef model::d2::point_xy<double> point_2d;
typedef model::segment<point_2d> segment_2d;
typedef model::polygon<point_2d, false> polygon_2d;
typedef vector<point_2d> pt_seq_t;
typedef vector<segment_2d> edge_seq_t;
typedef pt_seq_t::iterator pt_iter_t;
typedef pt_seq_t::const_iterator pt_const_iter_t;
typedef edge_seq_t::iterator edge_iter_t;
typedef edge_seq_t::const_iterator edge_const_iter_t;
typedef pair<point_2d, segment_2d> intersect_t;
typedef vector<intersect_t> intersect_seq_t;
typedef intersect_seq_t::iterator intersect_iter_t;

#endif