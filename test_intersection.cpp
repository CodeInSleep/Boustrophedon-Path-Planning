#include <iostream>
#include <string>
#include <vector>

#include <boost/geometry/geometry.hpp>
#include <boost/geometry/geometries/segment.hpp>
#include <boost/geometry/algorithms/intersection.hpp>

using namespace boost::geometry;
using namespace std;

int main() {
	typedef model::d2::point_xy<double> point_2d;
	typedef model::linestring<point_2d> linestring_2d;
    typedef boost::geometry::model::polygon<boost::geometry::model::d2::point_xy<double> > polygon;
	typedef std::deque<point_2d>::iterator deque_iter;

	double y_max = 10;
	double y_min = 0;
	double x = 2;

	polygon green;
	linestring_2d ls;

	boost::geometry::read_wkt(
        "POLYGON((0 0, 0 4, 4 4, 4 0, 0 0))", green);

	correct(green);
	append(ls, make<point_2d>(x, y_min));
	append(ls, make<point_2d>(x, y_max));

	deque<point_2d> output;

	intersection(green, ls, output);
	for (deque_iter it = output.begin(); it != output.end(); ++it) {
		cout << get<0>(*it) << ", " << get<1>(*it) << std::endl;
	}
}
