#include <algorithm> // for reverse, unique
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>

#include <boost/geometry/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>
#include <boost/geometry/geometries/adapted/c_array.hpp>
#include <boost/geometry/algorithms/assign.hpp>
#include "path_planning.hpp"

using namespace std;
std::string boolstr(bool v)
{
    return v ? "true" : "false";
}

struct drone_test_configuration {
    double altitude;
};

BoustrophedonCell::BoustrophedonCell(double leftEndPoint, segment_2d *startCeilSegmentPtr, segment_2d *startFloorSegmentPtr)
    : leftEndPoint(leftEndPoint) {
        ceilSegmentPtrs.push_back(startCeilSegmentPtr);
        floorSegmentPtrs.push_back(startFloorSegmentPtr);
}

// event types:
// BEGIN: 0
// END: 1
// IN: 2
// OUT: 3
// MIDDLE: 4
void BoustrophedonCell::update(Event event) {
    //TODO: check if the current cell is involved in the event (i.e. in contact with the event line)
}


std::vector<point_2d> gen_path(BoustrophedonCell cell) {
    std::vector<point_2d> generatedPath;

    return generatedPath;
}

pair<coord_seq_t::iterator, coord_seq_t::iterator> minmax(coord_seq_t &seq, int dim = 0) {
    coord_seq_t::iterator min_ele_iter = seq.begin();
    coord_seq_t::iterator max_ele_iter = seq.begin();

    for (coord_seq_t::iterator iter = seq.begin(); iter < seq.end(); ++iter) {
        double xCoord = get<0>(*iter);
        if (xCoord < get<0>(*min_ele_iter))
            min_ele_iter = iter;
        if (xCoord > get<0>(*max_ele_iter))
            max_ele_iter = iter;
    }

    return pair<coord_seq_t::iterator, coord_seq_t::iterator>(min_ele_iter, max_ele_iter);
}

coord_seq_t leftCircularShift(coord_seq_t &poly) {

    if (poly.size() < 3)
        throw invalid_argument("received a polygon with less than 3 vertices..");
    // left circular shift so that the leftmost point of
    // the polygon is in the first position, followed by
    // vertices in counterclockwise traversal. Assume poly
    // is in counterclockwise traversal.
    auto result = minmax(poly);
    vec_iter leftShift = result.first;
    vec_iter leftShiftCopy = leftShift;
    if (leftShift == poly.begin()) {
        leftShift = (get<0>(*leftShift) == get<0>(*(poly.end()-1))) ? poly.end()-1 : leftShift;
    }
    else if (get<0>(*leftShift) == (get<0>(*(leftShiftCopy--)))) {
        leftShift = leftShiftCopy;
    }
    cout << "left shift from: " << get<0>(*leftShift) << ", " << get<1>(*leftShift) << endl;
    coord_seq_t shifted_poly;
    
    vec_iter vtx_iter;
    for (vtx_iter = leftShift; vtx_iter < poly.end(); vtx_iter++) {
        shifted_poly.push_back(*vtx_iter);
    }
    for (vtx_iter = poly.begin(); vtx_iter < leftShift; vtx_iter++) {
        shifted_poly.push_back(*vtx_iter);
    }

    return shifted_poly;
}

void printVec(const vector<point_2d> &seq) {
    // loop through coordinates
    for (const point_2d &ele : seq) {
        cout << get<0>(ele) << ", " << get<1>(ele) << endl;
    }
}

SurveyArea::SurveyArea(double coords[][2], int num_of_pts) {
    // left circular shift so that the leftmost vertex is in the first position
    for (size_t i = 0; i < num_of_pts; ++i) {
        vertices.push_back(make<point_2d>(coords[i][0], coords[i][1]));
    }

    cout << "received: " << vertices.size() << " vertices" << endl;
    vertices = leftCircularShift(vertices);

    // ending point is the starting point
    vertices.push_back(make<point_2d>(coords[0][0], coords[0][1]));
    

    // make a boost polygon representation from vertices
    assign_points<polygon_2d, vector<point_2d> >(poly, vertices);
}


int main(void)
{
    double coords[4][2] = { {0, 400}, {400, 400}, {400, 0}, {0, 0} };
    int num_of_pts = sizeof(coords)/sizeof(coords[0]);
    double (*coords_ptr)[2] = coords;
    SurveyArea sa(coords_ptr, num_of_pts); 
    return 0;
}