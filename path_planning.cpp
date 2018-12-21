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
#include <boost/geometry/algorithms/intersection.hpp>
#include "Event.hpp"
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
void BoustrophedonCell::update(Event<point_2d> event) {
    //TODO: check if the current cell is involved in the event (i.e. in contact with the event line)
}

void BoustrophedonCell::update(Event<segment_2d> event) {
    //TODO: check if the current cell is involved in the event (i.e. in contact with the event line)
}



std::vector<point_2d> gen_path(BoustrophedonCell cell) {
    std::vector<point_2d> generatedPath;

    return generatedPath;
}

template <int D>
pair<pt_iter_t, pt_iter_t> minmax(pt_seq_t &seq) {
    pt_seq_t::iterator min_ele_iter = seq.begin();
    pt_seq_t::iterator max_ele_iter = seq.begin();

    for (pt_seq_t::iterator iter = seq.begin(); iter < seq.end(); ++iter) {
        double xCoord = get<D>(*iter);
        if (xCoord < get<D>(*min_ele_iter))
            min_ele_iter = iter;
        if (xCoord > get<D>(*max_ele_iter))
            max_ele_iter = iter;
    }

    return pair<pt_iter_t, pt_iter_t>(min_ele_iter, max_ele_iter);
}

pt_seq_t leftCircularShift(pt_seq_t &poly) {

    if (poly.size() < 3)
        throw invalid_argument("received a polygon with less than 3 vertices..");
    // left circular shift so that the leftmost point of
    // the polygon is in the first position, followed by
    // vertices in counterclockwise traversal. Assume poly
    // is in counterclockwise traversal.
    pair<pt_iter_t, pt_iter_t> result = minmax<0>(poly);
    pt_iter_t leftShift = result.first;
    pt_iter_t leftShiftCopy = leftShift;
    if (leftShift == poly.begin()) {
        leftShift = (get<0>(*leftShift) == get<0>(*(poly.end()-1))) ? poly.end()-1 : leftShift;
    }
    else if (get<0>(*leftShift) == (get<0>(*(leftShiftCopy--)))) {
        leftShift = leftShiftCopy;
    }
    cout << "left shift from: " << get<0>(*leftShift) << ", " << get<1>(*leftShift) << endl;
    pt_seq_t shifted_poly;
    
    pt_iter_t vtx_iter;
    for (vtx_iter = leftShift; vtx_iter < poly.end(); vtx_iter++) {
        shifted_poly.push_back(*vtx_iter);
    }
    for (vtx_iter = poly.begin(); vtx_iter < leftShift; vtx_iter++) {
        shifted_poly.push_back(*vtx_iter);
    }

    return shifted_poly;
}

template <typename T>
void printSeq(const vector<T> &seq) {
    // loop through coordinates
    for (typename vector<T>::const_iterator iter = seq.begin(); iter != seq.end(); ++iter) {
        cout << dsv(*iter);
    }

    cout << "\n";
}

template<>
void printSeq<intersect_t>(const vector<intersect_t> &seq) {
    for (vector<intersect_t>::const_iterator iter = seq.begin(); iter != seq.end(); ++iter) {
        cout << dsv((*iter).first) << ", " << dsv((*iter).second) << endl;
    }
}

Polygon::Polygon(double coords[][2], int num_of_pts) {
    // left circular shift so that the leftmost vertex is in the first position
    for (size_t i = 0; i < num_of_pts; ++i) {
        vertices.push_back(make<point_2d>(coords[i][0], coords[i][1]));
    }

    vertices = leftCircularShift(vertices);

    // ending point is the starting point
    vertices.push_back(vertices[0]);

    // create edges
    for (pt_const_iter_t iter = vertices.begin(); iter != vertices.end()-1; ++iter) {
        pt_const_iter_t nextPtIter = iter+1;
        edges.push_back(segment_2d(*iter, *nextPtIter));
    }

    cout << "vertices: " << endl;
    printSeq(vertices);
    cout << "edges: " << endl;
    printSeq(edges);

    // make a boost polygon representation from vertices
    assign_points<polygon_2d, vector<point_2d> >(poly, vertices);
}

intersect_seq_t Polygon::line_intersect(double xCoord) {
    pair<pt_iter_t, pt_iter_t> result = minmax<1>(vertices);

    // create vertical line that span the height of the polygon
    segment_2d scanLine = segment_2d(
        point_2d(xCoord, get<1>(*result.first)),
        point_2d(xCoord, get<1>(*result.second)));
    // cout << "scanLine: " << dsv(scanLine) << endl;

    intersect_seq_t intersects;
    for (edge_const_iter_t e_iter = edges.begin();e_iter != edges.end(); ++e_iter) {
        vector<point_2d> output;
        intersection(*e_iter, scanLine, output);

        // cout << "output: " << endl;
        // printSeq(output);

        if (output.size() != 0)
            intersects.push_back(intersect_t(output[0], *e_iter));
    }

    return intersects;
}

SurveyArea::SurveyArea(double coords[][2], int num_of_pts)
    : Polygon(coords, num_of_pts) {}

int main(void)
{
    double coords[4][2] = { {0, 400}, {400, 400}, {400, 0}, {0, 0} };
    int num_of_pts = sizeof(coords)/sizeof(coords[0]);
    double (*coords_ptr)[2] = coords;
    SurveyArea sa(coords_ptr, num_of_pts);

    intersect_seq_t intersects = sa.line_intersect(200);
    cout << "intersects: " << endl;
    printSeq(intersects);

    return 0;
}