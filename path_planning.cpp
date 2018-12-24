#include <algorithm> // for reverse, unique
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>

#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>
#include <boost/geometry/geometries/adapted/c_array.hpp>
#include <boost/geometry/algorithms/assign.hpp>
#include <boost/geometry/algorithms/intersection.hpp>
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

void BoustrophedonCell::update(IEvent *event) {
    //TODO: check if the current cell is involved in the event (i.e. in contact with the event line)
}

std::vector<point_2d> gen_path(BoustrophedonCell cell) {
    std::vector<point_2d> generatedPath;

    return generatedPath;
}

template <int D, bool increasing>
struct compare {
    bool operator()(const pt_iter_t &a_iter, const pt_iter_t &b_iter) {
        return (increasing) ? get<D>(*a_iter) < get<D>(*b_iter)
                        : get<D>(*a_iter) > get<D>(*b_iter);
    }
};

template <int D>
pair<double, double> minmax(pt_seq_t &seq) {
    double xMin = get<D>(seq[0]);
    double xMax = get<D>(seq[0]);

    // return the iterator(s) that point to the points in seq that have min/max
    // value in the dimension D.
    for (pt_seq_t::iterator iter = seq.begin(); iter < seq.end(); ++iter) {
        double xCoord = get<D>(*iter);
        if (xCoord < xMin)
            xMin = xCoord;

        if (xCoord > xMax)
            xMax = xCoord;
    }

    return pair<double, double>(xMin, xMax);
}

template <int D, int D_prime>
pt_seq_t leftCircularShift(pt_seq_t &pts) {
    // left circular shift so that the leftmost point of
    // the polygon is in the first position, followed by
    // vertices in counterclockwise traversal. Assume poly
    // is in counterclockwise traversal.
    pair<double, double> result = minmax<D>(pts);
    double minVal = result.first;

    vector<pt_iter_t> min_ele_iter_vec;
    for (pt_iter_t iter = pts.begin(); iter != pts.end(); ++iter) {
        if (get<D>(*iter) == minVal)
            min_ele_iter_vec.push_back(iter);
    }

    // sort the iterators in decreasing order according to the y-coordinates of
    // the points pointed to by the iterators. We do this so that in the return 
    // vector, we know which point pointed to by the iterators is on top and which
    // is on bottom.
    sort(min_ele_iter_vec.begin(), min_ele_iter_vec.end(), compare<D_prime, false>());

    // take the point with the smallest x-coordinate with the largest y-coordinate
    pt_const_iter_t leftShift = min_ele_iter_vec[0];
    // pt_iter_t leftShiftCopy = leftShift;
    // if (leftShift == poly.begin()) {
    //     leftShift = (get<0>(*leftShift) == get<0>(*(poly.end()-1))) ? poly.end()-1 : leftShift;
    // }
    // else if (get<0>(*leftShift) == (get<0>(*(leftShiftCopy--)))) {
    //     leftShift = leftShiftCopy;
    // }
    cout << "left shift from: " << get<0>(*leftShift) << ", " << get<1>(*leftShift) << endl;
    pt_seq_t shifted_poly;
    
    pt_const_iter_t vtx_iter;
    for (vtx_iter = leftShift; vtx_iter < pts.end(); vtx_iter++) {
        shifted_poly.push_back(*vtx_iter);
    }
    for (vtx_iter = pts.begin(); vtx_iter < leftShift; vtx_iter++) {
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
    if (num_of_pts < 3)
        throw invalid_argument("received a polygon with less than 3 vertices..");
    // left circular shift so that the leftmost vertex is in the first position
    for (size_t i = 0; i < num_of_pts; ++i) {
        vertices.push_back(make<point_2d>(coords[i][0], coords[i][1]));
    }

    // the dimension in which to perform left circular shift
    const int D = 0;

    // the dimension in which to compare points that have same value in D
    // to find the top left vertex
    const int D_prime = 1;
    vertices = leftCircularShift<D, D_prime>(vertices);

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
    pair<double, double> result = minmax<1>(vertices);

    // create vertical line that span the height of the polygon
    segment_2d scanLine = segment_2d(
        point_2d(xCoord, result.first),
        point_2d(xCoord, result.second));
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

vector<IEvent *> Polygon::generateEvents(string startEventName, string endEventName) {
    pair<double, double> result = minmax<0>(vertices);
    double minX = result.first;
    double maxX = result.second;

    vector<IEvent *> events;

    for (pt_const_iter_t iter = vertices.begin(); iter != vertices.end()-1; ++iter) {
        const point_2d *point_ptr = &(*iter);
        size_t vtx_idx = iter - vertices.begin();

        string eventName;
        enum IEvent::eventTypes eventType;
        double eventX = get<0>(*iter);
        if (eventX == minX) {
            eventName = startEventName;
            eventType = IEvent::BEGIN;
        }
        else if (eventX == maxX) {
            eventName = endEventName;
            eventType = IEvent::END;
        }
        else {
            eventName = "MIDDLE event";
            eventType = IEvent::MIDDLE;
        }

        size_t next_edge_idx = vtx_idx;
        size_t prev_edge_idx = (vtx_idx == 0) ? vertices.size()-1 : vtx_idx - 1;
        PtEvent *newPtEventPtr = new PtEvent(eventType, eventName, point_ptr, 
                &(edges[prev_edge_idx]), &(edges[next_edge_idx]));
        events.push_back(newPtEventPtr);

        // check if the event object is an edge
        if (iter - vertices.begin() > 0 && (get<0>(*iter) == get<0>(*(iter-1)))) {
            // pop last point event and add a new edge event
            events.pop_back();
            events.pop_back();
            // get the edge referred between current and previous points
            size_t event_edge_idx = prev_edge_idx;
            prev_edge_idx = (event_edge_idx == 0) ? edges.size()-1 : event_edge_idx - 1;
            const segment_2d *edge_ptr = &edges[event_edge_idx];

            SegEvent *newSegEventPtr = new SegEvent(eventType, eventName, edge_ptr,
                &edges[prev_edge_idx], &edges[next_edge_idx]);
            events.push_back(newSegEventPtr);
        }
    }

    return events;
}

vector<IEvent *> SurveyArea::generateEvents() {
    return Polygon::generateEvents("BEGIN event", "END event");
}

vector<IEvent *> Obstacle::generateEvents() {
    return Polygon::generateEvents("IN event", "OUT event");
}

int main(void)
{
    double coords[4][2] = {  {0, 0}, {400, 0}, {400, 400}, {0, 400} };
    int num_of_pts = sizeof(coords)/sizeof(coords[0]);
    double (*coords_ptr)[2] = coords;
    SurveyArea sa(coords_ptr, num_of_pts);

    // intersect_seq_t intersects = sa.line_intersect(200);
    // cout << "intersects: " << endl;
    // printSeq(intersects);

    vector<IEvent *> events = sa.generateEvents();
    // print events of survey area
    for (vector<IEvent *>::const_iterator iter = events.begin(); iter != events.end(); ++iter) {
        // dynamic cast to PtEvent or SegEvent
        PtEvent *ptEventPtr = dynamic_cast<PtEvent *>(*iter);
        SegEvent *segEventPtr = dynamic_cast<SegEvent *>(*iter);

        cout << "eventType: " << eventTypeNames[(**iter).eventType] << endl;
        cout << "eventName: " << (**iter).eventName << endl;

        if (ptEventPtr)
            cout << "eventObj(point): " << dsv(*ptEventPtr->eventObj) << endl;
        else if (segEventPtr)
            cout << "eventObj(segment): " << dsv(*segEventPtr->eventObj) << endl;
        cout << "prevEdge: " << dsv((**iter).prevEdge) << endl;
        cout << "nextEdge: " << dsv((**iter).nextEdge) << endl;
        cout << "\n";
        free(*iter);
    }

    return 0;
}