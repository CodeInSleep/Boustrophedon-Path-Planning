#include <algorithm>
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

typedef vector<BoustrophedonCell *> cell_seq_t;
typedef vector<BoustrophedonCell *>::iterator cell_iter_t;

std::string boolstr(bool v)
{
    return v ? "true" : "false";
}

// configuration of test drone
// struct drone_test_configuration {
//     double altitude;
// };

// function to generate path for a single B. cell
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

template <int D, bool increasing>
struct compareIntersects {
    bool operator()(const intersect_t &inter_a, const intersect_t &inter_b) {
        return (increasing) ? get<D>(inter_a.first) < get<D>(inter_b.first)
                        : get<D>(inter_a.first) > get<D>(inter_b.first);
    }
};

struct compare_xy {
    // sort all events from left to right, bottom to top.
    bool operator()(IEvent *event_a, IEvent *event_b) {
        return (event_a->getxCoord() == event_b->getxCoord()) ?  event_a->getMinY() < event_b->getMinY() :
                    event_a->getxCoord() < event_b->getxCoord();
    }
};



// function to find minimum and maximum item along dimension D
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

    // sort the iterators in decreasing order according to the D_prime dimension of
    // the points pointed to by the iterators. We do this so that in the return 
    // vector, we know which point is on top.
    sort(min_ele_iter_vec.begin(), min_ele_iter_vec.end(), compare<D_prime, false>());

    // take the point with the smallest x-coordinate with the largest y-coordinate
    pt_const_iter_t startPtr = min_ele_iter_vec[0];
    
    cout << "left shift from: " << get<0>(*startPtr) << ", " << get<1>(*startPtr) << endl;

    // store shifted result
    pt_seq_t shifted_poly;
    
    pt_const_iter_t vtx_iter;
    for (vtx_iter = startPtr; vtx_iter < pts.end(); vtx_iter++) {
        shifted_poly.push_back(*vtx_iter);
    }
    for (vtx_iter = pts.begin(); vtx_iter < startPtr; vtx_iter++) {
        shifted_poly.push_back(*vtx_iter);
    }

    return shifted_poly;
}

template <typename T>
void printSeq(const vector<T> &seq) {
    // loop through coordinates
    for (typename vector<T>::const_iterator iter = seq.begin(); iter != seq.end(); ++iter) {
        cout << typeid(seq[0]).name() << dsv(*iter);
    }

    cout << "\n";
}

template<>
void printSeq<intersect_t>(const vector<intersect_t> &seq) {
    for (vector<intersect_t>::const_iterator iter = seq.begin(); iter != seq.end(); ++iter) {
        cout << typeid(seq[0]).name() << dsv((*iter).first) << ", " << dsv(*((*iter).second)) << endl;
    }
}

// TODO: generalize with templating
void addIntersect(const segment_2d *l1, const segment_2d *l2, intersect_seq_t &intersects) {
    // test if there is an intersect. If there is one, add the pointer of the first argument
    // to the intersects vector,
    vector<point_2d> output;
    intersection(*l1, *l2, output);

    if (output.size() != 0) {
        intersect_t intersect = make_pair(output[0], l1);
        intersects.push_back(intersect);
    }
}

segment_2d scanLine(double xCoord, Polygon *poly) {
    pt_seq_t vertices = poly->getVertices();
    // return a scan line that is as long as the largest width of poly
    pair<double, double> mm = minmax<1>(vertices);

    return segment_2d(
        point_2d(xCoord, mm.first),
        point_2d(xCoord, mm.second)
    );
}


// given x-coordinate of scanLine, return the edges of the polygon intersected
// at xCoord and the corresponding points of intersection as pairs
intersect_seq_t line_intersect(Polygon *poly, IEvent *eventPtr) {
    // create vertical line that span the height of the polygon
    const segment_2d line = scanLine(eventPtr->getxCoord(), poly);

    edge_seq_t edges = poly->getEdges();
    intersect_seq_t intersects;
    for (edge_const_iter_t e_iter = edges.begin(); e_iter != edges.end(); ++e_iter) {
        SegEvent *segEventPtr = dynamic_cast<SegEvent *>(eventPtr);

        // not taking vertical edges of eventObj into consideration of intersection
        if (segEventPtr && segEventPtr->eventObj == &(*e_iter))
            continue;

        addIntersect(&(*e_iter), &line, intersects);
    }

    return intersects;
}


intersect_seq_t line_intersect(BoustrophedonCell *cell, IEvent *eventPtr) {
    intersect_seq_t intersects;

    const segment_2d line = scanLine(eventPtr->getxCoord(), cell->within());

    addIntersect(&(*(cell->getFloorSegPtrs().back())), &line, intersects);
    addIntersect(&(*(cell->getCeilSegPtrs().back())), &line, intersects);

    return intersects;
}


BoustrophedonCell::BoustrophedonCell(double leftStartPoint, 
    const segment_2d *startCeilSegPtr, const segment_2d *startFloorSegPtr, SurveyArea *sa)
    : terminated(false), start(leftStartPoint), surveyArea(sa) {
        ceilSegPtrs.push_back(startCeilSegPtr);
        floorSegPtrs.push_back(startFloorSegPtr);
}

void BoustrophedonCell::terminate(double rightEndPoint) {
    // if the cell has not been terminated, then terminate it. If it
    // has been terminated, we do nothing.
    if (!terminated) {
        terminated = true;
        end = rightEndPoint;
    }
}

bool BoustrophedonCell::affected(IEvent *event) {
    if (terminated)
        return false;

    // get the surveyArea that the cell is in
    SurveyArea *parent = within();
    pair<intersect_t, intersect_t> evtLine = parent->eventLine(event);

    // check if the latest floor and ceiling edges of the cells have intersects within
    // the eventLine, terminate them.
    intersect_seq_t cellIntersect = line_intersect(this, event);

    if (cellIntersect.size() != 2)
        return false;

    // check if the openings of opened B. cells are within the eventLine
    return (get<1>(cellIntersect[0].first) >= get<1>(evtLine.first.first) &&
        get<1>(cellIntersect[1].first) <= get<1>(evtLine.second.first));
}

Polygon::Polygon(double coords[][2], int num_of_pts, Polygon::polyTypes polyType)
    : polyType(polyType) {
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

vector<IEvent *> Polygon::generateEvents(string startEventName, string endEventName) {
    pair<double, double> result = minmax<0>(vertices);
    double minX = result.first;
    double maxX = result.second;

    vector<IEvent *> events;

    bool seenMax = false;
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
            seenMax = true;
        }
        else {
            if (seenMax) {
                ostringstream ceilStream;
                ceilStream << "Ceil " << (vertices.end() - iter - 1);
                eventName = ceilStream.str();
                eventType = IEvent::CEILING;
            } else {
                ostringstream floorStream;
                floorStream << "Floor " << (iter - vertices.begin());
                eventName = floorStream.str();
                eventType = IEvent::FLOOR;
            }
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

    // sort all events from left to right, bottom to top.
    sort(events.begin(), events.end(), compare_xy());
    return events;
}

// given xCoord and vector of obstacles generate the eventLine, which spans from the
// closest intersect (either with an obstacle or a border) above the event location to
// the close intersect (either with an obstacle or a border) below the event location
pair<intersect_t, intersect_t> SurveyArea::eventLine(IEvent *event) {
    // only called with IN, OUT events to determine the edges for the new cells; therefore,
    // it is guaranteed that there are a set of intersections above and below the eventObj.
    intersect_seq_t intersects;

    // get all intersecting points of the event line with the polygons and calculate how
    // wide the event line spans (within the border of survey area and without crossing obstacles)
    intersect_seq_t saIntersects = line_intersect(this, event);
    intersects.insert(intersects.end(), saIntersects.begin(), saIntersects.end());
    for (vector<Obstacle *>::iterator iter = obstacles.begin(); iter != obstacles.end(); ++iter) {
        intersect_seq_t obIntersects = line_intersect(*iter, event);
        intersects.insert(intersects.end(), obIntersects.begin(), obIntersects.end());
    }

    sort(intersects.begin(), intersects.end(), compareIntersects<1, true>());

    double eventMinY = event->getMinY();
    double eventMaxY = event->getMaxY();

    vector<intersect_t>::reverse_iterator lowIntersectIter;
    vector<intersect_t>::iterator highIntersectIter;
    
    // scanning for the intersect with the largest y coordinate from the bottom up
    for (highIntersectIter = intersects.begin(); highIntersectIter != intersects.end() &&
        get<1>((*highIntersectIter).first) < eventMaxY; ++highIntersectIter)

    // scanning for the intersect with the smallest y coordinate from the top down
    for (lowIntersectIter = intersects.rbegin(); lowIntersectIter != intersects.rend() &&
        get<1>((*lowIntersectIter).first) > eventMinY; ++lowIntersectIter)

    // special cases where the event has the maximum or minimum coordinates of intersects
    if (event->eventType == IEvent::CEILING)
        lowIntersectIter++;
    else if (event->eventType == IEvent::FLOOR)
        highIntersectIter++;

    return pair<intersect_t, intersect_t>(*lowIntersectIter, *highIntersectIter);
}

// find new openings to open new B. cells with. "opening lines" = eventLine - eventObjs
vector<BoustrophedonCell *> SurveyArea::openCells(IEvent *event, 
    vector<IEvent *>::iterator botEventIter, vector<IEvent *>::iterator topEventIter) {

    pair<intersect_t, intersect_t> evtLine = eventLine(event);

    // assume events are sorted from top to bottom and are within eventLine
    vector<BoustrophedonCell *> cellsOpened;    

    const segment_2d *topLinePtr = evtLine.second.second;
    const segment_2d *botLinePtr = evtLine.first.second;

    double leftStartPoint = get<0>(evtLine.first.first);
    BoustrophedonCell *newCell;
    vector<IEvent *>::iterator iter;

    for (iter = topEventIter; iter != botEventIter; ++iter) {
        if (iter == topEventIter) {
            newCell = new BoustrophedonCell(leftStartPoint, (*iter)->prevEdge, topLinePtr, this);
            cellsOpened.push_back(newCell);
            continue;
        } 

        // start ceiling is the previous event's next edge and start floor is the current event's 
        // previous edge
        newCell = new BoustrophedonCell(leftStartPoint, (*(iter-1))->nextEdge, (*iter)->prevEdge, this);
        cellsOpened.push_back(newCell);
    }

    newCell = new BoustrophedonCell(leftStartPoint, (*(iter-1))->prevEdge, botLinePtr, this);
    cellsOpened.push_back(newCell);

    return cellsOpened;
}

void SurveyArea::closeCells(IEvent *event) {
    for (cell_iter_t iter = cells.begin(); iter != cells.end(); ++iter) {
        BoustrophedonCell *cell = *iter;
        if (cell->affected(event)) {
            cell->terminate(event->getxCoord());
        }
    }
}

void SurveyArea::updateCells(IEvent *event) {
    // update cells that are affected by the event line
    for (cell_iter_t iter = cells.begin(); iter != cells.end(); ++iter) {
        BoustrophedonCell *cell = *iter;
        if (cell->affected(event)) {
            // update the cell with ceiling or floor edge
            if (event->eventType == IEvent::CEILING)
                cell->addCeilSegPtr(event->prevEdge);
            else if (event->eventType == IEvent::FLOOR)
                cell->addFloorSegPtr(event->nextEdge);
        }
    }
}

bool within(segment_2d seg_a, segment_2d seg_b) {
    // return true if seg_a is completely within seg_b
    return (get<1>(seg_b.first) <= get<1>(seg_a.first) && get<1>(seg_b.second) >= get<1>(seg_a.second));
}

// TODO: check that all obstacles are within the survey area (strong requirement)
SurveyArea::SurveyArea(double (*coords)[2], int num_of_pts)
    : Polygon(coords, num_of_pts, Polygon::OUTER) {}

vector<BoustrophedonCell *> SurveyArea::generateBCells() {
    // generate all events of the surveyArea and the obstacles contained within the survey area
    vector<IEvent *> allEvents;
    vector<IEvent *> saEvents = generateEvents();
    allEvents.insert(allEvents.end(), saEvents.begin(), saEvents.end());
    for (vector<Obstacle *>::iterator iter = obstacles.begin(); iter != obstacles.end(); ++iter) {
        vector<IEvent *> obEvents = (*iter)->generateEvents();
        allEvents.insert(allEvents.end(), obEvents.begin(), obEvents.end());
    }

    // loop through all events
    for (vector<IEvent *>::iterator event_iter = allEvents.begin(); event_iter != allEvents.end();) {
        IEvent *eventPtr = *event_iter;
        IEvent::eventTypes eventType = eventPtr->eventType;
        const segment_2d *prevEdge = eventPtr->prevEdge;
        const segment_2d *nextEdge = eventPtr->nextEdge;
        double xCoord = eventPtr->getxCoord();

        if (eventType == IEvent::BEGIN) {
            cout << "encountered BEGIN event" << endl;
            // begin a new cell with the begin event
            cells.push_back(new BoustrophedonCell(xCoord, prevEdge, nextEdge, this));
            event_iter++;
        } else if (eventType == IEvent::END) {
            cout << "encountered END event" << endl;
            // terminate all remaining cells that are opened
            for (vector<BoustrophedonCell *>::iterator iter = cells.begin(); iter != cells.end(); ++iter)
                (*iter)->terminate(xCoord);
            event_iter++;
        } else if (eventType == IEvent::IN) {
            cout << "encountered IN event" << endl;
            // subtract obstacles that are not events from 
            // close off any opened cells with floor and ceiling edges included
            // in scan line.
            closeCells(*event_iter);

            // test iterator that points to an event one past the events at xCoord
            vector<IEvent *>::iterator event_tmp_iter = event_iter;
            while (++event_tmp_iter != allEvents.end() && 
                (*event_tmp_iter)->getxCoord() == xCoord) {}
            vector<BoustrophedonCell *> cellsOpened = openCells(*event_iter, event_iter, event_tmp_iter);
            cells.insert(cells.end(), cellsOpened.begin(), cellsOpened.end());

            event_iter = event_tmp_iter;
        } else if (eventType == IEvent::OUT) {
            cout << "encountered OUT event" << endl;
            closeCells(*event_iter);
        } else {
            cout << "encountered MIDDLE event" << endl;
            updateCells(*event_iter);
        }
    }

    return cells;
}

SurveyArea::~SurveyArea() {
    for (vector<BoustrophedonCell *>::iterator iter = cells.begin(); iter != cells.end(); ++iter)
        free(*iter);
}

vector<IEvent *> SurveyArea::generateEvents() {
    return Polygon::generateEvents("BEGIN event", "END event");
}

Obstacle::Obstacle(double coords[][2], int num_of_pts)
    : Polygon(coords, num_of_pts, Polygon::INNER) {}

vector<IEvent *> Obstacle::generateEvents() {
    return Polygon::generateEvents("IN event", "OUT event");
}

