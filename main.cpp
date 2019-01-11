#include <iostream>

#include "path_planning.hpp"
#include "Event.hpp"

int main(void)
{
    double coords[][2] = {  {100, 0}, {400, 0}, {500, 400}, {250, 500}, {0, 400} };
    int num_of_pts = sizeof(coords)/sizeof(coords[0]);
    SurveyArea sa(coords, num_of_pts);

    vector<IEvent *> events;
    events = sa.generateEvents();

    for (vector<IEvent *>::iterator iter = events.begin(); iter != events.end(); ++iter) {
        cout << eventTypeNames[(*iter)->eventType] << endl;
        cout << (*iter)->eventName << endl;
        cout << "previous edge: " << dsv((*iter)->prevEdge) << endl;
        cout << "next edge: " << dsv((*iter)->nextEdge) << endl;
    }


    // vector<BoustrophedonCell *> cells = sa.generateBCells();

    // cout << "total of " << cells.size() << " number of cells.."<< endl;
    // // print out info about B. cells
    // for (vector<BoustrophedonCell *>::iterator iter = cells.begin(); iter != cells.end(); ++iter) {
    //     BoustrophedonCell *cell = *iter;
    //     cout << "Cell from (x = " << cell->fromX() << " to " << cell->tillX() << ")" << endl;
    //     cout << "ceilEdges: " << endl;
    //     vector<const segment_2d *> ceilEdges = cell->getCeilSegPtrs();
    //     vector<const segment_2d *> floorEdges = cell->getFloorSegPtrs();

    //     for (vector<const segment_2d *>::const_iterator ceil_iter = ceilEdges.begin(); ceil_iter != ceilEdges.end(); ++ceil_iter)
    //         cout << dsv(**ceil_iter) << endl;

    //     cout << "floorEdges: " << endl;
    //     for (vector<const segment_2d *>::const_iterator floor_iter = floorEdges.begin(); floor_iter != floorEdges.end(); ++floor_iter)
    //         cout << dsv(**floor_iter) << endl;
    // }
    return 0;
}