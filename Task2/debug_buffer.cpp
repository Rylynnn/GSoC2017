/*************************************************************************
    > File Name: debug_buffer.cpp
    > Author: Rylynnn
    > Mail: jingry0321@163.com
    > Created Time: Wed 07 Jun 2017 06:29:35 PM CST
 ************************************************************************/
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/geometries.hpp>


int main()
{
    typedef double coordinate_type;
    typedef boost::geometry::model::d2::point_xy<coordinate_type> point;
    typedef boost::geometry::model::polygon<point> polygon;

    // Declare strategies
    const double buffer_distance = 1.0;
    const int points_per_circle = 36;
    boost::geometry::strategy::buffer::distance_symmetric<coordinate_type> distance_strategy(buffer_distance);
    boost::geometry::strategy::buffer::join_round join_strategy(points_per_circle);
    boost::geometry::strategy::buffer::end_round end_strategy(points_per_circle);
    boost::geometry::strategy::buffer::point_circle circle_strategy(points_per_circle);
    boost::geometry::strategy::buffer::side_straight side_strategy;

    // Declare output
    boost::geometry::model::multi_polygon<polygon> result;

    // Declare/fill a linestring
    boost::geometry::model::linestring<point> ls;
    boost::geometry::read_wkt("LINESTRING(0 0,4 5,7 4,10 6)", ls);

    // Create the buffer of a linestring
    boost::geometry::buffer(ls, result,
                distance_strategy, side_strategy,
                join_strategy, end_strategy, circle_strategy);


    // Declare/fill a multi point
    boost::geometry::model::multi_point<point> mp;
    boost::geometry::read_wkt("MULTIPOINT((3 3),(4 4),(6 2))", mp);

    // Create the buffer of a multi point
    boost::geometry::buffer(mp, result,
                distance_strategy, side_strategy,
                join_strategy, end_strategy, circle_strategy);


    // Declare/fill a multi_polygon
    boost::geometry::model::multi_polygon<polygon> mpol;
    boost::geometry::read_wkt("MULTIPOLYGON(((0 1,2 5,5 3,0 1)),((1 1,5 2,5 0,1 1)))", mpol);

    // Create the buffer of a multi polygon
    boost::geometry::buffer(mpol, result,
                distance_strategy, side_strategy,
                join_strategy, end_strategy, circle_strategy);


    return 0;
}


