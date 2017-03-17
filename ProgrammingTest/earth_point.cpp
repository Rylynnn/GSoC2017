#include <iostream>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/earth_point.hpp>

namespace bg = boost::geometry;

int main(){
    bg::model::d3::earth_point<double> p1(108,56,23.18,'E',34,20,33.2,'N',385);
    bg::model::d3::earth_point<double> p2(108.94,34.3426,385,"Angle");
    bg::model::d3::earth_point<double> p3(1.90136,0.599931,385,"Redian");
    bg::model::d3::earth_point<double> p4(25,45,57.62,'W',54,23,38.26,'N',-5900);

    p1.set_lat(54.23406448,"Angle");
    p1.set_long(-0.449702,"Redian");
    p1.set_H(-986.7);

    std::cout<<p1.get_lat()<<' '<<p1.get_long()<<' '<<p1.get_H()<<' '<<p1.x()<<' '<<p1.y()<<' '<<p1.z()<<std::endl;

    bg::model::d3::earth_point<double> p5(109,9,41.37,'E',33,34,52.97,'N',1790);

    std::cout<<p1.shortest_distance_Vincenty(p5)<<std::endl;

    return 0;
}
