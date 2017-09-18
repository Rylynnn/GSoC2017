/*************************************************************************
    > File Name: test_andoyer_second.cpp
    > Author: Rylynnn
    > Mail: jingry0321@gmail.com
    > Created Time: 2017年09月03日 星期日 05时01分37秒
 ************************************************************************/
#include <ctime>
#include <cmath>

#include <iostream>

#include <boost/geometry/core/srs.hpp>
#include <boost/geometry/core/cs.hpp>
#include <boost/geometry/core/access.hpp>
#include <boost/geometry/core/radian_access.hpp>

#include <boost/geometry/util/math.hpp>

#include <boost/geometry/geometries/geometries.hpp>

#include <boost/geometry/algorithms/distance.hpp>
#include <boost/geometry/algorithms/comparable_distance.hpp>
#include <boost/geometry/algorithms/comparable_geographic_distance.hpp>

#include <boost/geometry/formulas/cartesine_distance.hpp>
#include <boost/geometry/formulas/compare_length_andoyer.hpp>
#include <boost/geometry/formulas/compare_length_haversine.hpp>
#include <boost/geometry/formulas/flat_approximation.hpp>
#include <boost/geometry/formulas/compare_length_vincenty.hpp>
#include <boost/geometry/formulas/compare_length_lambert.hpp>
#include <boost/geometry/formulas/compare_length_andoyer_second.hpp>
#define MAX 1000
#define EPS 1e-9

namespace bg = boost::geometry;

typedef bg::model::point
    <double, 2,
    bg::cs::geographic<bg::degree>
    > point_type;

double const earth_f = 1 / 298.257223563;
double const earth_e2 = earth_f * (2 - earth_f);
double const earth_r = 6317.0;
int main()
{
    double lon1, lat1, lon2, lat2, lon3, lat3, lon4, lat4;
    int num = 0, cnt = 0;
    freopen("testcase.in", "r", stdin);
    freopen("testresult_andoyer.out", "w", stdout);
    while(std::cin >> lon1 >> lat1 >> lon2 >> lat2 >> lon3 >> lat3 >> lon4 >> lat4){
        cnt++;
        bg::srs::spheroid<double> earth;

        int cg_distance_result, cg_distance_result_second;

        point_type P1(lon1, lat1);
        point_type P2(lon2, lat2);  
        point_type P3(lon3, lat3);
        point_type P4(lon4, lon4);
        
        cg_distance_result = bg::formula::compare_length_andoyer<double>
                             ::apply(P1, P2, P3, P4, earth);
        cg_distance_result_second = bg::formula::compare_length_andoyer_second<double>
                                    ::apply(P1, P2, P3, P4, earth);
        if(cg_distance_result_second != cg_distance_result){
            num++;
            std::cout << "Andoyer:" 
                      << cg_distance_result << std::endl;
            std::cout << "Andoyer_second: " << cg_distance_result_second << std::endl;
            std::cout << "P1: " << lon1 << ", " << lat1 
                      << " P2: " << lon2 << ", " << lat2 << std::endl;
            std::cout << "P3: " << lon3 << ", " << lat3 
                      << " P4: " << lon4 << ", " << lat4 << std::endl;
            
            std::cout << "consistency:" 
                      << (cg_distance_result_second == cg_distance_result)
                      << std::endl << std::endl;
        }
    }
    std::cout << "The number of testcase:" << cnt << std::endl;
    std::cout << "The number of different result:" << num << std::endl;
    return 0;
}
