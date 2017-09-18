#include <ctime>

#include <iostream>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/algorithms/distance.hpp>
#include <boost/geometry/algorithms/comparable_distance.hpp>

#include <boost/numeric/conversion/bounds.hpp>

int main()
{
		freopen("testcase_cd.out", "r", stdin);
		freopen("testresult_cd.out","w",stdout);
		double a, b, c, d;
		while(std::cin >> a >> b >> c >> d)
		{
			typedef boost::geometry::model::point
					<float, 2, 
					boost::geometry::cs::spherical_equatorial<boost::geometry::degree>
					> point_type;
			point_type P1(a, b);
			point_type P2(c, d);

			double const earth_radius = 6317.0;
															
			std::clock_t distance_start = std::clock();
			for(int i=0;i<1000;i++){
				double distance_result = boost::geometry::distance
						(P1, P2, 
						boost::geometry::strategy::distance::haversine<float>(earth_radius)
						);
			}
			std::clock_t distance_stop = std::clock();
			double secs_distance = double (distance_stop - distance_start)
									/ (double)CLOCKS_PER_SEC;

			std::clock_t comparable_distance_start = std::clock();
			for(int i=0;i<1000;i++)
			{
				double comparable_distance_result = boost::geometry::comparable_distance
						(P1, P2, 
						 boost::geometry::strategy::distance::haversine<float>(earth_radius)
						);
			}
			std::clock_t comparable_distance_stop = std::clock();
			double secs_comparable_distance = double(comparable_distance_stop - comparable_distance_start)
												/ (double)CLOCKS_PER_SEC;

		
			std::cout << "P!: " << a << ' ' << b << ' ' << std::endl;
			std::cout << "P2: " << c << ' ' << d << ' ' << std::endl;	
			std::cout << "distance:" << secs_distance << "s" << std::endl;
			std::cout << "comparable_distance:" << secs_comparable_distance << "s" << std::endl;
			std::cout << std::endl;
		}
						return 0;
}
