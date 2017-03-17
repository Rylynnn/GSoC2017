#include <cstddef>
#include <cmath>
#include <algorithm>
#include <string>
#include <iostream>
#define PI acos(-1)
#define earth_a (double)6378137.0
#define earth_f (double)1/298.257223563
#define earth_b ((double)1-earth_f) * earth_a
#define eps 1e-12

#include <boost/config.hpp>
#include <boost/mpl/int.hpp>

#include <boost/geometry/core/cs.hpp>
#include <boost/geometry/geometries/point.hpp>

namespace boost { namespace geometry
{

namespace model { namespace d3
{
template<typename CoordinateType, typename CoordinateSystem = cs::cartesian>
class earth_point : public model::point<CoordinateType, 3, CoordinateSystem>
{

private:
    double longitude, latitude, height;
public:

#ifndef BOOST_NO_CXX11_DEFAULTED_FUNCTIONS
    /// \constructor_default_no_init
    earth_point() = default;
#else
    /// \constructor_default_no_init
    inline earth_point()
    {}
#endif
    /// Constructor with longitude and latitude 
    inline earth_point(double const& long_degree, double const& long_minute, double const& long_second, string const& long_direct, double const& lat_degree, double const&  lat_minute, double const& lat_second, string const& lat_direct, double const& height)
    {
        this->longitude = (long_degree + long_minute/(double) 60+ long_second/(double)3600)* PI /(double)180;
        this->latitude = (lat_degree + lat_minute/(double)60 + lat_second/(double)3600)* PI /(double)180;

        double cartesian_x, cartesian_y, cartesian_z;
        cartesian_x = cos(latitude) * cos(longitude);
        cartesian_y = cos(latitude) * sin(longitude);
        cartesian_z = sin(latitude);

        this->template set<0>(cartesian_x);
        this->template set<1>(cartesian_y);
        this->template set<2>(cartesian_z);
    }

    inline earth_point(double const& longitude, double const& latitude, std::string const& NumericalType)
    {
        if(NumericalType == "Angle"){
            this->longitude = longitude * PI /(double)180;
            this->latitude = latitude * PI / (double)180;

            double cartesian_x, cartesian_y, cartesian_z;
            cartesian_x = cos(latitude) * cos(longitude);
            cartesian_y = cos(latitude) * sin(longitude);
            cartesian_z = sin(latitude);

            this->template set<0>(cartesian_x);
            this->template set<1>(cartesian_y);
            this->template set<2>(cartesian_z);
        }
        else if(NumericalType == "Redian"){
            this->longitude = longitude;
            this->latitude = latitude;

            double cartesian_x, cartesian_y, cartesian_z;
            cartesian_x = cos(latitude) * cos(longitude);
            cartesian_y = cos(latitude) * sin(longitude);
            cartesian_z = sin(latitude);

            this->template set<0>(cartesian_x);
            this->template set<1>(cartesian_y);
            this->template set<2>(cartesian_z);
        }
    }
    /// Constructor with x/y/z values
    inline earth_point(CoordinateType const& x, CoordinateType const& y, CoordinateType const& z)
        : model::point<CoordinateType, 3, CoordinateSystem>(x, y, z)
    {

    }

    /// Get x-value
    inline CoordinateType const& x() const
    { return this->template get<0>(); }

    /// Get y-value
    inline CoordinateType const& y() const
    { return this->template get<1>(); }

    /// Get z-value
    inline CoordinateType const& z() const
    { return this->template get<2>(); }

    /// Get longitude
    inline double const& get_long() const
    { return this->longitude; }

    // Get latitude
    inline double const& get_lat() const
    { return this->latitude; }

    /// Set x-value
    inline void x(CoordinateType const& v)
    { this->template set<0>(v); }

    /// Set y-value
    inline void y(CoordinateType const& v)
    { this->template set<1>(v); }

    /// Set z-value
    inline void z(CoordinateType const& v)
    { this->template set<2>(v); }

    /// Set longitude
    inline void set_long(double const& v)
    { this->longitude = v; }

    /// Set latitude
    inline void set_lat(double const& v)
    { this->latitude = v; }

    inline double sortest_distance_Vincenty(earth_point const& ep1, earth_point const& ep2)
    {
        double lambda_1, lambda_2, alfa_1, alfa_2, sin_alfa, cos2_alfa;
        double U1, U2, L, s, arc_lenth, sign, sign_;

        L = ep2.get_long() - ep1.get_long();
        U1 = atan((1 - earth_f) * tan(ep1.get_lat()));
        U2 = atan((1 - earth_f) * tan(ep2.get_lat()));

        double cosU1, sinU1, cosU2, sinU2, cc_U1U2, ss_U1U2, cs_U1U2, sc_U1U2;
        double sin_sign, cos_sign, sin_arc, cos_arc, cos_2arcm;
        double C;

        cosU1 = cos(U1);
        sinU1 = sin(U1);
        cosU2 = cos(U2);
        sinU2 = sin(U2);
        cc_U1U2 = cosU1 * cosU2;
        ss_U1U2 = sinU1 * sinU2;
        cs_U1U2 = cosU1 * sinU2;
        sc_U1U2 = sinU1 * cosU2;

        sign = L;
        int temp = 500;

        do
        {
            sin_sign = sin(sign), cos_sign = cos(sign);

            sin_arc = sqrt((cosU2 * sin_sign) * (cosU2 * sin_sign) + (cs_U1U2 - sc_U1U2 * cos_sign) * (cs_U1U2 - sc_U1U2 * cos_sign));
            if(sin_arc == 0)
            {
                return 0;
            }
            cos_arc = ss_U1U2 + cc_U1U2 * cos_sign;
            arc_lenth = atan(sin_arc / cos_arc);
            sin_alfa = (cc_U1U2 * sin_sign) / sin_arc;
            cos2_alfa = 1 - sin_alfa * sin_alfa;
            if(cos2_alfa == 0){
                cos_2arcm = 0;
            }
            else{
                cos_2arcm = cos_arc - (2 * ss_U1U2 / cos2_alfa);
            }

            C = earth_f / 6 * cos2_alfa * (4 + earth_f * (4 - 3 * cos2_alfa));
            sign_ = sign;
            sign = sign + (1 - C) * earth_f * sin_alfa * (arc_lenth + C * sin_arc * (cos_2arcm + C * cos_arc * (-1 + 2 * cos_2arcm * cos_2arcm)));
        }while(fabs(sign - sign_) > eps && temp > 0);

        double squ_u, A, B, delta_arc;

        squ_u = cos2_alfa * (earth_a * earth_a / earth_b * earth_b - 1) ;
        A = 1 + squ_u / 16384 * (4096 + squ_u * (-768 + squ_u * (320 - 175 * squ_u)));
        B = (squ_u / 1024) * (256 + squ_u * (-128 + squ_u * (74 - 47 * squ_u)));
        delta_arc = B * sin_arc * (cos_2arcm + 0.25 * B * (cos_arc * (-1 + 2 * cos_2arcm * cos_2arcm) - 1/6 * B * cos_2arcm * (-3 + 4 * sin_arc * sin_arc) * (-3 + 4 * cos_2arcm * cos_2arcm)));
        s = earth_b * A * (arc_lenth - delta_arc);

        return s;
    }
};


}} // namespace model::d3


// Adapt the earth_point to the concept
#ifndef DOXYGEN_NO_TRAITS_SPECIALIZATIONS
namespace traits
{

template <typename CoordinateType, typename CoordinateSystem>
struct tag<model::d3::earth_point<CoordinateType, CoordinateSystem> >
{
    typedef point_tag type;
};

template<typename CoordinateType, typename CoordinateSystem>
struct coordinate_type<model::d3::earth_point<CoordinateType, CoordinateSystem> >
{
    typedef CoordinateType type;
};

template<typename CoordinateType, typename CoordinateSystem>
struct coordinate_system<model::d3::earth_point<CoordinateType, CoordinateSystem> >
{
    typedef CoordinateSystem type;
};

template<typename CoordinateType, typename CoordinateSystem>
struct dimension<model::d3::earth_point<CoordinateType, CoordinateSystem> >
    : boost::mpl::int_<3>
{};

template<typename CoordinateType, typename CoordinateSystem, std::size_t Dimension>
struct access<model::d3::earth_point<CoordinateType, CoordinateSystem>, Dimension >
{
    static inline CoordinateType get(
        model::d3::earth_point<CoordinateType, CoordinateSystem> const& p)
    {
        return p.template get<Dimension>();
    }

    static inline void set(model::d3::earth_point<CoordinateType, CoordinateSystem>& p,
        CoordinateType const& value)
    {
        p.template set<Dimension>(value);
    }
};

} // namespace traits
#endif // DOXYGEN_NO_TRAITS_SPECIALIZATIONS

}} // namespace boost::geometry

namespace bg = boost::geometry;

int main(){
    bg::model::d3::earth_point<double> p1(34,20,33.2,108,56,23.18);
    std::cout<<p1.get_lat()<<" "<<p1.get_long()<<std::endl;
    return 0;
}
