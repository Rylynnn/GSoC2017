#ifndef BOOST_GEOMETRY_GEOMETRIES_EARTH_POINT_HPP
#define BOOST_GEOMETRY_GEOMETRIES_EARTH_POINT_HPP

#include <cstddef>
#include <cmath>
#include <algorithm>
#include <string>
#include <iostream>
#define PI acos(-1)
#define earth_a (double)6378137.0
#define earth_f (double)1/298.257223563
#define earth_b (double)6356752.314245
#define earth_e12 (earth_a * earth_a - earth_b * earth_b) / (earth_a * earth_a)
#define earth_e22 (earth_a * earth_a - earth_b * earth_b) / (earth_b * earth_b)
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
    inline earth_point(double const& long_degree, double const& long_minute, double const& long_second, std::string const& long_direct, double const& lat_degree, double const&  lat_minute, double const& lat_second, std::string const& lat_direct, double const& H)
    {
        this->height = H;
        this->longitude = (long_degree + long_minute/(double)60+ long_second/(double)3600) * PI /(double)180;
        this->latitude = (lat_degree + lat_minute/(double)60 + lat_second/(double)3600) * PI /(double)180;

        if(long_direct == "W")
            { this->longitude *= -1; }
        if(lat_direct == "S")
            { this->latitude *= -1; }
        double cartesian_x, cartesian_y, cartesian_z;
        cartesian_x = (earth_a / (sqrt(1 - earth_e12 * sin(latitude) * sin(latitude))) + height) * cos(latitude) * cos(longitude);
        cartesian_y = (earth_a / (sqrt(1 - earth_e12 * sin(latitude) * sin(latitude))) + height) * cos(latitude) * sin(longitude);
        cartesian_z = (earth_b / (earth_a * sqrt(1 - earth_e12 * sin(latitude) * sin(latitude))) + height) * sin(latitude);

        this->template set<0>(cartesian_x);
        this->template set<1>(cartesian_y);
        this->template set<2>(cartesian_z);
    }

    inline earth_point(double const& longi, double const& lati, double const& H, std::string const& NumericalType)
    {
        this->height = H;
        if(NumericalType == "Angle"){
            this->longitude = longi * PI /(double)180;
            this->latitude = lati * PI / (double)180;

            double cartesian_x, cartesian_y, cartesian_z;
            cartesian_x = (earth_a / (sqrt(1 - earth_e12 * sin(latitude) * sin(latitude))) + height) * cos(latitude) * cos(longitude);
            cartesian_y = (earth_a / (sqrt(1 - earth_e12 * sin(latitude) * sin(latitude))) + height) * cos(latitude) * sin(longitude);
            cartesian_z = (earth_b / (earth_a * sqrt(1 - earth_e12 * sin(latitude) * sin(latitude))) + height) * sin(latitude);
            this->template set<0>(cartesian_x);
            this->template set<1>(cartesian_y);
            this->template set<2>(cartesian_z);
        }
        else if(NumericalType == "Redian"){
            this->longitude = longi;
            this->latitude = lati;

            double cartesian_x, cartesian_y, cartesian_z;
            cartesian_x = (earth_a / (sqrt(1 - earth_e12 * sin(latitude) * sin(latitude))) + height) * cos(latitude) * cos(longitude);
            cartesian_y = (earth_a / (sqrt(1 - earth_e12 * sin(latitude) * sin(latitude))) + height) * cos(latitude) * sin(longitude);
            cartesian_z = (earth_b / (earth_a * sqrt(1 - earth_e12 * sin(latitude) * sin(latitude))) + height) * sin(latitude);

            this->template set<0>(cartesian_x);
            this->template set<1>(cartesian_y);
            this->template set<2>(cartesian_z);
        }
    }
    /// Constructor with x/y/z values
    inline earth_point(CoordinateType const& x, CoordinateType const& y, CoordinateType const& z)
        : model::point<CoordinateType, 3, CoordinateSystem>(x, y, z)
    {
        double theta;
        theta = atan(this->template get<2>() * earth_a / (sqrt(this->template get<0>() * this->template get<0>() + this->template get<1>() * this->template get<1>()) * earth_b));
        this->latitude = atan((this->template get<2>() + earth_e22 * earth_b * sin(theta) * sin(theta) * sin(theta)) / (sqrt(this->template get<0>() * this->template get<0>() + this->template get<1>() * this->template get<1>()) - earth_e12 * earth_a * cos(theta) * cos(theta) * cos(theta)));
        this->longitude = atan(this->template get<1>() / this->template get<0>());
        this->height = sqrt(this->template get<0>() * this->template get<0>() + this->template get<1>() * this->template get<1>()) / cos(latitude) - earth_a / (sqrt(1 - earth_e12 * sin(latitude) * sin(latitude)));
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
    {    
        return this->longitude; 
    }

    /// Get latitude
    inline double const& get_lat() const
    {
        return this->latitude;
    }

    /// Get height
    inline double const& get_H() const
    {
        return this->height;
    }

    /// Set x-value
    inline void x(CoordinateType const& v)
    {
        this->template set<0>(v);
        double theta;
        theta = atan(this->template get<2>() * earth_a / (sqrt(this->template get<0>() * this->template get<0>() + this->template get<1>() * this->template get<1>()) * earth_b));
        this->latitude = atan((this->template get<2>() + earth_e22 * earth_b * sin(theta) * sin(theta) * sin(theta)) / (sqrt(this->template get<0>() * this->template get<0>() + this->template get<1>() * this->template get<1>()) - earth_e12 * earth_a * cos(theta) * cos(theta) * cos(theta)));
        this->longitude = atan(this->template get<1>() / this->template get<0>());
        this->height = sqrt(this->template get<0>() * this->template get<0>() + this->template get<1>() * this->template get<1>()) / cos(latitude) - earth_a / (sqrt(1 - earth_e12 * sin(latitude) * sin(latitude)));
    }

    /// Set y-value
    inline void y(CoordinateType const& v)
    { 
        this->template set<1>(v);
        double theta;
        theta = atan(this->template get<2>() * earth_a / (sqrt(this->template get<0>() * this->template get<0>() + this->template get<1>() * this->template get<1>()) * earth_b));
        this->latitude = atan((this->template get<2>() + earth_e22 * earth_b * sin(theta) * sin(theta) * sin(theta)) / (sqrt(this->template get<0>() * this->template get<0>() + this->template get<1>() * this->template get<1>()) - earth_e12 * earth_a * cos(theta) * cos(theta) * cos(theta)));
        this->longitude = atan(this->template get<1>() / this->template get<0>());
        this->height = sqrt(this->template get<0>() * this->template get<0>() + this->template get<1>() * this->template get<1>()) / cos(latitude) - earth_a / (sqrt(1 - earth_e12 * sin(latitude) * sin(latitude)));
    }

    /// Set z-value
    inline void z(CoordinateType const& v)
    { 
        this->template set<2>(v);
        double theta;
        theta = atan(this->template get<2>() * earth_a / (sqrt(this->template get<0>() * this->template get<0>() + this->template get<1>() * this->template get<1>()) * earth_b));
        this->latitude = atan((this->template get<2>() + earth_e22 * earth_b * sin(theta) * sin(theta) * sin(theta)) / (sqrt(this->template get<0>() * this->template get<0>() + this->template get<1>() * this->template get<1>()) - earth_e12 * earth_a * cos(theta) * cos(theta) * cos(theta)));
        this->longitude = atan(this->template get<1>() / this->template get<0>());
        this->height = sqrt(this->template get<0>() * this->template get<0>() + this->template get<1>() * this->template get<1>()) / cos(latitude) - earth_a / (sqrt(1 - earth_e12 * sin(latitude) * sin(latitude)));
    }

    /// Set longitude
    inline void set_long(double const& v, std::string const& NumericalType)
    { 
        if(NumericalType == "Redian")
        { this->longitude = v; }
        else if(NumericalType == "Angle")
        { this->longitude = v * PI / 180; }
        double cartesian_x, cartesian_y, cartesian_z;
        cartesian_x = (earth_a / (sqrt(1 - earth_e12 * sin(latitude) * sin(latitude))) + height) * cos(latitude) * cos(longitude);
        cartesian_y = (earth_a / (sqrt(1 - earth_e12 * sin(latitude) * sin(latitude))) + height) * cos(latitude) * sin(longitude);
        cartesian_z = (earth_b / (earth_a * sqrt(1 - earth_e12 * sin(latitude) * sin(latitude))) + height) * sin(latitude);

        this->template set<0>(cartesian_x);
        this->template set<1>(cartesian_y);
        this->template set<2>(cartesian_z);
    }

    /// Set latitude
    inline void set_lat(double const& v, std::string const& NumericalType)
    { 
        if(NumericalType == "Redian")
        { this->latitude = v; }
        else if(NumericalType == "Angle")
        { this->latitude = v * PI / 180; } 
        double cartesian_x, cartesian_y, cartesian_z;
        cartesian_x = (earth_a / (sqrt(1 - earth_e12 * sin(latitude) * sin(latitude))) + height) * cos(latitude) * cos(longitude);
        cartesian_y = (earth_a / (sqrt(1 - earth_e12 * sin(latitude) * sin(latitude))) + height) * cos(latitude) * sin(longitude);
        cartesian_z = (earth_b / (earth_a * sqrt(1 - earth_e12 * sin(latitude) * sin(latitude))) + height) * sin(latitude);

        this->template set<0>(cartesian_x);
        this->template set<1>(cartesian_y);
        this->template set<2>(cartesian_z);
    }

    /// Set height
    inline void set_H(double const& v)
    { 
        this->height = v;
        double cartesian_x, cartesian_y, cartesian_z;
        cartesian_x = (earth_a / (sqrt(1 - earth_e12 * sin(latitude) * sin(latitude))) + height) * cos(latitude) * cos(longitude);
        cartesian_y = (earth_a / (sqrt(1 - earth_e12 * sin(latitude) * sin(latitude))) + height) * cos(latitude) * sin(longitude);
        cartesian_z = (earth_b / (earth_a * sqrt(1 - earth_e12 * sin(latitude) * sin(latitude))) + height) * sin(latitude);

        this->template set<0>(cartesian_x);
        this->template set<1>(cartesian_y);
        this->template set<2>(cartesian_z);
    }

    inline double shortest_distance_Vincenty(earth_point const& ep2)
    {
        //std::cout<<this->longitude<<' '<<this->latitude<<"  "<<ep2.get_long()<<' '<<ep2.get_lat()<<std::endl;     
        double lambda_1, lambda_2, alfa_1, alfa_2, sin_alfa, cos2_alfa;
        double tanU1, tanU2, L, s, arc_lenth, sign, sign_;

        L = ep2.get_long() - this->longitude;
        if (L < -PI) L += 2 * PI;
        if (L > PI) L -= 2 * PI;
        tanU1 = (1 - earth_f) * tan(this->latitude);
        tanU2 = (1 - earth_f) * tan(ep2.get_lat());
        double cosU1, sinU1, cosU2, sinU2, cc_U1U2, ss_U1U2, cs_U1U2, sc_U1U2;
        double sin_sign, cos_sign, sin_arc, cos_arc, cos_2arcm;
        double C;

        cosU1 = 1 / sqrt(1 + tanU1 * tanU1);
        sinU1 = tanU1 * cosU1;
        cosU2 = 1 / sqrt(1 + tanU2 * tanU2);
        sinU2 = tanU2 * cosU2;
        //std::cout<<cosU1<<' '<<cosU2<<' '<<sinU1<<' '<<sinU2<<std::endl;
        cc_U1U2 = cosU1 * cosU2;
        ss_U1U2 = sinU1 * sinU2;
        cs_U1U2 = cosU1 * sinU2;
        sc_U1U2 = sinU1 * cosU2;

        sign = L;
        int temp = 10;

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
            C = earth_f / 16 * cos2_alfa * (4 + earth_f * (4 - 3 * cos2_alfa));
            
            sign_ = sign;
            sign = L + (1 - C) * earth_f * sin_alfa * (arc_lenth + C * sin_arc * (cos_2arcm + C * cos_arc * (-1 + 2 * cos_2arcm * cos_2arcm)));
            //std::cout<<sign<<' '<<sign_<<' '<<std::endl;
            temp--;
        }while(fabs(sign - sign_) > eps && temp > 0);
        if(temp == 0){
            std::cout<<"Formula failed"<<std::endl;
        }
        double squ_u, A, B, delta_arc;

        squ_u = cos2_alfa * ((earth_a/earth_b) * (earth_a/earth_b) - 1);
        A = 1 + squ_u / 16384 * (4096 + squ_u * (-768 + squ_u * (320 - 175 * squ_u)));
        B = (squ_u / 1024) * (256 + squ_u * (-128 + squ_u * (74 - 47 * squ_u)));
        delta_arc = B * sin_arc * (cos_2arcm + 0.25 * B * (cos_arc * (-1 + 2 * cos_2arcm * cos_2arcm) - 1/6 * B * cos_2arcm * (-3 + 4 * sin_arc * sin_arc) * (-3 + 4 * cos_2arcm * cos_2arcm)));
        s = earth_b * A * (arc_lenth - delta_arc);
        s /= 1e3;
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

#endif // BOOST_GEOMETRY_GEOMETRIES_earth_point_HPP