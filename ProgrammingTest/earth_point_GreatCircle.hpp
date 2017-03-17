#ifndef BOOST_GEOMETRY_GEOMETRIES_EARTH_POINT_HPP
#define BOOST_GEOMETRY_GEOMETRIES_EARTH_POINT_HPP

#include <cstddef>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <string>
#define PI acos(-1)
#define earth_r (double)6371000

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
    inline earth_point(double const& long_degree, double const& long_minute, double const& long_second, std::string const& long_direct, double const& lat_degree, double const&  lat_minute, double const& lat_second, std::string const& lat_direct)
    {
        this->longitude = (long_degree + long_minute/(double) 60+ long_second/(double)3600)* PI /(double)180;
        this->latitude = (lat_degree + lat_minute/(double)60 + lat_second/(double)3600)* PI /(double)180;
        
        if(long_direct == "W")
            { this->longitude *= -1; }
        if(lat_direct == "S")
            { this->latitude *= -1; }
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
    {}   

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

    /// Set x-value
    inline void x(CoordinateType const& v)
    {
        this->template set<0>(v);
    }

    /// Set y-value
    inline void y(CoordinateType const& v)
    { 
        this->template set<1>(v);
    }

    /// Set z-value
    inline void z(CoordinateType const& v)
    { 
        this->template set<2>(v);
    }

    /// Set longitude
    inline void set_long(double const& v, std::string const& NumericalType)
    { 
        if(NumericalType == "Redian")
        { this->longitude = v; }
        else if(NumericalType == "Angle")
        { this->longitude = v * PI / 180; }
        double cartesian_x, cartesian_y, cartesian_z;
        cartesian_x = cos(latitude) * cos(longitude);
        cartesian_y = cos(latitude) * sin(longitude);
        cartesian_z = sin(latitude);

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
        cartesian_x = cos(latitude) * cos(longitude);
        cartesian_y = cos(latitude) * sin(longitude);
        cartesian_z = sin(latitude);

        this->template set<0>(cartesian_x);
        this->template set<1>(cartesian_y);
        this->template set<2>(cartesian_z);
    }

    inline double shortest_distance_GreatCircle(earth_point const& ep2)
    {
        double s;
        s = earth_r * acos(this->template get<0>() * ep2.x() + this->template get<1>() * ep2.y() + this->template get<2>() - ep2.z());
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