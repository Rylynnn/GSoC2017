# Test Case
## Get Test
### Test code
```
namespace bg = boost::geometry;

int main(){
    bg::model::d3::earth_point<double> p1(108,56,23.18,'E',34,20,33.2,'N',385);
    bg::model::d3::earth_point<double> p2(108.94,34.3426,385,"Angle");
    bg::model::d3::earth_point<double> p3(1.90136,0.599931,385,"Redian");
    bg::model::d3::earth_point<double> p4(25,45,57.62,'W',54,23,38.26,'N',-5900);
    std::cout<<p1.get_lat()<<' '<<p1.get_long()<<' '<<p1.get_H()<<' '<<p1.x()<<' '<<p1.y()<<' '<<p1.z()<<std::endl;
    std::cout<<p2.get_lat()<<' '<<p2.get_long()<<' '<<p2.get_H()<<' '<<p2.x()<<' '<<p2.y()<<' '<<p2.z()<<std::endl;
    std::cout<<p3.get_lat()<<' '<<p3.get_long()<<' '<<p3.get_H()<<' '<<p3.x()<<' '<<p3.y()<<' '<<p3.z()<<std::endl;
    std::cout<<p4.get_lat()<<' '<<p4.get_long()<<' '<<p4.get_H()<<' '<<p4.x()<<' '<<p4.y()<<' '<<p4.z()<<std::endl;
    return 0;
}
```

### Output
```
0.599391 1.90136 385 -1.71123e+06 4.9868e+06 217.757
0.599391 1.90136 385 -1.71125e+06 4.98679e+06 217.757
0.599931 1.90136 385 -1.71061e+06 4.98496e+06 217.929
0.949354 -0.449702 -5900 3.34854e+06 -1.6163e+06 -4796.12
```

## Set Test
### Test code
```
namespace bg = boost::geometry;

int main(){
    bg::model::d3::earth_point<double> p1(108,56,23.18,'E',34,20,33.2,'N',385);
    p1.set_lat(54.23406448,"Angle");
    std::cout<<p1.get_lat()<<' '<<p1.get_long()<<' '<<p1.get_H()<<' '<<p1.x()<<' '<<p1.y()<<' '<<p1.z()<<std::endl;
    p1.set_long(-0.449702,"Redian");
    std::cout<<p1.get_lat()<<' '<<p1.get_long()<<' '<<p1.get_H()<<' '<<p1.x()<<' '<<p1.y()<<' '<<p1.z()<<std::endl;
    p1.set_H(-986.7);
    std::cout<<p1.get_lat()<<' '<<p1.get_long()<<' '<<p1.get_H()<<' '<<p1.x()<<' '<<p1.y()<<' '<<p1.z()<<std::endl;
    return 0;
}
```

### Output
```
0.946563 1.90136 385 -1.21272e+06 3.53405e+06 313.204
0.946563 -0.449702 385 3.36485e+06 -1.62417e+06 313.204
0.946563 -0.449702 -986.7 3.36413e+06 -1.62382e+06 -799.809
5.89675e-05 -0.449702 -2.64261e+06 3.36413e+06 -1.62382e+06 217.757
```

## Distance Test
### Test code
```
namespace bg = boost::geometry;

int main(){
    bg::model::d3::earth_point<double> p1(108,56,23.18,'E',34,20,33.2,'N',385);
    bg::model::d3::earth_point<double> p2(25,45,57.62,'W',54,23,38.26,'N',-5900);

    std::cout<<p1.shortest_distance_GreatCircle(p2)<<std::endl;

    return 0;
}
```

### Output
```
7.80642e+06
```
