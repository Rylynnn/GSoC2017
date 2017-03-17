# GSoC2017 for Boost.Geometry
## Programming Compentency Test
### Source Code
All the source code are in the file called "ProgrammingTest", and the README.md of it is the a rest of testcase which I used to test the correct of my code, all the functions in code is tested by myself so it may have little bugs. Every time you find bugs, it is very kind of you to send it to me using issue pages.

#### earth_point_GreatCircle
earth_point_GreatCircle.hpp is the point on the earth with spherical model. I used [Great circle distance](https://en.wikipedia.org/wiki/Great-circle_distance)to calculate the point on the earth. 

I have three constructed functions, which can construct point with longitude and latitude and support transform to cartisian coordinate. The Great circle distance need cartisian coordinates to calculate with the algorithm complexity O(1).

#### earth_point

earth_point.hpp is the point on the earth with ellipsoidal model. I used [Vincenty formulea](https://en.wikipedia.org/wiki/Vincenty%27s_formulae)to calculate the point on the earth.

I have three constructed functions, which can construct point with geodesic coordinate system(B,L,H)and cartisian coordinate(x,y,z). But it just support geodesic coordinate transform to cartician, the implement of cartisian coordinate to geodesic must use Newton iteration method. The Vincenti formulea need geodesic coordinates to calculate. Newton's method has been successfully used to give rapid convergence for all pairs of input, and the algorithm complexity of this algorithm is hard to estimate.


## Issue Pages
There are the records of my ToDo list of project in GSoC2017 of Boost.Geometry and your issue about my project. :)