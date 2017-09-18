/*************************************************************************
    > File Name: test_task_cd.cpp
    > Author: Rylynnn
    > Mail: jingry0321@163.com
    > Created Time: Tue 13 Jun 2017 03:34:06 AM CST
 ************************************************************************/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <iostream>
using namespace std;
int main()
{
	freopen("tc_begin_equtor_cd.in", "w", stdout);
    double a, b, c, d, e, f, g, h;
    for (int i=90; i>0; i--){
        a = 0;
        b = 0;
        c = 0;
		d = i;
        for (int j = 180; j>0; j--){
		    e = 10;
            f = 0;
            g = 10;
            h = j;
			cout << a << ' ' << b << " " 
                 << c << " " << d << ' ' 
                 << e << ' ' << f << ' ' 
                 << g << ' ' << h << endl;	
        }
	}
    return 0;
}
