a = 1e-3 ;
lc = a ;

Point(1) = {0,0,0,lc};
Point(2) = {a,0,0,lc};
Point(3) = {0,a,0,lc};

Line(4) = {1,2};
Line(5) = {2,3};
Line(6) = {3,1};

Line Loop(7) = {4,5,6};
Plane Surface(8) = {7};
