cl = 1.0;
clr = 10.0;
factor = 0.01;
z_coord = 0.0;
distance = 1.0;
outer = 160;
Point(1)={0.0,0.0,z_coord,   cl};
Point(2)={outer,0.0,z_coord, clr};
Point(3)={0.0,outer,z_coord, clr};
Point(4)={-outer,0.0,z_coord, clr};
Point(5)={0.0,-outer,z_coord, clr};
Circle(1)={2,1,3};
Circle(2)={3,1,4};
Circle(3)={4,1,5};
Circle(4)={5,1,2};
Line Loop(5)={1,2,3,4};
PI = 3.14159265359;
currentC = 0;

Function Inclusion
	
	ss = 3.0;

 	
	nNodes = 100;
	step = 0.5*PI/nNodes;
	For t In {0:nNodes-1}
		theta = step*t;
		fac = 1.0+epsilon*Cos(ss*(theta-ww));
		xx = r*fac*Cos(theta);
		yy = r*fac*Sin(theta);
		//Printf("Pos: %g %g", xx,yy);
		pA[t] = newp; Point(pA[t])={x_c+xx,y_c+yy,z_coord,factor*cl};
	EndFor
	For t In {0:nNodes-1}
		theta = step*t+0.5*PI;
		fac = 1.0+epsilon*Cos(ss*(theta-ww));
		xx = r*fac*Cos(theta);
		yy = r*fac*Sin(theta);
		//Printf("Pos: %g %g", xx,yy);
		pB[t] = newp; Point(pB[t])={x_c+xx,y_c+yy,z_coord,factor*cl};
	EndFor
	For t In {0:nNodes-1}
		theta=step*t+1.0*PI;
		fac = 1.0+epsilon*Cos(ss*(theta-ww));
		xx = r*fac*Cos(theta);
		yy = r*fac*Sin(theta);
		//Printf("Pos: %g %g", xx,yy);
		pC[t] = newp; Point(pC[t])={x_c+xx,y_c+yy,z_coord,factor*cl};
	EndFor
	For t In {0:nNodes-1}
		theta = step*t+1.5*PI;
		fac = 1.0+epsilon*Cos(ss*(theta-ww));
		xx = r*fac*Cos(theta);
		yy = r*fac*Sin(theta);
		//Printf("Pos: %g %g", xx,yy);
		pD[t] = newp; Point(pD[t])={x_c+xx,y_c+yy,z_coord,factor*cl};
	EndFor
	
	s1 = newc; Spline(s1)={pA[], pB[0] };
	s2 = newc; Spline(s2)={pB[], pC[0] };
	s3 = newc; Spline(s3)={pC[], pD[0] };
	s4 = newc; Spline(s4)={pD[], pA[0] };

	c5 = newc; Line Loop(c5)={s1,s2,s3,s4};
	currentC = c5;

Return


ww = 0.0;
r = 0.8;
epsilon = 0.0;
shift = 0.0;
x_c = 0.0;
y_c = 0.0;

factor = 0.1;
Call Inclusion;
innerBoundaries[0] = currentC;

r = 1.0;
factor = 0.002;
Call Inclusion;
innerBoundaries[1] = currentC;

r = 2.0;
factor = 0.2;
Call Inclusion;
innerBoundaries[2] = currentC;

r = 5.0;
factor = 1.0;
Call Inclusion;
innerBoundaries[3] = currentC;


Plane Surface(1)={5,innerBoundaries[3]};
Plane Surface(2)={innerBoundaries[3],innerBoundaries[2]};
Plane Surface(3)={innerBoundaries[2],innerBoundaries[1]};
Plane Surface(4)={innerBoundaries[1],innerBoundaries[0]};
Plane Surface(5)={innerBoundaries[0]};

Physical Surface(6)={1,2,3,4,5};
Mesh.Algorithm = 4;
