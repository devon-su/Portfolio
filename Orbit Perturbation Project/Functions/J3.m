function acc_J3 = J3(rvect, DCM, mu_earth, r_earth)
% individual J3 pert
J3 =  -2.53265648533e-6;

r = norm(rvect);
rx = rvect(1);
ry = rvect(2);
rz = rvect(3);

aI = ((-5*J3*mu_earth*r_earth^3*rx)/(2*r^7)) * (3*rz - (7*rz^3)/(r^2));
aJ = ((-5*J3*mu_earth*r_earth^3*ry)/(2*r^7)) * (3*rz - (7*rz^3)/(r^2));
aK = ((-5*J3*mu_earth*r_earth^3)/(2*r^7)) * (6*rz^2 - (7*rz^4)/(r^2) - (3/5)*r^2);

aJ3 = [aI; aJ; aK];
aJ3 = DCM*aJ3;

acc_J3 = aJ3;
end