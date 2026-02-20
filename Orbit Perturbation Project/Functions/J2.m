function acc_J2 = J2(rvect, DCM, mu_earth, r_earth)
% individual J2 pert
J2 = 1.08262668355e-3;
r = norm(rvect);

aI = ((-3*J2*mu_earth*(r_earth^2)*rvect(1))/(2*(r^5)))*(1 - (5*(rvect(3)^2))/(r^2));
aJ = ((-3*J2*mu_earth*(r_earth^2)*rvect(2))/(2*(r^5)))*(1 - (5*(rvect(3)^2))/(r^2));
aK = ((-3*J2*mu_earth*(r_earth^2)*rvect(3))/(2*(r^5)))*(3 - (5*(rvect(3)^2))/(r^2));

aJ2 = [aI; aJ; aK];
aJ2 = DCM*aJ2;

acc_J2 = aJ2;
end