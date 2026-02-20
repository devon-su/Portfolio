function acc_J4 = J4(rvect, DCM, mu_earth, r_earth)
% individual J4 pert
J4 = -1.61962159137e-6;
r = norm(rvect);

aI = ((-15*J4*mu_earth*(r_earth^4)*rvect(1))/(8*(r^7)))*(1- ((14*(rvect(3)^2))/(r^2)) + ((21*(rvect(3)^4))/(r^4)));
aJ = ((-15*J4*mu_earth*(r_earth^4)*rvect(2))/(8*(r^7)))*(1- ((14*(rvect(3)^2))/(r^2)) + ((21*(rvect(3)^4))/(r^4)));
aK = ((-15*J4*mu_earth*(r_earth^4)*rvect(3))/(8*(r^7)))*(5- ((70*(rvect(3)^2))/(3*(r^2))) + ((21*(rvect(3)^4))/(r^4)));

aJ4 = [aI; aJ; aK];
aJ4 = DCM*aJ4;

acc_J4 = aJ4;
end