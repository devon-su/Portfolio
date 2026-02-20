function acc_J5 = J5(rvect, DCM, mu_earth, r_earth)
% individual J5 pert
J5 = -2.27296082869e-7;
r = norm(rvect);

aI = ((3*J5*mu_earth*(r_earth^5)*rvect(1)*rvect(3))/(8*(r^9)))*(35- 210*((rvect(3)^2)/(r^2))+ 231*((rvect(3)^4)/(r^4)));
aJ = ((3*J5*mu_earth*(r_earth^5)*rvect(2)*rvect(3))/(8*(r^9)))*(35- 210*((rvect(3)^2)/(r^2))+ 231*((rvect(3)^4)/(r^4)));
aK = ((3*J5*mu_earth*(r_earth^5)*(rvect(3)^2))/(8*(r^9)))*(105- 315*((rvect(3)^2)/(r^2))+ 231*((rvect(3)^4)/(r^4))) - (15*J5*mu_earth*(r_earth^5))/(8*(r^7));

aJ5 = [aI; aJ; aK];
aJ5 = DCM*aJ5;

acc_J5 = aJ5;
end