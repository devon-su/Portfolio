function acc_J6 = J6(rvect, DCM, mu_earth, r_earth)
% individual J6 pert
J6 = 5.40681239107e-7;
r = norm(rvect);

aI = ((-J6*mu_earth*(r_earth^6)*rvect(1))/(16*(r^9)))*(35 - 945*((rvect(3)^2)/(r^2)) + 3465*((rvect(3)^4)/(r^4)) -3003*((rvect(3)^6)/(r^6)));
aJ = ((-J6*mu_earth*(r_earth^6)*rvect(2))/(16*(r^9)))*(35 - 945*((rvect(3)^2)/(r^2)) + 3465*((rvect(3)^4)/(r^4)) -3003*((rvect(3)^6)/(r^6)));
aK = ((-J6*mu_earth*(r_earth^6)*rvect(3))/(16*(r^9)))*(245 - 2205*((rvect(3)^2)/(r^2)) + 4851*((rvect(3)^4)/(r^4)) -3003*((rvect(3)^6)/(r^6)));

aJ6 = [aI; aJ; aK];
aJ6 = DCM*aJ6;

acc_J6 = aJ6;
end