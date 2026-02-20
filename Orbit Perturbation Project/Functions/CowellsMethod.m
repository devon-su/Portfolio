function dstatedt = CowellsMethod(time, state, mu, cd, area, mass)
w_earth = [0; 0; 72.9211e-6];    
r_earth = 6378;

% vvect
dx = state(4);
dy = state(5);
dz = state(6);

r = state(1:3)';
R = r;

alt = R - r_earth;
rho = expDrag(alt);

vrel = ([dx dy dz] - cross(w_earth, r))*10^3;

a_drag = -(1/2)*cd*(area/mass)*rho*(norm(vrel)^2)*vrel/norm(vrel); %ECI (m/s^2)
a_drag = a_drag*10^-3;      %km/s^2

ddx = -mu*state(1)/R^3 + a_drag(1);
ddy = -mu*state(2)/R^3 + a_drag(2); 
ddz = -mu*state(3)/R^3 + a_drag(3);

dstatedt = [dx; dy; dz; ddx; ddy; ddz];
end