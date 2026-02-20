function [dstate] = gausVoP(time, state, Cd, area, mass, drag, j2, Cr, jd_epoch, SRP, nbody)
w_earth = [0; 0; 72.9211e-6];
mu_earth = 398600;
mu_sun = 132.712e9;
r_earth = 6378;

h = state(1);
ecc = state(2);
theta = state(3);
raan = state(4);
inc = state(5);
omega = state(6);

% From COES into perifocal frame
[rvect, vvect] = COEs2rv(h, ecc, theta, raan, inc, omega, 'idk', mu_earth);
r = norm(rvect);
    
% Build Rotation Matrix
R_hat = rvect/r;
N_hat = cross(rvect,vvect)/norm(cross(rvect,vvect));
T_hat = cross(N_hat,R_hat);
DCM = [R_hat T_hat N_hat];
DCM = DCM';

% drag pert
if drag == 1
    alt = r - 6378;
    rho = (1e9)*expDrag(alt);
    Vrel = vvect - cross(w_earth,rvect);
    acc_drag = -0.5*Cd*area*(1/mass).*rho.*(norm(Vrel)).*(Vrel);
    acc_drag = DCM*acc_drag;
else
    acc_drag = 0;
end

jd = jd_epoch + time/(60*60*24);

[~, ~, r_S] = solar_position(jd); % positional stuff
r_S_SC = r_S - rvect;

% SRP pert
if SRP == 1
    SF = 1367;
    Psr = SF/3*10^8;

    thetaB = acos(r_earth/r);
    thetaA = acos(r_earth/r_S);
    theta = acos(dot(r_S, rvect)/(norm(r_S)*r));

    if (thetaA + thetaB) < theta
        F = 0;
    else
        F = 1;
    end

    acc_SRP = -(1000)*Psr*Cr*area*(1/mass)*(r_S)*(1/norm(r_S))*F;
    acc_SRP = DCM*acc_SRP;
else
    acc_SRP = 0;
end

% j2 pert 
if j2 > 1
    acc_J2 = J2_total(rvect, DCM, mu_earth, r_earth, j2);
else
    acc_J2 = 0;
end

% N-body
if nbody == 1
    q = dot(rvect, (2*r_S - rvect))/(norm(r_S)^2);
    F = q*(q^2 - 3*q + 3)/(1+ ((1-q)^(1.5)))';
    acc_nbody = mu_sun*(F*r_S - rvect)/(norm(r_S_SC)^3);
    acc_nbody = DCM*acc_nbody;
else
    acc_nbody = 0;
end

p = acc_drag + acc_J2 + acc_SRP + acc_nbody;
R = p(1);
T = p(2);
N = p(3);
    
dh = r*T;
    
decc = (h/mu_earth)*sin(theta)*R + (1/(mu_earth*h))*((h^2 + mu_earth*r)*cos(theta) + mu_earth*ecc*r)*T;
    
acc_2b = h/r^2;
dtheta_pert = (1/(ecc*h))*((h^2/mu_earth)*cos(theta)*R - (r + (h^2/mu_earth))*sin(theta)*T);
dtheta =  acc_2b + dtheta_pert;
draan = (r/(h*sin(inc)))*sin(omega + theta)*N;
dinc = (r/h)*cos(omega + theta)*N;
domega = -dtheta_pert- ((r*sin(omega + theta))/(h*tan(inc)))*N;
    
dstate = [dh; decc; dtheta; draan; dinc; domega];
end