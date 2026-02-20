%% Aero 452 Project
clc; clear; close all

global mu_earth r_earth
mu_earth = 398600;
r_earth = 6378;
%% Epoch (UTC):
% 08 October 2024 09:29:31
% Eccentricity: 0.0001055
% inclination: 0.0461°
% perigee height: 35782 km
% apogee height: 35791 km
% right ascension of ascending node: 87.6205°
% argument of perigee: 187.0055°
% revolutions per day: 1.00271054
% mean anomaly at epoch: 217.8251°
% orbit number at epoch: 1821

% initial steps
% propogate to desired time then start maneuvers
% put chaser 100 km awway from target

% ecc.target = 0.0001055;
ecc.target = 1e-9;
inc.target = 0.0461;
r.target = (35782 + 35791)/2 + r_earth; % circularize
RAAN.target = 87.6205;
AoP.target = 187.0055;
T.target = 2*pi*sqrt(r.target^3/mu_earth);
theta.target = 217.8251;
a.target = r.target;

h.target = sqrt(r.target*mu_earth*(1 + ecc.target* cosd(theta.target)));
[rvect_ECI.target, vvect_ECI.target] = COEs2rv(theta.target, ecc.target, a.target, AoP.target, inc.target, RAAN.target, 'idk', mu_earth);
hvect.target = cross(rvect_ECI.target,vvect_ECI.target);

[Q_eci2lvlh] = rv2DCM(rvect_ECI.target', vvect_ECI.target'); % rot matrix 

rvect_LVLH.target = Q_eci2lvlh*rvect_ECI.target';

% target stuff
rho0 = [0; 100; 0]; % bc 100 km away

rvect_LVLH.chaser = rvect_LVLH.target - rho0;
rvect_ECI.chaser = Q_eci2lvlh'*rvect_LVLH.chaser;

v_ECI.chaser = sqrt(mu_earth/norm(rvect_ECI.chaser));
v_ECIdire.chaser = cross(hvect.target, rvect_ECI.chaser)/norm(cross(hvect.target,rvect_ECI.chaser));
vvect_ECI.chaser = v_ECI.chaser * v_ECIdire.chaser;

vvect_LVLH.target = Q_eci2lvlh * vvect_ECI.target';
vvect_LVLH.chaser = Q_eci2lvlh * vvect_ECI.chaser';

[rhovect, drhovect, ddrhovect] = rva_relative(rvect_LVLH.chaser, vvect_LVLH.chaser, rvect_LVLH.target, vvect_LVLH.target);
disp("The chaser is " + norm(rhovect) + " km away!") % HAS TO BE 100 OR ELSE WE'RE NOT STARTING 100 KM AWAY

tspan = [0 T.target];
options = odeset('RelTol', 1e-8, 'AbsTol',1e-8);

state0.target = [rvect_ECI.target, vvect_ECI.target];
[tnew.target, orbit.target] = ode45(@airplane, tspan, state0.target, options, mu_earth);

figure()
earth_sphere(gca)
hold on
plot3(orbit.target(:, 1), orbit.target(:, 2), orbit.target(:, 3), "--k")
plot3(orbit.target(end, 1), orbit.target(end, 2), orbit.target(end, 3), "diamond", "Color", "b")
plot3(rvect_ECI.chaser(1), rvect_ECI.chaser(2), rvect_ECI.chaser(3), "o", "Color", "k");
ylim padded 
xlim padded 
zlim padded
xlabel("Downrange [km]")
ylabel("Altitude [km]")
title("ECI frame: At Start of Mission")
legend("", "Target Orbit", "Target position", "Chaser position", "Location", "best")
% set(gca,'Color','k')
grid on

vbar0 = rvect_LVLH.target(2) - rvect_LVLH.chaser(2);
rbar0 = rvect_LVLH.target(1) - rvect_LVLH.chaser(1);

figure()
plot(0, 0, "diamond", "LineWidth", 2, "Color", "b")
hold on
plot(vbar0, rbar0, "o", "Linewidth", 2, "Color", "k")
ylim padded 
xlim padded 
xlabel("Downrange [km]")
ylabel("Altitude [km]")
title("LVLH Frame: Relative Distance at Start of Mission")
legend("Target", "Chaser", "Location", "best")
grid on

%% Hop 1
period = T.target;
time_mission = 21; % change val to optimize delta V 
time = time_mission*3600; % [hrs] -> [sec]
pos = [0; 30; 0]; % where we want to end up with hop

[hop1.transferburn1, hop1.burn1, hop1.transferburn2, hop1.delV] = cw_rv(rhovect, pos, drhovect, period, time);

tspan = [0 time];
delrho = hop1.transferburn1; % new relative vel @ burn

state = [rhovect' delrho' rvect_ECI.target vvect_ECI.target];
[timenew, statenew] = ode45(@LinearEOM, tspan, state, options, mu_earth);

hop1.rECI_target = [statenew(:,7), statenew(:,8), statenew(:,9)];
hop1.vECI_target = [statenew(:,10), statenew(:,11), statenew(:,12)];
hop1.rhovect = [statenew(:,1), statenew(:,2), statenew(:,3)];
hop1.drhovect = [statenew(:,4), statenew(:,5), statenew(:,6)]; 

hop1.rECI_chaser = hop1.rECI_target - hop1.rhovect;

for i = 1:length(timenew)
    hop1.DCM = rv2DCM(hop1.rECI_target(i, :)', hop1.vECI_target(i, :)');
    hop1.vECI_chaser = (hop1.DCM' * hop1.drhovect(i, :)') + hop1.vECI_target(i, :);
end

hop1.vECI_chaser(end, 1:3) % new rel vel after hop more complicated bc omega term

rho_hop1 = hop1.rhovect(end,1:3);
disp("The chaser is " + norm(rho_hop1) + " km away!") % must be 30 to match criteria

figure()
earth_sphere(gca)
hold on
plot3(orbit.target(:, 1), orbit.target(:, 2), orbit.target(:, 3), "--k") % orbit
plot3(hop1.rECI_target(:,1), hop1.rECI_target(:,2), hop1.rECI_target(:,3), "--g") % not a full orbit which matches with the 21 hrs we chose
plot3(hop1.rECI_target(end,1), hop1.rECI_target(end,2), hop1.rECI_target(end,3), "diamond", "Color", "b") % new loc of target
plot3(hop1.rECI_chaser(end, 1), hop1.rECI_chaser(end, 2), hop1.rECI_chaser(end, 3), "square", "Color", "k")
xlim padded
ylim padded
zlim padded
xlabel("Downrange [km]")
ylabel("Altitude [km]")
title("ECI Frame: After Hop 1")
legend("", "Target Orbit", "Target path during hop", "Target position", "Chaser position", "Location", "best")
grid on

figure()
plot(0, 0, "diamond", "Linewidth", 2, "Color", "b")
hold on
plot(hop1.rhovect(:,2), hop1.rhovect(:,1), "--r", "LineWidth", 2)
plot(hop1.rhovect(end, 2), hop1.rhovect(end,1), "o", "Linewidth", 2, "Color", "k")
ylim padded 
xlim padded 
xlabel("Downrange [km]")
ylabel("Altitude [km]")
title("LVLH Frame: After Hop 1")
legend("Target location", "Hop trajectory", "Chaser Final location", "Location", "northwest")
grid on

disp("Total Delta-V used: " + hop1.delV*1000 + " m/s.")
disp("Time elapsed: " + time_mission + " hours.")

%% Football Hold 1
n.target = (2*pi)/T.target;
a = hop1.rhovect(end, 2);
fb1.xdot0 = a*n.target/2;
fb1.dv = [fb1.xdot0 0 0];
dv0fb1 = hop1.drhovect + fb1.dv;

fb1.delV_onto = norm(dv0fb1);

dr0 = hop1.rhovect(end, 1:3);
dv0 = fb1.dv;

tspan = [0 T.target]; % football hold for 1 period aka 24 hrs

% after hop
rA = hop1.rECI_target(end,1:3)';
vA = hop1.vECI_target(end,1:3)';
rB = hop1.rECI_chaser(end,1:3)';
vB = hop1.vECI_chaser(end,1:3)';

[fb1.rhovect, fb1.drhovect, fb1.ddrhovect] = rva_relative(rA, vA, rB, vB);

fb1.state = [dr0'; dv0'; rA; vA];

[fb1.timenew, fb1.statenew] = ode45(@LinearEOM, tspan, fb1.state, options, mu_earth);

fb1.rECI_target = [fb1.statenew(:,7), fb1.statenew(:,8), fb1.statenew(:,9)];
fb1.vECI_target = [fb1.statenew(:,10), fb1.statenew(:,11), fb1.statenew(:,12)];
fb1.rhovect = [fb1.statenew(:,1), fb1.statenew(:,2), fb1.statenew(:,3)];
fb1.drhovect = [fb1.statenew(:,4), fb1.statenew(:,5), fb1.statenew(:,6)]; 

figure()
plot(0, 0, "diamond", "Linewidth", 2, "Color", "b")
hold on
plot(fb1.rhovect(:, 2), fb1.rhovect(:, 1), "--r", "LineWidth", 2)
plot(hop1.rhovect(end, 2), hop1.rhovect(end,1), "o", "Linewidth", 2, "Color", "k")
ylim padded 
xlim padded 
xlabel("Downrange [km]")
ylabel("Altitude [km]")
title("LVLH Frame: After Football Hold 1")
legend("Target location", "Football Hold Trajectory", "Chaser Final location", "Location", "best")
grid on

%% Functions
function [xx,yy,zz] = earth_sphere(varargin)
% code I found to plot the earth
%
%
%
% EARTH_SPHERE Generate an earth-sized sphere.
%   [X,Y,Z] = EARTH_SPHERE(N) generates three (N+1)-by-(N+1)
%   matrices so that SURFACE(X,Y,Z) produces a sphere equal to 
%   the radius of the earth in kilometers. The continents will be
%   displayed.
%
%   [X,Y,Z] = EARTH_SPHERE uses N = 50.
%
%   EARTH_SPHERE(N) and just EARTH_SPHERE graph the earth as a 
%   SURFACE and do not return anything.
%
%   EARTH_SPHERE(N,'mile') graphs the earth with miles as the unit rather
%   than kilometers. Other valid inputs are 'ft' 'm' 'nm' 'miles' and 'AU'
%   for feet, meters, nautical miles, miles, and astronomical units
%   respectively.
%
%   EARTH_SPHERE(AX,...) plots into AX instead of GCA.
% 
%  Examples: 
%    earth_sphere('nm') produces an earth-sized sphere in nautical miles
%
%    earth_sphere(10,'AU') produces 10 point mesh of the Earth in
%    astronomical units
%
%    h1 = gca;
%    earth_sphere(h1,'mile')
%    hold on
%    plot3(x,y,z)
%      produces the Earth in miles on axis h1 and plots a trajectory from
%      variables x, y, and z
%   Clay M. Thompson 4-24-1991, CBM 8-21-92.
%   Will Campbell, 3-30-2010
%   Copyright 1984-2010 The MathWorks, Inc. 
% Input Handling
[cax,args,nargs] = axescheck(varargin{:}); % Parse possible Axes input
error(nargchk(0,2,nargs)); % Ensure there are a valid number of inputs
% Handle remaining inputs.
% Should have 0 or 1 string input, 0 or 1 numeric input
j = 0;
k = 0;
n = 50; % default value
units = 'km'; % default value
for i = 1:nargs
    if ischar(args{i})
        units = args{i};
        j = j+1;
    elseif isnumeric(args{i})
        n = args{i};
        k = k+1;
    end
end
if j > 1 || k > 1
    error('Invalid input types')
end
% Calculations
% Scale factors
Scale = {'km' 'm'  'mile'            'miles'           'nm'              'au'                 'ft';
         1    1000 0.621371192237334 0.621371192237334 0.539956803455724 6.6845871226706e-009 3280.839895};
% Identify which scale to use
try
    myscale = 6378.1363*Scale{2,strcmpi(Scale(1,:),units)};
catch %#ok<*CTCH>
    error('Invalid units requested. Please use m, km, ft, mile, miles, nm, or AU')
end
     
% -pi <= theta <= pi is a row vector.
% -pi/2 <= phi <= pi/2 is a column vector.
theta = (-n:2:n)/n*pi;
phi = (-n:2:n)'/n*pi/2;
cosphi = cos(phi); cosphi(1) = 0; cosphi(n+1) = 0;
sintheta = sin(theta); sintheta(1) = 0; sintheta(n+1) = 0;
x = myscale*cosphi*cos(theta);
y = myscale*cosphi*sintheta;
z = myscale*sin(phi)*ones(1,n+1);
% Plotting
if nargout == 0
    cax = newplot(cax);
    % Load and define topographic data
    load('topo.mat','topo','topomap1');
    % Rotate data to be consistent with the Earth-Centered-Earth-Fixed
    % coordinate conventions. X axis goes through the prime meridian.
    % http://en.wikipedia.org/wiki/Geodetic_system#Earth_Centred_Earth_Fixed_.28ECEF_or_ECF.29_coordinates
    %
    % Note that if you plot orbit trajectories in the Earth-Centered-
    % Inertial, the orientation of the contintents will be misleading.
    topo2 = [topo(:,181:360) topo(:,1:180)]; %#ok<NODEF>
    
    % Define surface settings
    props.FaceColor= 'texture';
    props.EdgeColor = 'none';
    props.FaceLighting = 'phong';
    props.Cdata = topo2;
    % Create the sphere with Earth topography and adjust colormap
    surface(x,y,z,props,'parent',cax)
    colormap(topomap1)
% Replace the calls to surface and colormap with these lines if you do 
% not want the Earth's topography displayed.
%     surf(x,y,z,'parent',cax)
%     shading flat
%     colormap gray
    
    % Refine figure
    axis equal
    xlabel(['X [' units ']'])
    ylabel(['Y [' units ']'])
    zlabel(['Z [' units ']'])
    view(127.5,30)
else
    xx = x; yy = y; zz = z;
end
end

function [r, v] = COEs2rv(theta, ecc, a, AoP, inc, RAAN, h, mu)
        theta = theta * pi/180;
        AoP = AoP * pi/180;
        inc = inc * pi/180;
        RAAN = RAAN * pi/180;
        
        C_eci2pqw = Cz(AoP)*Cx(inc)*Cz(RAAN);
        C_pqw2eci = C_eci2pqw';    
    if h == 'idk'
        r_p = a*(1-ecc);

        h_PQW = sqrt(mu*r_p*(1+ecc));
        r_PQW = h_PQW^2/mu / (1+ecc*cos(theta)).*[cos(theta); sin(theta); 0];
        v_PQW = mu/h_PQW.*[-sin(theta); ecc + cos(theta); 0];
    
        r = C_pqw2eci * r_PQW;
        v = C_pqw2eci * v_PQW;
    
        r = r';
        v = v';
    elseif a == 'idk'
        r_PQW = ((h^2)/mu)*(1/(1+ecc*cos(theta)))*[cos(theta); sin(theta); 0];
        v_PQW = (mu/h)*[-sin(theta); ecc + cos(theta); 0];

        r = C_pqw2eci * r_PQW;
        v = C_pqw2eci * v_PQW;
    
        % r = r';
        % v = v';
    end
end

function R = Cx(theta)
R = [1 0 0;
        0 cos(theta) sin(theta);
        0 -sin(theta) cos(theta)];
end

function R = Cz(theta)
R = [cos(theta) sin(theta) 0;
        -sin(theta) cos(theta) 0;
        0 0 1];
end

function [Q_eci2lvlh] = rv2DCM(r, v)
    h = cross(r, v); 

    i_hat = r / norm(r); 
    k_hat = h / norm(h);
    j_hat = cross(k_hat, i_hat);

    Q_eci2lvlh = [i_hat(1), i_hat(2), i_hat(3);
                  j_hat(1), j_hat(2), j_hat(3);
                  k_hat(1), k_hat(2), k_hat(3)];
end

function d_state = airplane(time, state, mu_earth)

% EQ of motion for a 2 body system

x = state(1);
y = state(2);
z = state(3);

dx = state(4);
dy = state(5);
dz = state(6);

r = norm([x y z]);

ddx = -mu_earth*x/r^3;
ddy = -mu_earth*y/r^3;
ddz = -mu_earth*z/r^3;

d_state = [dx; dy; dz; ddx; ddy; ddz]; % Must be a column vector
end

function [rhovect, drhovect, ddrhovect] = rva_relative(rA, vA, rB, vB)
    mu = 398600;

    hA = cross(rA,vA);

    
    i_hat = rA / norm(rA); 
    k_hat = hA / norm(hA); 
    j_hat = cross(k_hat,i_hat);
    
    QXx = [i_hat'; j_hat'; k_hat'];
    
    Omega = hA / (norm(rA)^2);
    Omega_dot = -2*dot(rA,vA) / norm(rA)^2 * Omega;

    aA = -mu * rA / norm(rA)^3;
    aB = -mu * rB / norm(rB)^3;

    r_rel = rB - rA;
    v_rel = vB - vA - cross(Omega,r_rel);

    a_rel = aB - aA - cross(Omega_dot,r_rel) - cross(Omega,cross(Omega,r_rel))- 2 * cross(Omega,v_rel);

    rhovect = QXx*r_rel;
    drhovect = QXx*v_rel;
    ddrhovect = QXx*a_rel;
end

function [dv0_plus,DV_0,DV_f,DV_total] = cw_rv(dr,drf,dv0,period,t)
% CIRCULAR ORBIS ONLY
n = 2*pi / period;
% Matrix Method
phiRR = [4 - 3*cos(n*t),         0,     0;
         6*(sin(n*t) - (n*t)),   1,     0;
         0,                      0,     cos(n*t)];

phiRV = [(1/n)*sin(n*t),         (2/n)*(1-cos(n*t)) ,           0;
         (2/n)*(cos(n*t) - 1),   (1/n)*(4*sin(n*t)-(3*n*t)),    0;
         0,                      0,                             (1/n)*sin(n*t)];

phiVR = [3*n*sin(n*t),         0,           0;
         6*n*(cos(n*t) - 1),   0,           0;
         0,                    0,           -n*sin(n*t)];

phiVV = [cos(n*t),          2*sin(n*t) ,            0;
         -2*sin(n*t),       4*cos(n*t)-3,           0;
         0,                      0,                 cos(n*t)];

dv0_PLUS_start_burn = inv(phiRV)*(drf + (-phiRR*dr));
dvf_MINUS_off_burn = phiVR*dr + phiVV*dv0_PLUS_start_burn;

dv0_plus = dv0_PLUS_start_burn;
dv0_minus = dv0;
dvf_minus = dvf_MINUS_off_burn;

DV_0 = dv0_plus - dv0_minus;

DV_f = -dvf_minus;

DV_total = norm(DV_0) + norm(DV_f);
end

function [dydt] = LinearEOM(time, state, mu)
% state(1:3) - delr
% state(4:6) - delv
% state(7:9) - pos of target
% state(10:12) - vel of target

    Rvect = state(7:9);
    R = norm(Rvect);
    vvect = state(10:12);
        
    hvect = cross(Rvect, vvect);
    h = norm(hvect);

    % chaser/deputy stuff
    vc_x = state(4); % velocities
    vc_y = state(5);
    vc_z = state(6);

    ac_x = ((2*mu)/(R^3) + (h^2 / R^4)) * state(1) - (2 * dot(vvect, Rvect)*(h/R^4)) * state(2) + 2*(h/R^2)*vc_y; % accelerations
    ac_y = ((h^2)/(R^4) - (mu / R^3)) * state(2) + (2 * dot(vvect, Rvect)*(h/R^4)) * state(1) - 2*(h/R^2)*vc_x;
    ac_z = -mu/(R^3)*state(3);

    % target/chief stuff
    vt_x = vvect(1); % velocities
    vt_y = vvect(2);
    vt_z = vvect(3);

    at_x = -mu * state(7) / R^3; % accelerations
    at_y = -mu * state(8) / R^3;
    at_z = -mu * state(9) / R^3;

    dydt = [vc_x vc_y vc_z ac_x ac_y ac_z vt_x vt_y vt_z at_x at_y at_z]'; % dydt(1:6) chaser stuff, dydt(7:12) target stuff
end