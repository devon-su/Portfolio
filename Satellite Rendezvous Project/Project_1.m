%% Aero 452 Project
clc; clear; close all
addpath("C:\Users\Devon\Desktop\meow\AERO\Year 4\Aero 452\Functions")

global mu_earth r_earth
mu_earth = 398600;
r_earth = 6378;
%% Initial Setup
% Epoch (UTC):
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

disp("~~~~~~~~~~~~~~~~~~~~ Mission Start Summary ~~~~~~~~~~~~~~~~~~~~")
disp("The chaser is " + norm(rhovect) + " km away!") % HAS TO BE 100 OR ELSE WE'RE NOT STARTING 100 KM AWAY
disp("Total Delta-V used: " + 0 + " m/s.")
disp("Time elapsed: " + 0 + " hours.")
disp(" ")

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

for i = 1:length(timenew) % new rel vel after hop more complicated bc omega term
    hop1.DCM = rv2DCM(hop1.rECI_target(i, :)', hop1.vECI_target(i, :)');
    hop1.vECI_chaser = (hop1.DCM' * hop1.drhovect(i, :)') + hop1.vECI_target(i, :);
end

rho_hop1 = hop1.rhovect(end,1:3);

% ECI plot
figure()
earth_sphere(gca)
hold on
plot3(hop1.rECI_target(:,1), hop1.rECI_target(:,2), hop1.rECI_target(:,3), "--g") % not a full orbit which matches
plot3(hop1.rECI_target(end,1), hop1.rECI_target(end,2), hop1.rECI_target(end,3), "diamond", "Color", "b") % new loc of target
plot3(hop1.rECI_chaser(end, 1), hop1.rECI_chaser(end, 2), hop1.rECI_chaser(end, 3), "square", "Color", "k")
xlim padded
ylim padded
zlim padded
xlabel("Downrange [km]")
ylabel("Altitude [km]")
title("ECI Frame: After Hop 1")
legend("", "Target path during hop", "Target position", "Chaser position", "Location", "best")
grid on

% LVLH plot
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

delV_mission = hop1.delV;

disp("~~~~~~~~~~~~~~~~~~~~ Hop 1 Summary ~~~~~~~~~~~~~~~~~~~~")
disp("The chaser is " + norm(hop1.rhovect(end, 1:3)) + " km away!")
disp("Total Delta-V used: " + delV_mission*1000 + " m/s.")
disp("Time elapsed: " + time_mission + " hours.")
disp(" ")

%% Football Hold 1
n.target = (2*pi)/T.target;
a = hop1.rhovect(end, 2);
fb1.xdot0 = a*n.target/2;
fb1.dv = [fb1.xdot0 0 0];
dv0fb1 = hop1.drhovect + fb1.dv;

fb1.time = T.target*3;

fb1.delV_onto = norm(dv0fb1); % small bc only one burn onto football trajectory

dr0 = hop1.rhovect(end, 1:3);
dv0 = fb1.dv;

tspan = [0 fb1.time]; % football hold for 1 period aka 24 hrs

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

fb1.rECI_chaser = fb1.rECI_target - fb1.rhovect;
for i = 1:length(fb1.timenew)
    fb1.DCM = rv2DCM(fb1.rECI_target(i, :)', fb1.vECI_target(i, :)');
    fb1.vECI_chaser = (fb1.DCM' * fb1.drhovect(i, :)') + fb1.vECI_target(i, :);
end

time_mission = time_mission + fb1.time/60/60;

% ECI plot
figure()
earth_sphere(gca)
hold on
plot3(fb1.rECI_target(:,1), fb1.rECI_target(:,2), fb1.rECI_target(:,3), "--g") % full orbit bc 1 period football hold
plot3(fb1.rECI_target(end,1), fb1.rECI_target(end,2), fb1.rECI_target(end,3), "diamond", "Color", "b") % new loc of target
plot3(fb1.rECI_chaser(end, 1), fb1.rECI_chaser(end, 2), fb1.rECI_chaser(end, 3), "square", "Color", "k")
xlim padded
ylim padded
zlim padded
xlabel("Downrange [km]")
ylabel("Altitude [km]")
title("ECI Frame: After Football Hold 1")
legend("", "Target path during football hold", "Target position", "Chaser position", "Location", "best")
grid on

% LVLH plot
figure()
plot(0, 0, "diamond", "Linewidth", 2, "Color", "b")
hold on
plot(fb1.rhovect(:, 2), fb1.rhovect(:, 1), "--r", "LineWidth", 2)
plot(fb1.rhovect(end, 2), fb1.rhovect(end,1), "o", "Linewidth", 2, "Color", "k")
ylim padded 
xlim padded 
xlabel("Downrange [km]")
ylabel("Altitude [km]")
title("LVLH Frame: After Football Hold 1")
legend("Target location", "Football Hold Trajectory", "Chaser Final location", "Location", "best")
grid on

delV_mission = delV_mission + fb1.delV_onto/1000;

disp("~~~~~~~~~~~~~~~~~~~~ Football Hold 1 Summary ~~~~~~~~~~~~~~~~~~~~")
disp("The chaser is " + round(norm(fb1.rhovect(end, 1:3))) + " km away!")
disp("Total Delta-V used: " + delV_mission*1000 + " m/s.")
disp("Time elapsed: " + time_mission + " hours.")
disp(" ")

%% Hop 2
hop2.pos = [0; 1; 0]; % want to be 1 km away from target
hop2.time = 18; % change to optimize delV

t = hop2.time*60*60;

[hop2.rhovect, hop2.drhovect, hop2.ddrhovect] = rva_relative(fb1.rECI_chaser(end,1:3)', fb1.vECI_chaser(end,1:3)', fb1.rECI_target(end,1:3)', fb1.vECI_target(end,1:3)'); 

hop2.dr = fb1.rhovect(end, 1:3)';
hop2.dv0 = fb1.drhovect(end, 1:3);

[hop2.dv0_plus, hop2.DV_0, hop2.DV_f, hop2.DV_total] = cw_rv(hop2.dr, hop2.pos, hop2.dv0, period, t);

tspan = [0 t];
hop2.dv = hop2.dv0_plus;
hop2.state = [hop2.dr; hop2.dv; fb1.rECI_target(end, 1:3)'; fb1.vECI_target(end, 1:3)'];

[hop2.timenew, hop2.statenew] = ode45(@LinearEOM, tspan, hop2.state, options, mu_earth);

hop2.rECI_target = [hop2.statenew(:,7), hop2.statenew(:,8), hop2.statenew(:,9)];
hop2.vECI_target = [hop2.statenew(:,10), hop2.statenew(:,11), hop2.statenew(:,12)];
hop2.rhovect = [hop2.statenew(:,1), hop2.statenew(:,2), hop2.statenew(:,3)];
hop2.drhovect = [hop2.statenew(:,4), hop2.statenew(:,5), hop2.statenew(:,6)];

hop2.rECI_chaser = hop2.rECI_target - hop2.rhovect;
for i = 1:length(hop2.timenew)
    hop2.DCM = rv2DCM(hop2.rECI_target(i, :)', hop2.vECI_target(i, :)');
    hop2.vECI_chaser = (hop2.DCM' * hop2.drhovect(i, :)') + hop2.vECI_target(i, :);
end

delV_mission = delV_mission + hop2.DV_total;
time_mission = time_mission + hop2.time;

% ECI plot
figure()
earth_sphere(gca)
hold on
plot3(hop2.rECI_target(:,1), hop2.rECI_target(:,2), hop2.rECI_target(:,3), "--g") % not full orbit based on time chosen
plot3(hop2.rECI_target(end,1), hop2.rECI_target(end,2), hop2.rECI_target(end,3), "diamond", "Color", "b") % new loc of target
plot3(hop2.rECI_chaser(end, 1), hop2.rECI_chaser(end, 2), hop2.rECI_chaser(end, 3), "square", "Color", "k")
xlim padded
ylim padded
zlim padded
xlabel("Downrange [km]")
ylabel("Altitude [km]")
title("ECI Frame: After Hop 2")
legend("", "Target path during hop 2", "Target position", "Chaser position", "Location", "best")
grid on

% LVLH plot
figure()
plot(0, 0, "diamond", "Linewidth", 2, "Color", "b")
hold on
plot(hop2.rhovect(:, 2), hop2.rhovect(:, 1), "--r", "LineWidth", 2)
plot(hop2.rhovect(end, 2), hop2.rhovect(end,1), "o", "Linewidth", 2, "Color", "k")
ylim padded 
xlim padded 
xlabel("Downrange [km]")
ylabel("Altitude [km]")
title("LVLH Frame: After Hop 2")
legend("Target location", "Hop 2 trajectory", "Chaser Final location", "Location", "best")
grid on

disp("~~~~~~~~~~~~~~~~~~~~ Hop 2 Summary ~~~~~~~~~~~~~~~~~~~~")
disp("The chaser is " + norm(hop2.rhovect(end, 1:3)) + " km away!")
disp("Total Delta-V used: " + delV_mission*1000 + " m/s.")
disp("Time elapsed: " + time_mission + " hours.")
disp(" ")

%% V-bar hold 1
t_hold = 24*2;

vb1.dr = [0; 1; 0];
vb1.pos = vb1.dr;
vb1.dv0 = hop2.drhovect(end, 1:3);
vb1.dv = vb1.dv0;

t = t_hold*60*60;

[vb1.dv0_plus, vb1.DV_0, vb1.DV_f, vb1.DV_total] = cw_rv(vb1.dr, vb1.pos, vb1.dv0, period, t);

tspan = [0 t];
dv = vb1.dv0_plus;

vb1.state = [vb1.dr; dv; hop2.rECI_chaser(end, 1:3)'; hop2.vECI_chaser(end, 1:3)'];
[vb1.timenew, vb1.statenew] = ode45(@circularLinearEOM, tspan, vb1.state, options, T.target, mu_earth);

vb1.rECI_target = [vb1.statenew(:,7), vb1.statenew(:,8), vb1.statenew(:,9)];
vb1.vECI_target = [vb1.statenew(:,10), vb1.statenew(:,11), vb1.statenew(:,12)];
vb1.rhovect = [vb1.statenew(:,1), vb1.statenew(:,2), vb1.statenew(:,3)];
vb1.drhovect = [vb1.statenew(:,4), vb1.statenew(:,5), vb1.statenew(:,6)];

vb1.rECI_chaser = vb1.rECI_target - vb1.rhovect;
for i = 1:length(vb1.timenew)
    vb1.DCM = rv2DCM(vb1.rECI_target(i, :)', vb1.vECI_target(i, :)');
    vb1.vECI_chaser = (vb1.DCM' * vb1.drhovect(i, :)') + vb1.vECI_target(i, :);
end

delV_mission = delV_mission + vb1.DV_total/1000;
time_mission = time_mission + t_hold;

% ECI plot
figure()
earth_sphere(gca)
hold on
plot3(vb1.rECI_target(:,1), vb1.rECI_target(:,2), vb1.rECI_target(:,3), "--g") % full orbit bc 24 hour 
plot3(vb1.rECI_target(end,1), vb1.rECI_target(end,2), vb1.rECI_target(end,3), "diamond", "Color", "b") % new loc of target
plot3(vb1.rECI_chaser(end, 1), vb1.rECI_chaser(end, 2), vb1.rECI_chaser(end, 3), "square", "Color", "k")
xlim padded
ylim padded
zlim padded
xlabel("Downrange [km]")
ylabel("Altitude [km]")
title("ECI Frame: After V-bar Hold 1")
legend("", "Target path during V-bar hold", "Target position", "Chaser position", "Location", "best")
grid on

% LVLH plot
figure()
plot(0, 0, "diamond", "Linewidth", 2, "Color", "b")
hold on
plot(vb1.rhovect(:, 2), vb1.rhovect(:, 1), "--r", "LineWidth", 2)
plot(vb1.rhovect(end, 2), vb1.rhovect(end, 1), "o", "Linewidth", 2, "Color", "k")
ylim padded 
xlim padded 
xlabel("Downrange [km]")
ylabel("Altitude [km]")
title("LVLH Frame: After V-bar hold 1")
legend("Target location", "V-bar hold trajectory", "Chaser Final location", "Location", "best")
grid on

disp("~~~~~~~~~~~~~~~~~~~~ V-bar Hold 1 Summary ~~~~~~~~~~~~~~~~~~~~")
disp("The chaser is " + norm(vb1.rhovect(1, 1:3)) + " km away!")
disp("Total Delta-V used: " + delV_mission*1000 + " m/s.")
disp("Time elapsed: " + time_mission + " hours.")
disp(" ")

%% Hop 3
hop3.pos = [0; 0.3; 0];
hop3.time = 24;

t = hop3.time*60*60;

[hop3.rhovect, hop3.drhovect, hop3.ddrhovect] = rva_relative(hop2.rECI_chaser(end,1:3)', hop2.vECI_chaser(end,1:3)', hop2.rECI_target(end,1:3)', hop2.vECI_target(end,1:3)'); 

hop3.dr = hop2.rhovect(end, 1:3)';
hop3.dv0 = hop2.drhovect(end, 1:3);

[hop3.dv0_plus, hop3.DV_0, hop3.DV_f, hop3.DV_total] = cw_rv(hop3.dr, hop3.pos, hop3.dv0, period, t);

tspan = [0 t];
hop3.dv = hop3.dv0_plus;
hop3.state = [hop3.dr; hop3.dv; hop2.rECI_target(end, 1:3)'; hop2.vECI_target(end, 1:3)'];

[hop3.timenew, hop3.statenew] = ode45(@LinearEOM, tspan, hop3.state, options, mu_earth);

hop3.rECI_target = [hop3.statenew(:,7), hop3.statenew(:,8), hop3.statenew(:,9)];
hop3.vECI_target = [hop3.statenew(:,10), hop3.statenew(:,11), hop3.statenew(:,12)];
hop3.rhovect = [hop3.statenew(:,1), hop3.statenew(:,2), hop3.statenew(:,3)];
hop3.drhovect = [hop3.statenew(:,4), hop3.statenew(:,5), hop3.statenew(:,6)];

hop3.rECI_chaser = hop3.rECI_target - hop3.rhovect;
for i = 1:length(hop3.timenew)
    hop3.DCM = rv2DCM(hop3.rECI_target(i, :)', hop3.vECI_target(i, :)');
    hop3.vECI_chaser = (hop3.DCM' * hop3.drhovect(i, :)') + hop3.vECI_target(i, :);
end

delV_mission = delV_mission + hop3.DV_total;
time_mission = time_mission + hop3.time;

% ECI plot
figure()
earth_sphere(gca)
hold on
plot3(hop3.rECI_target(:,1), hop3.rECI_target(:,2), hop3.rECI_target(:,3), "--g") % not full orbit based on time chosen
plot3(hop3.rECI_target(end,1), hop3.rECI_target(end,2), hop3.rECI_target(end,3), "diamond", "Color", "b") % new loc of target
plot3(hop3.rECI_chaser(end, 1), hop3.rECI_chaser(end, 2), hop3.rECI_chaser(end, 3), "square", "Color", "k")
xlim padded
ylim padded
zlim padded
xlabel("Downrange [km]")
ylabel("Altitude [km]")
title("ECI Frame: After Hop 3")
legend("", "Target path during hop 3", "Target position", "Chaser position", "Location", "best")
grid on

% LVLH plot
figure()
plot(0, 0, "diamond", "Linewidth", 2, "Color", "b")
hold on
plot(hop3.rhovect(:, 2)*1000, hop3.rhovect(:, 1)*1000, "--r", "LineWidth", 2)
plot(hop3.rhovect(end, 2)*1000, hop3.rhovect(end,1)*1000, "o", "Linewidth", 2, "Color", "k")
ylim padded 
xlim padded 
xlabel("Downrange [m]")
ylabel("Altitude [m]")
title("LVLH Frame: After Hop 3")
legend("Target location", "Hop 3 trajectory", "Chaser Final location", "Location", "best")
grid on

disp("~~~~~~~~~~~~~~~~~~~~ Hop 3 Summary ~~~~~~~~~~~~~~~~~~~~")
disp("The chaser is " + norm(hop3.rhovect(end, 1:3)) + " km away!")
disp("Total Delta-V used: " + delV_mission*1000 + " m/s.")
disp("Time elapsed: " + time_mission + " hours.")
disp(" ")

%% V-bar hold 2
t_hold = 24*2 + 12.3957 - 3.2 - 11.8;

vb2.dr = [0; 0.3; 0];
vb2.pos = vb2.dr;
vb2.dv0 = hop3.drhovect(end, 1:3);
vb2.dv = vb2.dv0;

t = t_hold*60*60;

[vb2.dv0_plus, vb2.DV_0, vb2.DV_f, vb2.DV_total] = cw_rv(vb2.dr, vb2.pos, vb2.dv0, period, t);

tspan = [0 t];
dv = vb2.dv0_plus;

vb2.state = [vb2.dr; dv; hop3.rECI_chaser(end, 1:3)'; hop3.vECI_chaser(end, 1:3)'];
[vb2.timenew, vb2.statenew] = ode45(@circularLinearEOM, tspan, vb2.state, options, period, mu_earth);

vb2.rECI_target = [vb2.statenew(:,7), vb2.statenew(:,8), vb2.statenew(:,9)];
vb2.vECI_target = [vb2.statenew(:,10), vb2.statenew(:,11), vb2.statenew(:,12)];
vb2.rhovect = [vb2.statenew(:,1), vb2.statenew(:,2), vb2.statenew(:,3)];
vb2.drhovect = [vb2.statenew(:,4), vb2.statenew(:,5), vb2.statenew(:,6)];

vb2.rECI_chaser = vb2.rECI_target - vb2.rhovect;
for i = 1:length(vb2.timenew)
    vb2.DCM = rv2DCM(vb2.rECI_target(i, :)', vb2.vECI_target(i, :)');
    vb2.vECI_chaser = (vb2.DCM' * vb2.drhovect(i, :)') + vb2.vECI_target(i, :);
end

delV_mission = delV_mission + vb2.DV_total/1000;
time_mission = time_mission + t_hold;

% ECI plot
figure()
earth_sphere(gca)
hold on
plot3(vb2.rECI_target(:,1), vb2.rECI_target(:,2), vb2.rECI_target(:,3), "--g") % full orbit bc 24 hour 
plot3(vb2.rECI_target(end,1), vb2.rECI_target(end,2), vb2.rECI_target(end,3), "diamond", "Color", "b") % new loc of target
plot3(vb2.rECI_chaser(end, 1), vb2.rECI_chaser(end, 2), vb2.rECI_chaser(end, 3), "square", "Color", "k")
xlim padded
ylim padded
zlim padded
xlabel("Downrange [km]")
ylabel("Altitude [km]")
title("ECI Frame: After V-bar Hold 2")
legend("", "Target path during V-bar hold", "Target position", "Chaser position", "Location", "best")
grid on

% LVLH plot
figure()
plot(0, 0, "diamond", "Linewidth", 2, "Color", "b")
hold on
plot(vb2.rhovect(:, 2)*1000, vb2.rhovect(:, 1)*1000, "--r", "LineWidth", 2)
plot(vb2.rhovect(end, 2)*1000, vb2.rhovect(end,1)*1000, "o", "Linewidth", 2, "Color", "k")
ylim padded 
xlim padded 
xlabel("Downrange [m]")
ylabel("Altitude [m]")
title("LVLH Frame: After V-bar hold 2")
legend("Target location", "V-bar hold trajectory", "Chaser Final location", "Location", "best")
grid on

disp("~~~~~~~~~~~~~~~~~~~~ V-bar Hold 2 Summary ~~~~~~~~~~~~~~~~~~~~")
disp("The chaser is " + norm(vb2.rhovect(end, 1:3))*1000 + " m away!")
disp("Total Delta-V used: " + delV_mission*1000 + " m/s.")
disp("Time elapsed: " + time_mission + " hours.")
disp(" ")

%% eric code
disp("~~~~~~~~~~~~~~~~~~~~ V-bar Approach ~~~~~~~~~~~~~~~~~~~~")

% Define initial distance of chaser and target final distance
initial_distance = 0.3; % initial distance in km (300 m in front of the target)
target_distance = 0.02; % target distance in km (20 m in front of the target)
% approach_speed = -0.00001; % nominal approach speed in km/s (negative to approach target)

% Extract final states from the last hold
r_final_chaser = vb2.rhovect(end, 1:3)'; % Chaser's relative position vector
v_final_chaser = vb2.drhovect(end, 1:3)'; % Chaser's relative velocity vector
t_approach = 24*60*60;
approach_speed = abs((target_distance - r_final_chaser(1)) / t_approach); % time to reach 20m in front of the target
% Set the initial distance in the LVLH frame
r_final_chaser(1) = initial_distance; % Chaser starts 300 m in front of the target
v_final_chaser(2:3) = 0; % No motion in radial (altitude) or out-of-plane direction
v_final_chaser(1) = approach_speed; % Set motion along the V-bar direction (downrange)

% Calculate the required acceleration to maintain a nominal approach
n_target = (2 * pi) / period;
a_approach = [-2 * n_target * approach_speed; 0; 0]; % Acceleration for V-bar approach

% Update the chaser's initial state for the approach maneuver
% approach_speed = abs((target_distance - r_final_chaser(1)) / t_approach); % time to reach 20m in front of the target
% t_approach = 5;

% Set initial state for the approach
state_approach = [r_final_chaser; v_final_chaser; vb2.rECI_target(end, 1:3)'; vb2.vECI_target(end, 1:3)'];

% Propagate the approach using ODE
tspan_approach = [0 t_approach];
[time_vbar, state_vbar] = ode45(@LinearEOM, tspan_approach, state_approach, options, mu_earth);

% Extract results from ODE propagation
vbar.rhovect = state_vbar(:, 1:3);
vbar.drhovect = state_vbar(:, 4:6);
vbar.rECI_target = state_vbar(:, 7:9);
vbar.vECI_target = state_vbar(:, 10:12);
vbar.rECI_chaser = vbar.rECI_target - vbar.rhovect;

% Plot ECI frame results
figure()
earth_sphere(gca)
hold on
plot3(vbar.rECI_target(:, 1), vbar.rECI_target(:, 2), vbar.rECI_target(:, 3), "--g")
plot3(vbar.rECI_target(end, 1), vbar.rECI_target(end, 2), vbar.rECI_target(end, 3), "diamond", "Color", "b")
plot3(vbar.rECI_chaser(end, 1), vbar.rECI_chaser(end, 2), vbar.rECI_chaser(end, 3), "square", "Color", "k")
xlim padded
ylim padded
zlim padded
xlabel("Downrange [km]")
ylabel("Altitude [km]")
title("ECI Frame: After V-bar Approach")
legend("", "Target path", "Target position", "Chaser position", "Location", "best")
grid on

% Plot LVLH frame results 
figure()
plot(0, 0, "diamond", "LineWidth", 2, "Color", "b") % Target stays at 0
hold on
plot(vbar.rhovect(:, 1) * 1000, zeros(size(vbar.rhovect(:, 1))), "--r", "LineWidth", 2) % Only plot x (downrange)
plot(vbar.rhovect(end, 1) * 1000, 0, "o", "LineWidth", 2, "Color", "k") % Chaser final location
xlim padded
ylim([-5 5]) % Set tight limits to emphasize downrange movement only
xlabel("Downrange [m]")
ylabel("Altitude [m]")
title("LVLH Frame: After V-bar Approach")
legend("Target location", "V-bar Approach Trajectory", "Chaser Final location", "Location", "best")
grid on

% Display approach summary
disp("The chaser started " + (initial_distance * 1000) + " m in front of the target.")
disp("The chaser is now " + round(norm(vbar.rhovect(end, 1:3)) * 1000, 2) + " m away from the target!")
% disp("Total Delta-V used during V-bar Approach: " + delV_mission * 1000 + " m/s.")
disp("Time elapsed for V-bar Approach: " + t_approach / 3600 + " hours.")
disp("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

%% all maneuver plot
% figure()
% plot(0, 0, "diamond", "LineWidth", 2, "Color", "m")
% hold on
% plot(vbar0, rbar0, "o", "Linewidth", 2, "Color", "b")
% plot(hop1.rhovect(:,2), hop1.rhovect(:,1), "--k", "LineWidth", 2)
% plot(fb1.rhovect(:, 2), fb1.rhovect(:, 1), "--k", "LineWidth", 2)
% plot(hop2.rhovect(:, 2), hop2.rhovect(:, 1), "--k", "LineWidth", 2)
% plot(vb1.rhovect(:, 2), vb1.rhovect(:, 1), "--k", "LineWidth", 2)
% plot(hop3.rhovect(:, 2), hop3.rhovect(:, 1), "--k", "LineWidth", 2)
% plot(vb2.rhovect(:, 2), vb2.rhovect(:, 1), "--k", "LineWidth", 2)
% 
% ylim padded 
% xlim padded 
% xlabel("Downrange [km]")
% ylabel("Altitude [km]")
% title("LVLH Frame: Relative Distance at Start of Mission")
% legend("Target", "Chaser", "Location", "best")
% grid on

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

    
    i = rA / norm(rA); 
    k = hA / norm(hA); 
    j = cross(k,i);
    
    DCM = [i'; j'; k'];
    
    Omega = hA / (norm(rA)^2);
    Omega_dot = -2*dot(rA,vA) / norm(rA)^2 * Omega;

    aA = -mu * rA / norm(rA)^3;
    aB = -mu * rB / norm(rB)^3;

    r_rel = rB - rA;
    v_rel = vB - vA - cross(Omega,r_rel);

    a_rel = aB - aA - cross(Omega_dot,r_rel) - cross(Omega,cross(Omega,r_rel))- 2 * cross(Omega,v_rel);

    rhovect = DCM*r_rel;
    drhovect = DCM*v_rel;
    ddrhovect = DCM*a_rel;
end

function [dv0_plus,DV_0,DV_f,DV_total] = cw_rv(dr,drf,dv0,period,t)
% CIRCULAR ORBIS ONLY
n = 2*pi / period;
% Matrix Method
phiRR = phi_rr(n, t);
phiRV = phi_rv(n, t);
phiVR = phi_vr(n, t);
phiVV = phi_vv(n, t);

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

function rr_matrix = phi_rr(n, t)
% t - time in seconds
% n - mean motion in ???
    rr_matrix = [4-3*cos(n*t) 0 0; 
                 6*(sin(n*t)-n*t) 1 0;
                 0 0 cos(n*t)];
end

function rv_matrix = phi_rv(n, t)
% t - time in seconds
% n - mean motion in ???
    rv_matrix = [1/n*sin(n*t) 2/n*(1-cos(n*t)) 0; 
                 2/n*(cos(n*t)-1) 1/n*(4*sin(n*t) - 3*n*t) 0;
                 0 0 1/n*sin(n*t)];
end

function vr_matrix = phi_vr(n, t)
% t - time in seconds
% n - mean motion in ???
    vr_matrix = [3*n*sin(n*t) 0 0;
                 6*n*(cos(n*t) -1) 0 0;
                 0 0 -n*sin(n*t)];
end

function vv_matrix = phi_vv(n, t)
% t - time in seconds
% n - mean motion in ???
    vv_matrix = [cos(n*t) 2*sin(n*t) 0;
                 -2*sin(n*t) (4*cos(n*t) -3) 0
                 0 0 cos(n*t)];
end

function dydt = circularLinearEOM(time,state,n,mu)
y = state; 
rvect = state(7:9);
R = norm(rvect);

% chaser
dydt(1) = y(4);
dydt(2) = y(5);
dydt(3) = y(6);

dydt(4) = 3*n^2*y(1) + 2*n*y(5) - 2*n*y(5);
dydt(5) = -2*n*y(4); 
dydt(6) = -n^2*y(3);

% target
dydt(7) = y(10);
dydt(8) = y(11);
dydt(9) = y(12);
dydt(10) = -mu * y(7) / R^3;
dydt(11) = -mu * y(8) / R^3;
dydt(12) = -mu * y(9) / R^3;

dydt = dydt';
end
