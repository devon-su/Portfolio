%% Aero 452 Project 2 -- Group 6
clc; clear; close all
addpath("Functions\")
addpath("Project 2")

global mu_earth r_earth
mu_earth = 398600;
r_earth = 6378;
cd = 2.2;
cr = 1.2;
w_earth = [0 0 72.9211e-6]; % [rad/s]
muSun = 132.712e9; % km3/s2
Psr = 4.57*10^-6;

%% LEO analysis
LEO = TLE_init("Tiros1TLE.txt", mu_earth);

LEO.mass = 122.5; % [kg]
LEO.area = pi*(107/100/1000/2)^2; % [km2]

% propogate to 1/1/2025 00:00:00
jd = juliandate(2025, 1, 1, 0, 0, 0);
LEO.jd = juliandate(LEO.epoch);

time_till = jd - LEO.jd;

COEs_state = [LEO.h; LEO.ecc; LEO.theta; LEO.RAAN; LEO.inc; LEO.omega];

[r_begin, v_begin] = COEs2rv(LEO.h, LEO.ecc, LEO.theta, LEO.RAAN, LEO.inc, LEO.omega, "idk", mu_earth);

options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
[~, state_begin] = ode45(@airplane, [0 time_till], [r_begin, v_begin], options, mu_earth);

%
t = 3*24*60*60; % [days] -> [sec]
t_span = [0 t];

% vop.options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'Events', @stop);
% gausVoP(time, state, Cd, area, mass, drag, j2, Cr, jd, SRP, nbody)
% [dh; decc; dtheta; draan; dinc; domega];
% [vop.time, vop.state] = ode45(@gausVoP, t_span, COEs_state, vop.options, cd, LEO.area, LEO.mass, 1, 6, cr, LEO.jd, 1, 1); 
vop.options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12,'Events', @eventDeOrbit);
[vop.time, vop.state] = ode45(@vop_ODE, t_span, COEs_state, vop.options, w_earth, r_earth, mu_earth, muSun, cd, LEO.area, LEO.mass, LEO.jd, cr, Psr);

for i = 1:length(vop.time)
    [r_temp, v_temp] = COEs2rv(vop.state(i, 1), vop.state(i, 2), vop.state(i, 3), vop.state(i, 4), vop.state(i, 5), vop.state(i, 6), "idk", mu_earth);
    r(i, 1:3) = r_temp;
    v(i, 1:3) = v_temp;
    h = vop.state(i, 1);
    ecc = vop.state(i, 2);
    a = (h^2)/(mu_earth*(1-ecc^2));
    ra(i) = a + a*ecc;
    rp(i) = 2*a - ra(i);
end

time = vop.time/24/60/60;

% trajectory of LEO object
figure()
earth_sphere(gca)
hold on
plot3(r(:,1), r(:,2), r(:,3), '.')
plot3(r(1,1), r(1,2), r(1,3), '*', 'LineWidth', 5)
plot3(r(end,1),r(end,2),r(end,3), '*', 'LineWidth', 5)
legend("", "Orbital Path", "Start Position", "End Position", 'Location', 'southoutside');
xlabel("X [Km]")
ylabel("Y [Km]")
zlabel("Z [Km]")

% r_p and r_a
rr_a = ra - r_earth;
rr_p = rp - r_earth;

figure
plot(time, rr_a, 'LineWidth', 2)
hold on
plot(time, rr_p, 'LineWidth', 2)
ylim padded 
xlim tight 
xlabel('Time [days]'); 
ylabel('Altitude [km]'); 
grid on 
legend('Apogee', 'Perigee', 'Location', 'best')

% COEs
figure()
tiledlayout("vertical")
nexttile
plot(time, vop.state(:,1) - LEO.h)
ylim padded 
xlim tight  
ylabel('h - h_0')
grid on 

nexttile
plot(time, vop.state(:,2) - LEO.ecc)
ylim padded 
xlim tight  
ylabel('ecc - ecc_0');  
grid on 

nexttile
plot(time, rad2deg(vop.state(:,4) - LEO.RAAN))
ylim padded 
xlim tight 
ylabel('\Omega - \Omega_0'); 
grid on 

nexttile
plot(time, rad2deg(vop.state(:,5) - LEO.inc)) 
ylim padded 
xlim tight 
ylabel('inc - inc_0');  
grid on 

nexttile
plot(time, rad2deg(vop.state(:,6) - LEO.omega))
ylim padded 
xlim tight 
xLab = xlabel('Time [days]'); 
yLab = ylabel('\omega - \omega_0'); 
grid on 

%% GEO analysis
GEO = TLE_init("intellsat3TLE.txt", mu_earth);

GEO.mass = 150; % kg
GEO.area = pi*(142/100/1000/2)^2; % [km2]

GEO.jd = juliandate(GEO.epoch);

t = 100*24*60*60; % [days] -> [sec]
t_span = [0 t];

COEs_state = [GEO.h; GEO.ecc; GEO.theta; GEO.RAAN; GEO.inc; GEO.omega];
[r_begin, v_begin] = COEs2rv(GEO.h, GEO.ecc, GEO.theta, GEO.RAAN, GEO.inc, GEO.omega, "idk", mu_earth);

% vop.options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'Events', @stop);
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
vop.options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12,'Events', @eventDeOrbit);
% gausVoP(time, state, Cd, area, mass, drag, j2, Cr, jd, SRP, nbody)
% [dh; decc; dtheta; draan; dinc; domega];
% [vop.time, vop.state] = ode45(@gausVoP, t_span, COEs_state, vop.options, cd, GEO.area, GEO.mass, 1, 6, cr, jd, 0, 1); 

r0 = 1e4*[3.12335, -2.8766, 0.0194];
v0 = [2.0738, 2.256, 0.0193];

[h, ecc, theta, RAAN, inc, AoP, a, T] = COEs(r0, v0, mu_earth);

COEs_state_brandon = [h, ecc, theta, RAAN, inc, AoP];

[vop.time, vop.state] = ode45(@vop_ODE, t_span, COEs_state, vop.options, w_earth, r_earth, mu_earth, muSun, cd, GEO.area, GEO.mass, GEO.jd, cr, Psr);
% [vop.time, vop.state] = ode45(@vop_ODE, t_span, COEs_state_brandon, vop.options, w_earth, r_earth, mu_earth, muSun, cd, GEO.area, GEO.mass, GEO.jd, cr, Psr);

for i = 1:length(vop.time)
    [r_temp, v_temp] = COEs2rv(vop.state(i, 1), vop.state(i, 2), vop.state(i, 3), vop.state(i, 4), vop.state(i, 5), vop.state(i, 6), "idk", mu_earth);
    r(i, 1:3) = r_temp;
    v(i, 1:3) = v_temp;
    h = vop.state(i, 1);
    ecc = vop.state(i, 2);
    a = (h^2)/(mu_earth*(1-ecc^2));
    ra(i) = a + a*ecc;
    rp(i) = 2*a - ra(i);
end

time = vop.time/24/60/60;

% trajectory of GEO object
figure()
earth_sphere(gca)
hold on
plot3(r(:,1), r(:,2), r(:,3), '.')
plot3(r(1,1), r(1,2), r(1,3), '*', 'LineWidth', 5)
plot3(r(end,1),r(end,2),r(end,3), '*', 'LineWidth', 5)
legend("", "Orbital Path", "Start Position", "End Position", 'Location', 'southoutside');
xlabel("X [Km]")
ylabel("Y [Km]")
zlabel("Z [Km]")

% r_p and r_a
rr_a = ra - r_earth;
rr_p = rp - r_earth;

figure
plot(time, rr_a, 'LineWidth', 2)
hold on
plot(time, rr_p, 'LineWidth', 2)
ylim padded 
xlim tight 
xlabel('Time [days]'); 
ylabel('Altitude [km]'); 
grid on 
legend('Apogee', 'Perigee', 'Location', 'best')

% COEs
figure()
tiledlayout("vertical")
nexttile
plot(time, vop.state(:,1) - GEO.h)
ylim padded 
xlim tight  
title('h - h_0')
grid on 

nexttile
plot(time, vop.state(:,2) - GEO.ecc)
ylim padded 
xlim tight  
title('ecc - ecc_0');  
grid on 

nexttile
plot(time, rad2deg(vop.state(:,4) - GEO.RAAN))
ylim padded 
xlim tight 
title('\Omega - \Omega_0'); 
grid on 

nexttile
plot(time, rad2deg(vop.state(:,5) - GEO.inc)) 
ylim padded 
xlim tight 
title('inc - inc_0');  
grid on 

nexttile
plot(time, rad2deg(vop.state(:,6) - GEO.omega))
ylim padded 
xlim tight 
xlabel('Time [days]'); 
title('\omega - \omega_0'); 
grid on 

%% Lamberts
[r_pert, v_pert] = COEs2rv(vop.state(end, 1), vop.state(end, 2), vop.state(end, 3), vop.state(end, 4), vop.state(end, 5), vop.state(end, 6), "idk", mu_earth);

dt = 24*60*60; 
tm = -1;

[GEO.rv, GEO.vv] = COEs2rv(GEO.h, GEO.ecc, GEO.theta, GEO.RAAN, GEO.inc, GEO.omega, "idk", mu_earth);
GEO.state = [GEO.rv, GEO.vv];

[~, state_coast] = ode45(@airplane, t_span, GEO.state, options, mu_earth);

rf = [state_coast(end,1), state_coast(end,2), state_coast(end,3)];
vf = [state_coast(end,4), state_coast(end,5), state_coast(end,6)];

[v1vec, v2vec] = lambert(r_pert, rf, dt, tm, mu_earth);
deltaV = norm(v1vec - v_pert) + norm(v2vec - vf);

state_coasting = [r_pert, v_pert];
[~, state_lambert] = ode45(@airplane, [0 dt], state_coasting, options, mu_earth);

figure()
earth_sphere(gca)
hold on
plot3(r(:,1), r(:,2), r(:,3), 'LineWidth', 1)
plot3(state_coast(:,1),state_coast(:,2),state_coast(:,3),'LineWidth',1)
plot3(state_lambert(:,1),state_lambert(:,2),state_lambert(:,3),'LineWidth',4)
plot3(state_lambert(1,1),state_lambert(1,2),state_lambert(1,3),'*','LineWidth',5)
plot3(state_lambert(end,1),state_lambert(end,2),state_lambert(end,3),'*','LineWidth',5)
legend("","Perturbed orbit","Osculating Orbit","Lambert Trajectory","Initial Position","Final Position",'Location','southoutside');

%%