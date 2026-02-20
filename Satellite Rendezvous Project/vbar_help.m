clc; clear; close all
mu_earth = 398600;

period = 8.6165e4;
t_hold = 24;

hop2.rhovect = [0; 1; 0]; % 1 km away from target
hop2.drhovect = 1.0e-3*[-0.526 0.0011 0];

hop2.rECI_chaser = 1.0e+4*[-0.0883   -4.2156   -0.0001];
hop2.vECI_chaser = [3.0740   -0.0644   -0.0025];

vb1.dr = hop2.rhovect;
vb1.pos = vb1.dr;
vb1.dv0 = hop2.drhovect;
vb1.dv = vb1.dv0;

t = t_hold*60*60;

[vb1.dv0_plus, vb1.DV_0, vb1.DV_f, vb1.DV_total] = cw_rv(vb1.dr, vb1.pos, vb1.dv0, period, t);

tspan = [0 t];
dv = vb1.dv0_plus;

vb1.state = [vb1.dr; dv; hop2.rECI_chaser'; hop2.vECI_chaser'];

options = odeset('RelTol', 1e-8, 'AbsTol',1e-8);
[vb1.timenew, vb1.statenew] = ode45(@circularLinearEOM, tspan, vb1.state, options, period, mu_earth);

vb1.rECI_target = [vb1.statenew(:,7), vb1.statenew(:,8), vb1.statenew(:,9)];
vb1.vECI_target = [vb1.statenew(:,10), vb1.statenew(:,11), vb1.statenew(:,12)];
vb1.rhovect = [vb1.statenew(:,1), vb1.statenew(:,2), vb1.statenew(:,3)];
vb1.drhovect = [vb1.statenew(:,4), vb1.statenew(:,5), vb1.statenew(:,6)];

vb1.rECI_chaser = vb1.rECI_target - vb1.rhovect;
for i = 1:length(vb1.timenew)
    vb1.DCM = rv2DCM(vb1.rECI_target(i, :)', vb1.vECI_target(i, :)');
    vb1.vECI_chaser = (vb1.DCM' * vb1.drhovect(i, :)') + vb1.vECI_target(i, :);
end

% delV_mission = delV_mission + vb1.DV_total/1000;
% time_mission = time_mission + t_hold;

figure()
plot(0, 0, "diamond", "Linewidth", 2, "Color", "b")
hold on
plot(vb1.rhovect(:, 2), vb1.rhovect(:, 1), "--r", "LineWidth", 2)
plot(vb1.rhovect(end, 2), vb1.rhovect(end,1), "o", "Linewidth", 2, "Color", "k")
ylim padded 
xlim padded 
xlabel("Downrange [km]")
ylabel("Altitude [km]")
title("LVLH Frame: After V-bar hold 1")
legend("Target location", "V-bar hold trajectory", "Chaser Final location", "Location", "best")
grid on

%%
t_hold = 24;

vb1.dr = hop2.rhovect(end, 1:3)';
vb1.pos = vb1.dr;
vb1.dv0 = hop2.drhovect(end, 1:3);
vb1.dv = vb1.dv0;

t = t_hold*60*60;

[vb1.dv0_plus, vb1.DV_0, vb1.DV_f, vb1.DV_total] = cw_rv(vb1.dr, vb1.pos, vb1.dv0, period, t);

tspan = [0 t];
dv = vb1.dv0_plus;

vb1.state = [vb1.dr; dv; hop2.rECI_chaser(end, 1:3)'; hop2.vECI_chaser(end, 1:3)'];
[vb1.timenew, vb1.statenew] = ode45(@circularLinearEOM, tspan, vb1.state, options, period, mu_earth);

vb1.rECI_target = [vb1.statenew(:,7), vb1.statenew(:,8), vb1.statenew(:,9)];
vb1.vECI_target = [vb1.statenew(:,10), vb1.statenew(:,11), vb1.statenew(:,12)];
vb1.rhovect = [vb1.statenew(:,1), vb1.statenew(:,2), vb1.statenew(:,3)];
vb1.drhovect = [vb1.statenew(:,4), vb1.statenew(:,5), vb1.statenew(:,6)];

vb1.rECI_chaser = vb1.rECI_target - vb1.rhovect;
for i = 1:length(vb1.timenew)
    vb1.DCM = rv2DCM(vb1.rECI_target(i, :)', vb1.vECI_target(i, :)');
    vb1.vECI_chaser = (vb1.DCM' * vb1.drhovect(i, :)') + vb1.vECI_target(i, :);
end

% LVLH plot
figure()
plot(0, 0, "diamond", "Linewidth", 2, "Color", "b")
hold on
plot(vb1.rhovect(:, 2), vb1.rhovect(:, 1), "--r", "LineWidth", 2)
plot(vb1.rhovect(end, 2), vb1.rhovect(end,1), "o", "Linewidth", 2, "Color", "k")
ylim padded 
xlim padded 
xlabel("Downrange [km]")
ylabel("Altitude [km]")
title("LVLH Frame: After V-bar hold 1")
legend("Target location", "V-bar hold trajectory", "Chaser Final location", "Location", "best")
grid on
