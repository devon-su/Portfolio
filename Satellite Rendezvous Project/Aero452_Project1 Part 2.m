% Constants
mu = 3.986e14; % m^3/s^2
R_earth = 6378e3; % m
h = 35786e3; % m
r = R_earth + h; % m
n = sqrt(mu / r^3); % rad/s

% Orbital data
Ecc = 0.0001055; % ecc
Inc = deg2rad(0.0461); % rad
RAAN = deg2rad(87.6205); %
ArgPerigee = deg2rad(187.0055); % rad
MeanAnomaly = deg2rad(217.8251); % rad

n = sqrt(mu / r^3); % rad/s
T = 2 * pi / n; % Orbital period in seconds

[TargetR, TargetV] = orbitalElementsToECI(mu, r, Ecc, Inc, RAAN, ArgPerigee, MeanAnomaly);

% Chaser starts 1 km behind the target in the same circular orbit
InitialRelDist = 1000; % 1 km behind in the velocity direction
FinalRelDist = 300; % Final distance of 300 meters ( SETUP for next )
ChaserR = TargetR - (TargetV / norm(TargetV)) * InitialRelDist; % 1km behind the target
ChaserV = TargetV; % given

% Nominal approach speed 
Vc = 0.01; % assumption

% Simulation Parameters
dt = 60; % Time step (seconds)
TotTime = 10 * 24 * 3600; % 10 days in seconds, but we stop earlier
time = 0:dt:TotTime;

% Initialize arrays for position, velocity, relative distance and LVLH frame
ChaserPos = zeros(length(time), 3);
TargetPos = zeros(length(time), 3);
RelDist = zeros(length(time), 1);
RelPosLVLH = zeros(length(time), 3); % To store relative position in LVLH frame

% Set initial conditions
ChaserPos(1, :) = ChaserR';
TargetPos(1, :) = TargetR';
RelDist(1) = InitialRelDist; % Start with 1 km relative distance

% propogate orbit for both chaser and target
for i = 2:length(time)
    % updating orbit
    ThetaTarget = MeanAnomaly + n * time(i);
    ThetaChaser = MeanAnomaly + n * time(i) + (Vc * time(i) / r); % Chaser gradually catching up
    
    % Update positions
    TargetPos(i, :) = [r * cos(ThetaTarget), r * sin(ThetaTarget), 0];
    ChaserPos(i, :) = [r * cos(ThetaChaser), r * sin(ThetaChaser), 0];
    
    % Relative distance between chaser and target
    RelDist(i) = norm(ChaserPos(i, :) - TargetPos(i, :));
    
    % ECI to LVLH
    RelPosECI = ChaserPos(i, :) - TargetPos(i, :);
    RelPosLVLH(i, :) = ECI2LVLH(TargetPos(i, :), TargetV, RelPosECI);
   
    % Stop adjusting once the relative distance is close to the desired final distance (300 meters)
    if RelDist(i) <= FinalRelDist
        Vc = 0; % Stop further velocity adjustments when within 300 meters
        RelDist(i:end) = FinalRelDist; % Maintain the relative distance at 300 meters
        break;
    end
    
    % Vc is nominal approach speed, delta-v applied each time step
    DeltaVy = -2 * n * Vc * dt; 
    Vc = Vc + DeltaVy; % Update velocity
    
    % Adjust chaser's velocity along the V-bar axis
    ChaserV = ChaserV + DeltaVy * (ChaserV / norm(ChaserV));
end

% Trimming unused values
ChaserPos = ChaserPos(1:i, :);
TargetPos = TargetPos(1:i, :);
RelDist = RelDist(1:i);
RelPosLVLH = RelPosLVLH(1:i, :);
time = time(1:i);

% Plot relative motion in ECI
figure
plot3(ChaserPos(:, 1) / 1000, ChaserPos(:, 2) / 1000, ChaserPos(:, 3) / 1000, 'r', 'DisplayName', 'Chaser Orbit in ECI')
hold on
plot3(TargetPos(:, 1) / 1000, TargetPos(:, 2) / 1000, TargetPos(:, 3) / 1000, 'b', 'DisplayName', 'Target Orbit in ECI')
title('Relative Motion in ECI Frame')
xlabel('X position (km)')
ylabel('Y position (km)')
zlabel('Z position (km)')
legend
grid on
axis equal

% Plot relative motion in LVLH frame (2D: Radial vs In-Track)
figure
plot(RelPosLVLH(:, 1) / 1000, RelPosLVLH(:, 2) / 1000, 'g', 'DisplayName', 'Chaser Motion in LVLH')
title('Relative Motion in LVLH Frame (Radial vs In-Track)')
xlabel('Radial (km)')
ylabel('In-Track (km)')
legend
grid on
axis equal

% Plot relative motion in LVLH frame (2D: Radial vs Cross-Track)
figure
plot(RelPosLVLH(:, 1) / 1000, RelPosLVLH(:, 3) / 1000, 'g', 'DisplayName', 'Chaser Motion in LVLH')
title('Relative Motion in LVLH Frame (Radial vs Cross-Track)')
xlabel('Radial (km)')
ylabel('Cross-Track (km)')
legend
grid on
axis equal

% Plot approach distance over time
figure;
plot(time / 3600, RelDist);
title('Approach Distance Over Time')
xlabel('Time (hours)')
ylabel('Relative Distance (m)')
grid on

final_distance = RelDist(end);
disp(['Final relative distance between chaser and target: ', num2str(final_distance), ' meters.']);

%% Function
function RelPosLVLH = ECI2LVLH(TargetPos, TargetVel, RelPosECI)
    % LVLH frame: X (radial), Y (along-track), Z (cross-track)
    % Inputs:
    % TargetPos: Target position in ECI (m)
    % TargetVel: Target velocity in ECI (m/s)
    % RelPosECI: Relative position of chaser to target in ECI (m)
    % Output:
    % RelPosLVLH: Relative position of chaser to target in LVLH frame (m)
    
    % Define unit vectors for LVLH frame
    Z_LVLH = -TargetPos / norm(TargetPos); % Radial direction (toward Earth)
    Y_LVLH = cross(TargetPos, TargetVel); % Cross-track direction
    Y_LVLH = Y_LVLH / norm(Y_LVLH); % Normalize
    X_LVLH = cross(Y_LVLH, Z_LVLH); % Along-track direction
    
    % Rotation matrix from ECI to LVLH
    R_ECI_to_LVLH = [X_LVLH; Y_LVLH; Z_LVLH];
    
    % Transform relative position to LVLH frame
    RelPosLVLH = R_ECI_to_LVLH * RelPosECI';
    RelPosLVLH = RelPosLVLH'; % Return as row vector
end

function [r_eci, v_eci] = orbitalElementsToECI(mu, a, Ecc, inc, RAAN, ArgPerigee, true_anomaly)
    % Converts orbital elements to ECI position and velocity vectors
    % Inputs:
    % mu - Gravitational parameter (m^3/s^2)
    % a - Semi-major axis (m)
    % Ecc - Eccentricity
    % inc - Inc (radians)
    % RAAN - Right Ascension of Ascending Node (radians)
    % ArgPerigee - Argument of perigee (radians)
    % true_anomaly - True anomaly at epoch (radians)
    %
    % Outputs:
    % r_eci - Position vector in ECI (m)
    % v_eci - Velocity vector in ECI (m/s)
    
    % Step 1: Calculate the distance and speed in the orbital plane
    p = a * (1 - Ecc^2); % Semi-latus rectum
    r_orbit = p / (1 + Ecc * cos(true_anomaly)); % Distance from focus to object
    v_orbit = sqrt(mu * (2 / r_orbit - 1 / a)); % Orbital velocity
    
    % Step 2: Position in Perifocal coordinates (PQW frame)
    r_pqw = [r_orbit * cos(true_anomaly); r_orbit * sin(true_anomaly); 0];
    v_pqw = [-sin(true_anomaly); Ecc + cos(true_anomaly); 0] * sqrt(mu / p);
    
    % Step 3: Rotation matrices to transform from PQW to ECI
    R3_W = [cos(-RAAN), sin(-RAAN), 0; -sin(-RAAN), cos(-RAAN), 0; 0, 0, 1]; % Rotation about z-axis (RAAN)
    R1_i = [1, 0, 0; 0, cos(-inc), sin(-inc); 0, -sin(-inc), cos(-inc)]; % Rotation about x-axis (Inc)
    R3_w = [cos(-ArgPerigee), sin(-ArgPerigee), 0; -sin(-ArgPerigee), cos(-ArgPerigee), 0; 0, 0, 1]; % Rotation about z-axis (arg of perigee)
    
    % Complete rotation from PQW to ECI
    Q_pqw2eci = (R3_W * R1_i * R3_w)';
    
    % Step 4: Transform position and velocity vectors to ECI
    r_eci = Q_pqw2eci * r_pqw;
    v_eci = Q_pqw2eci * v_pqw;
end

