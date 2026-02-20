function v2 =  hgibbs (r1, r2, r3, t1, t2, t3)
% Herrick Gibbs orbit determination method using 3 vectors at 3 different times
% ------------------------------------------------------------------------------
% time in datetime format
% rvect in ECI
% All angles in RADS
% ------------------------------------------------------------------------------

% Gravitational Parameter
mu = 398600; % (km^3/s^2)

dt21 = seconds(t2-t1);
dt31 = seconds(t3-t1);
dt32 = seconds(t3-t2);

cross_r23 = cross(r2, r3);
a_cop =  asin(dot(cross_r23/norm(cross_r23), r1/norm(r1)));

% Check for coplanar tolerance - if not coplanar then Herrick Gibbs fails
if abs(a_cop) > deg2rad(1)
    warning('Not coplanar!')
end

theta12  = acos(dot(r1/norm(r1), r2/norm(r2)));
theta23 = acos(dot(r2/norm(r2), r3/norm(r3)));

% Checking for vector angles - if greater than 5 degrees calculations will not be 100% accurate
if theta12 > deg2rad(5) || theta23 > deg2rad(5)
    warning('Angle between vectors is greater than 5 degrees!')
end

term1 = -dt32*(1/(dt21*dt31) + mu/(12*norm(r1)*norm(r1)*norm(r1)));
term2 = (dt32-dt21)*(1.0/(dt21*dt32) + mu/(12*norm(r2)*norm(r2)*norm(r2)));
term3 =  dt21*(1/(dt32*dt31) + mu/(12*norm(r3)*norm(r3)*norm(r3)));

% Velocity Vector estimation using Taylor Series expansion
v2 =  term1*r1 + term2*r2 + term3*r3;
end