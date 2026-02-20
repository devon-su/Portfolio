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