function [el,beta] = r2elaz(rho,lat,lon,lst)
% Convert r vector into elevation and azimuth (deg)
% MAKE SURE YOUR rho = SLANT RANGE IN ECI, OR SET LST = LON IF IN ECEF
% VERIFIED

   rhoecef = R3(lst-lon)*rho';
    %rhoecef = R3(lst)*rho'
    rhosez = R2(90-lat)*R3(lon)*rhoecef
    el = asin(rhosez(3)/norm(rhosez)); % Elevation, rad
    el = mod(el,2*pi);
    if el ~= pi/2
       if atan2(rhosez(1),rhosez(2)) < pi/2
           % 0 to 180 deg for beta
           beta = acos(-rhosez(1)/sqrt(rhosez(1)^2 + rhosez(2)^2));
           beta = mod(beta,2*pi);
       else
           % 180 to 360 deg for beta
           beta = acos(-rhosez(1)/sqrt(rhosez(1)^2 + rhosez(2)^2));
           beta = 2*pi - mod(beta,2*pi);
       end
    else
        % error('Elevation is 90 degrees');
        beta = 0;
    end
    el = el*180/pi;
    beta = beta*180/pi;
    
end

