function [r, v] = COEs2rv(h, ecc, theta, RAAN, inc, AoP, a, mu)
% coes in rad
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