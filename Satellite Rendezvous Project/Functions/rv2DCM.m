function [Q_eci2lvlh] = rv2DCM(r, v)
% Creates a rotation matrix (ECI -> LVLH) given an r and v vector

    % CALCULATE ANGULAR MOMENTUM VECTOR OF S/C A
    h = cross(r, v); % km2/s

    % CALCULATE UNIT VECTORS OF LVLH FRAME
    i_hat = r / norm(r); 
    k_hat = h / norm(h); % mag of h never changes (assume no perturbations)
    j_hat = cross(k_hat, i_hat);
    
    % CALCULATE ORTHOGONAL DIRECTION COSINE MATRIX (QXx) using Eq. 7.11
    Q_eci2lvlh = [i_hat(1), i_hat(2), i_hat(3);
                  j_hat(1), j_hat(2), j_hat(3);
                  k_hat(1), k_hat(2), k_hat(3)];
end