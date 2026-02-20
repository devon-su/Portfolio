function rot_matrix = C_sez2ecef(a, b)
% Rotation matrix to convert from SEZ to ECEF

    rot_matrix = [sin(a)*cos(b) -sin(b) cos(a)*cos(b); 
              sin(a)*sin(b) cos(b) cos(a)*sin(b); 
              -cos(a) 0 sin(a)];
end

