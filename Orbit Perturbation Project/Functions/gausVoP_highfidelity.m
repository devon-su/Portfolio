function [dstate] = gausVoP(time, state, Cd, area, mass, drag, j2)
w_earth = [0; 0; 72.9211e-6];
mu_earth = 398600;
r_earth = 6378;

h = state(1);
ecc = state(2);
theta = state(3);
raan = state(4);
inc = state(5);
omega = state(6);

% From COES into perifocal frame
[rvect, vvect] = COEs2rv(h, ecc, theta, raan, inc, omega, 'idk', mu_earth);
r = norm(rvect);
    
% Build Rotation Matrix
R_hat = rvect/r;
N_hat = cross(rvect,vvect)/norm(cross(rvect,vvect));
T_hat = cross(N_hat,R_hat);
DCM = [R_hat T_hat N_hat];
DCM = DCM';
    
% drag pert
if drag == 1
    t0 = datetime("2024-11-09 00:00:00");
    JD0 = juliandate(t0);
    JD = JD0 + time/(24*3600);

    utc = datetime(JD, 'convertfrom', 'juliandate');
    [year,month,d] = ymd(utc);
    dayofyear = day(utc,"dayofyear");
    [hour,min,sec] = hms(utc);
    UTsec = hour*3600 + min*60 + sec;
    utc_long = [year month d hour min sec];

    lla = eci2lla(rvect'*10^3, utc_long);
    lat = lla(1);
    long = lla(2);
    alt2 = abs(lla(3));

    [~, rho_MSISE] = atmosnrlmsise00(alt2, lat, long, year, dayofyear, UTsec);
    
    rho = rho_MSISE(6)*1000^3; %kg/km^3

    Vrel = vvect - cross(w_earth,rvect);
    acc_drag = -0.5*Cd*area*(1/mass).*rho.*(norm(Vrel)).*(Vrel);
    acc_drag = DCM*acc_drag;
else
    acc_drag = 0;
end

% j2 pert 
if j2 > 1
    acc_J2 = J2_total(rvect, DCM, mu_earth, r_earth, j2);
else
    acc_J2 = 0;
end

p = acc_drag + acc_J2; % + acc_SRP;
R = p(1);
T = p(2);
N = p(3);
    
dh = r*T;
    
decc = (h/mu_earth)*sin(theta)*R + (1/(mu_earth*h))*((h^2 + mu_earth*r)*cos(theta) + mu_earth*ecc*r)*T;
    
acc_2b = h/r^2;
dtheta_pert = (1/(ecc*h))*((h^2/mu_earth)*cos(theta)*R - (r + (h^2/mu_earth))*sin(theta)*T);
dtheta =  acc_2b + dtheta_pert;
draan = (r/(h*sin(inc)))*sin(omega + theta)*N;
dinc = (r/h)*cos(omega + theta)*N;
domega = -dtheta_pert- ((r*sin(omega + theta))/(h*tan(inc)))*N;
    
dstate = [dh; decc; dtheta; draan; dinc; domega];
end