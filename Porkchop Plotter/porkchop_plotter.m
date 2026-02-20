function porkchop_plotter(depart_date1, depart_date2, arrival_date1, arrival_date2, depart_planet, arrive_planet)

% Inputs
% depart_date - departure date(s)
% arrival_date - arrival date(s)
% depart_planet and arrive_planet - planet identifier for departure and arrival respectively
% 1 = Mercury
% 2 = Venus
% 3 = Earth
% 4 = Mars
% 5 = Jupiter
% 7 = Uranus
% 8 = Neptune
% 9 = Pluto

% Outputs
% dV_short - delta V for short way lamberts [km/s]
% dV_long - delta V for long way lamberts [km/s]
% TOF - time of flight [days]
% V_inf_short 
% V_inf_long
% C3_short
% C3_long

mu_sun = 1.32712440018e11; % km^3/s^2

v_inf_cutoff = 100;
c3_cutoff = 10000;

% 1 day step size, using Julian date for simplicity
dep_array = juliandate(depart_date1:days(1):depart_date2);
arr_array = juliandate(arrival_date1:days(1):arrival_date2);

% preallocate
C3_short = zeros(length(arr_array), length(dep_array));
C3_long = zeros(length(arr_array), length(dep_array));

V_inf_short = zeros(length(arr_array), length(dep_array));
V_inf_long = zeros(length(arr_array), length(dep_array));

TOF = zeros(length(arr_array), length(dep_array));

dV_short = zeros(length(arr_array), length(dep_array));
dV_long = zeros(length(arr_array), length(dep_array));

% lamberts nested for loop
for i = 1:length(arr_array)
    for j = 1:length(dep_array)
    tof = arr_array(i) - dep_array(j);

    % earth and mars pos from keplers
    [~, r_dep, v_dep, ~] = AERO557planetcoe_and_sv(depart_planet, dep_array(j));
    [~, r_arr, v_arr, ~] = AERO557planetcoe_and_sv(arrive_planet, arr_array(i));

    % lamberts -- short way
    [v_dep_short, v_arr_short, ~, flag_short] = izzo_lambert(r_dep', r_arr', tof, 0, mu_sun);
    
    if flag_short < 0
        v_dep_short = 1e10*ones(3,1);
        v_arr_short = 1e10*ones(3,1);
    end

    c3_short = norm(v_dep_short - v_dep')^2;
    v_inf_short = norm(v_arr_short - v_arr');
        
    % long way
    [v_dep_long, v_arr_long, ~, flag_long] = izzo_lambert(r_dep', r_arr', -tof, 0, mu_sun);
    
    if flag_long < 0
        v_dep_long = 1e10*ones(3,1);
        v_arr_long = 1e10*ones(3,1);
    end

    c3_long = norm(v_dep_long - v_dep')^2;
    v_inf_long = norm(v_arr_long - v_arr');
    
    % cutoff
    if c3_short > c3_cutoff
        c3_short = NaN;
    end

    if c3_long > c3_cutoff
        c3_long = NaN;
    end

    if v_inf_short > v_inf_cutoff
        v_inf_short = NaN;
    end

    if v_inf_long > v_inf_cutoff
        v_inf_long = NaN;
    end

    % index everything
    C3_short(i, j) = c3_short;
    C3_long(i, j) = c3_long;

    V_inf_short(i, j) = v_inf_short;
    V_inf_long(i, j) = v_inf_long;

    TOF(i, j) = tof;
    
    % delta V
    dV_short(i, j) = v_inf_short + sqrt(c3_short);
    dV_long(i, j) = v_inf_long + sqrt(c3_long);
    end
end

departure = dep_array - juliandate(depart_date1);
arrival = arr_array - juliandate(arrival_date1);

% Plot -- Configured for Earth to Mars 
figure(1)
% c3
contour(departure, arrival, C3_short, 'ShowText', 'on', 'color', 'r')
hold on
contour(departure, arrival, C3_long, 'ShowText', 'on', 'color', 'r')
% v_inf
contour(departure, arrival, V_inf_short, 0:5:40, 'ShowText', 'on', 'color', 'b')
contour(departure, arrival, V_inf_long, 0:5:40, 'ShowText', 'on', 'color', 'b')
% time of flight
contour(departure, arrival, TOF, 0:60:600, 'ShowText', 'on', "color", "#808080")

grid on
xlabel("Departure: Days past " + string(depart_date1))
ylabel("Arrival: Days past " + string(arrival_date1))
legend("C3 [km2/s2]", "", "V_\infty at Mars [km/s]", "", "Time of flight [days]", "Location", "northwest")
title('Porkchop Plot')

figure(2)
% delta V
contour(departure, arrival, dV_short, 0:1:30, 'ShowText', 'on', 'color', 'b')
hold on
contour(departure, arrival, dV_long, 0:1:30, 'ShowText', 'on', 'color', 'r')
% time of flight
contour(departure, arrival, TOF, 60:60:480, 'ShowText', 'on', "color", "#808080")

grid on
xlabel("Departure: Days past " + string(depart_date1))
ylabel("Arrival: Days past " + string(arrival_date1))
legend("Type I \DeltaV [km/s]", "Type II \DeltaV [km/s]", "Time of flight [days]", "Location", "northwest")
title('\DeltaV')
end

