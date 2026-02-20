function acc_J2total = J2_total(rvect, DCM, mu_earth, r_earth, j2)
if j2 == 6
    acc_J2 = J2(rvect, DCM, mu_earth, r_earth);
    acc_J3 = J3(rvect, DCM, mu_earth, r_earth);
    acc_J4 = J4(rvect, DCM, mu_earth, r_earth);
    acc_J5 = J5(rvect, DCM, mu_earth, r_earth);
    acc_J6 = J6(rvect, DCM, mu_earth, r_earth);
    acc_J2total = acc_J2 + acc_J3 + acc_J4 + acc_J5 + acc_J6;
elseif j2 == 5
    acc_J2 = J2(rvect, DCM, mu_earth, r_earth);
    acc_J3 = J3(rvect, DCM, mu_earth, r_earth);
    acc_J4 = J4(rvect, DCM, mu_earth, r_earth);
    acc_J5 = J5(rvect, DCM, mu_earth, r_earth);
    acc_J2total = acc_J2 + acc_J3 + acc_J4 + acc_J5;
elseif j2 == 4
    acc_J2 = J2(rvect, DCM, mu_earth, r_earth);
    acc_J3 = J3(rvect, DCM, mu_earth, r_earth);
    acc_J4 = J4(rvect, DCM, mu_earth, r_earth);
    acc_J2total = acc_J2 + acc_J3 + acc_J4;
elseif j2 == 3
    acc_J2 = J2(rvect, DCM, mu_earth, r_earth);
    acc_J3 = J3(rvect, DCM, mu_earth, r_earth);
    acc_J2total = acc_J2 + acc_J3;
elseif j2 == 2
    acc_J2 = J2(rvect, DCM, mu_earth, r_earth);
    acc_J2total = acc_J2;
end
end
