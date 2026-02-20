function [v, a_minE, t_minE] = lambert_minE (r1, r2, dm, nrev)
        mu_earth = 398600;

        cos_dnu = dot(r1,r2)/(norm(r1)*norm(r2));
             
        c = sqrt(norm(r1)^2 + norm(r2)^2 - 2*norm(r1)*norm(r2)*cos_dnu);
        s = 1/2 * (norm(r1) + norm(r2) + c);
        a_minE = 0.5*s;

        alphae = pi;
        betae = 2*asin(sqrt((s-c)/s));
        
        if dm == "long"
            t_minE = sqrt(a_minE^3/mu_earth)*(2*nrev*pi + alphae - (betae-sin(betae)));
        else
            t_minE = sqrt(a_minE^3/mu_earth)*(2*nrev*pi + alphae + (betae-sin(betae)));
        end

        rcrossr = cross(r1,r2);
        pmin = norm(r1)*norm(r2)/c*(1 - cos_dnu);

        if dm == "short"
            sin_dnu = norm(rcrossr)/(norm(r1)*norm(r2));
        else
            sin_dnu = -norm(rcrossr)/(norm(r1)*norm(r2));            
        end
        for i = 1:3
            v(i) = sqrt(mu_earth*pmin)/(norm(r1)*norm(r2)*sin_dnu)*(r2(i) - (1-norm(r2)/pmin*(1-cos_dnu))*r1(i));
        end
end         


