function [r2vec, v2vec] = AERO557doubleR(qhat1, qhat2, qhat3, R1, R2, R3, tau1, tau3, r_guess)
% R1, R2, R3 - sites ECI
% qhat1, qhat2, qhat3 - L from RA and dec
% r_guess from Gauss

        % initial guess for the iteration 
        r1 = r_guess;
        r2 = 1.01*r1;
        c1 = dot(2*qhat1,R1);
        c2 = dot(2*qhat2,R2);
        count = 1;
        Err = 1;
        while Err > .001
            [F1,F2,f,g,r3vec,r2vec]=AERO557drfunc(c1,c2,r1,r2,tau1,tau3,qhat1,qhat2,qhat3,R1,R2,R3);
            dr1 = .005*r1;
            dr2 = .005*r2;
            [F1r1dr1,F2r1dr1,~,~]=AERO557drfunc(c1,c2,r1+dr1,r2,tau1,tau3,qhat1,qhat2,qhat3,R1,R2,R3);
            [F1r2dr2,F2r2dr2,~,~]=AERO557drfunc(c1,c2,r1,r2+dr2,tau1,tau3,qhat1,qhat2,qhat3,R1,R2,R3);
            
            dF1dr1 = (F1r1dr1 - F1)/dr1;
            dF2dr1 = (F2r1dr1 - F2)/dr1;
            dF1dr2 = (F1r2dr2 - F1)/dr2;
            dF2dr2 = (F2r2dr2 - F2)/dr2;
            
            %run the function with the dr conditions
            
            delta = dF1dr1*dF2dr2 - dF2dr1*dF1dr2;
            delta1 = dF2dr2*F1 - dF1dr2*F2;
            delta2 = dF1dr1*F2 - dF2dr1*F1;
            dr1 = -delta1/delta;
            dr2 = -delta2/delta;
            
            Err = (abs(dr1) + abs(dr2))/2;
            r1 = r1 + dr1;
            r2 = r2 + dr2;
            count=count+1;
            
        end
        
        %after convergence
        v2vec = (r3vec - f*r2vec)/g;
 
        %% double R function
    function [F1,F2,f,g,r3vec,r2vec]=AERO557drfunc(c1,c2,r1,r2,tau1,tau3,qhat1,qhat2,qhat3,R1,R2,R3)
        q1 = (-c1+sqrt(c1^2 - 4*(dot(R1,R1)-r1^2)))/2;
        q2 = (-c2+sqrt(c2^2 - 4*(dot(R2,R2)-r2^2)))/2;
        
        r1vec = R1 + q1*qhat1;
        r1 = norm(r1vec);
        r2vec = R2 + q2*qhat2;
        r2 = norm(r2vec);
        
        what = cross(r1vec,r2vec)/(norm(r1vec)*norm(r2vec));
        q3 = dot(-R3,what)/dot(qhat3,what);
        r3vec = R3 + q3*qhat3;
        r3 = norm(r3vec);
        
        %difference angles
        costheta21 = dot(r2vec,r1vec)/(norm(r2vec)*norm(r1vec));
        t21 = acos(costheta21);
        costheta31 = dot(r3vec,r1vec)/(norm(r3vec)*norm(r1vec));
        t31 = acos(costheta31);
        costheta32 = dot(r3vec,r2vec)/(norm(r3vec)*norm(r2vec));
        t32 = acos(costheta32);
        
        % if what(3)>0
        sintheta21 = sqrt(1-costheta21^2);
        sintheta31 = sqrt(1-costheta31^2);
        sintheta32 = sqrt(1-costheta32^2);
        % elseif what(3)<0
        %     sintheta21 = -sqrt(1-costheta21^2);
        %     sintheta31 = -sqrt(1-costheta31^2);
        %     sintheta32 = -sqrt(1-costheta32^2);
        % end
        
        theta31 = acos(costheta31);
        if theta31>pi
            c1 = (r2*sintheta32)/(r1*sintheta31);
            c3 = (r2*sintheta21)/(norm(r3vec)*sintheta31);
            p = (c1*r1 + c3*r3 - r2)/(c1 + c3 -1);
        elseif theta31<=pi
            c1 = (r1*sintheta31)/(r2*sintheta32);
            c3 = (r1*sintheta21)/(r3*sintheta32);
            p = (c3*r3 - c1*r2 + r1)/(-c1 + c3 +1);
        end
        
        ecostheta1 = p/r1 -1;
        ecostheta2 = p/r2 -1;
        ecostheta3 = p/r3 -1;
        
        theta21 = acos(costheta21);
        
        if theta21 ~= pi
            esintheta2 = (-costheta21*ecostheta2 + ecostheta1)/sintheta21;
        else
            esintheta2 = (costheta32*ecostheta2 - ecostheta3)/sintheta31;
        end
        
        e2 = ecostheta2^2 + esintheta2^2;
        a = p/(1-e2);
        
        if sqrt(e2)<1
            n = sqrt(398600/a^3);
            S = r2/p*sqrt(1-e2)*esintheta2;
            C = r2/p*(e2 + ecostheta2);
            sinE32 = r3/sqrt(a*p)*sintheta32 - r3/p*(1 - costheta32)*S;
            cosE32 = 1 - (r2*r3)/(a*p)*(1-costheta32);
            sinE21 = r1/sqrt(a*p)*sintheta21 + r1/p*(1 - costheta21)*S;
            cosE21 = 1 - (r2*r1)/(a*p)*(1-costheta21);
            M32 = acos(cosE32) + 2*S*sin(acos(cosE32)/2)^2 - C*sinE32;
            M12 = -acos(cosE21) + 2*S*sin(acos(cosE21)/2)^2 + C*sinE21;
            F1 = tau1 - M12/n;
            F2 = tau3 - M32/n;
            f = 1 - a/r2*(1 - cosE32);
            g = tau3 - sqrt(a^3/398600)*(acos(cosE32) - sinE32);
        else
            n = sqrt(398600/-a^3);
            Sh = r2/p*sqrt(e2-1)*esintheta2;
            Ch = r2/p*(e2 + ecostheta2);
            sinhF32 = r3/sqrt(-a*p)*sintheta32 - r3/p*(1 - costheta32)*Sh;
            F32 = log(sinhF32 + sqrt(sinhF32 + 1));
            sinhF21 = r1/sqrt(-a*p)*sintheta21 + r1/p*(1 - costheta32)*Sh;
            F21 = log(sinhF21 + sqrt(sinhF21^2 + 1));
            M32 = -F32 + 2*Sh*sinh(F32/2)^2 + Ch*sinhF32;
            M12 = F21 + 2*Sh*sinh(F21/2)^2 + Ch*sinhF21;
            F1 = tau1 - M12/n;
            F2 = tau3 - M32/n;
            f = 1 - (-a)/r2*(1 - cosh(F32));
            g = tau3 - sqrt((-a)^3/398600)*(F32 - sinhF32);
        end
        
    end %drfunc
    end %main
