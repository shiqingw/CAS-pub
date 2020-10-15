function jacob = giveJ(pr,pn,in,ti,mi,ai,ji,si,ni,vi,theta)
% This funciton returns the jacobian for the newton solver
% of the form,
%
% | 1  0  0  0.... 0 |  
% | J1 J2 J3 0 ....0 |           
% | 0  J1 J2 J3 0 .0 |           
% |      . . .       |    
% |       . . .      |         
% | 0 0 ...  J1 J2 J3|        
% | 0 0 ...  0  -1 1 | 


    J1 = -2.*vi.*theta^2.*pr(in-1)' - 2.*vi.*theta.*(1-theta).*pn(in-1)' - ni.*theta;
    J2 = 2.*(mi.*ti.*theta + si.*theta^2 + vi.*theta^2).*pr(in)' + ...
        (mi.*ti.*(1-2.*theta) + 2.*si.*theta.*(1-theta) + 2.*vi.*theta.*(1-theta)).*pn(in)'+...
        ai.*mi + ni.*theta + ji.*theta;
    J3 = -2.*si.*theta^2.*pr(in+1)' - 2.*si.*theta.*(1-theta).*pn(in+1)' - ji.*theta;

    jacob = diag([J1,-1],-1) + diag([1,J2,1]) + diag([0,J3],1);


