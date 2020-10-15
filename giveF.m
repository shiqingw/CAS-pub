function F = giveF(pr,pn,in,rmesh,ti,mi,ai,ji,si,ni,vi,theta,BHCP)
% This function returns the function to be solved for the newton solver

% Defining funciton in matrix form

% Left boundary: P1^(n+1) - BHCP = 0
% Right boundary: flux is zero. Therefore,
%                 PN+1 - PN = 0
    
% | (  (  0  0.... 0 | | P1  |^(n+1)   % | 0  0  0  0.... 0 | | P1  |^(n)  
% | A1 A2 A3 0 ....0 | | P2  |           | B1 B2 B3 0 ....0 | | P2  |
% | 0  A1 A2 A3 0 .0 | | P3  |           | 0  B1 B2 B3 0 .0 | | P3  |
% |      . . .       | | ..  |       +   |      . . .       | | ..  |  =
% |       . . .      | | ..  |           |       . . .      | | ..  |
% | 0 0 ...  A1 A2 A3| | PN  |           | 0 0 ...  B1 B2 B3| | PN  |
% | 0 0 ...  0  -1 1 | | PN+1|           | 0 0 ...  0  0  0 | | PN+1|
%
%           A(P^(n+1)) * P^(n+1)       +              B(P^n) * P^n  

% | cont |
% | 0    |
% | ..   |
% | ..   |                          = C 
% | ..   |
% | ..   |
% | 0    |


    A1 = -vi.*theta^2.*pr(in-1)' - 2.*vi.*theta.*(1-theta).*pn(in-1)' - ni.*theta;
    A2 = (mi.*ti.*theta + si.*theta^2 + vi.*theta^2).*pr(in)' + ...
        (mi.*ti.*(1-2.*theta) + 2.*si.*theta.*(1-theta) + 2.*vi.*theta.*(1-theta)).*pn(in)'+...
        ai.*mi + ni.*theta + ji.*theta;
    A3 = -si.*theta^2.*pr(in+1)' - 2.*si.*theta.*(1-theta).*pn(in+1)' - ji.*theta;

    B1 = -vi.*(1-theta)^2.*pn(in-1)' - ni.*(1-theta);
    B2 = (-mi.*ti.*(1-theta) + si.*(1-theta)^2 + vi.*(1-theta)^2).*pn(in)' ...
        - ai.*mi + ni.*(1-theta) + ji.*(1-theta);
    B3 = -si.*(1-theta)^2.*pn(in+1)' - ji.*(1-theta);
    


    A = diag([A1,-1],-1) + diag([1,A2,1]) + diag([0,A3],1);
    B = diag([B1,0],-1) + diag([0,B2,0]) + diag([0,B3],1);
    
    C = [BHCP, zeros(1,length(rmesh(2:end)))]';
    
    F = A*pr + B*pn - C;  % F = A(P^(n+1)) * P^(n+1) 
                          %     + B(P^n) * P^n - C 

    F(1) = (-xi/(area_ann*dr))*pr(1)^2 + ( (q_pump_initial/area_ann) +(beta*xi/(area_ann*dr)) -Cl )*pr(2) ...
            + (xi/area_ann*dr)*pr(2)*pr(1) - (beta*xi/(area_ann*dr))*pr(2) - 2*pr(end)*Cl^2;
    
    f1 = pr(end-1) - pr(end) - (pr(1)-beta)/Cl;
    f2 = 


% % Defining function as F(P^n, P^(n+1))=0
%
%     F = (mi.*ti.*theta).*pr(in)'.^2 + ((mi.*ti + 2.*si.*theta.*(1-theta) +...
%         2.*vi.*theta.*(1-theta)).*pn(in)' + ai.*mi + ni.*theta + ji.*theta).*pr(in)'...
%         + (mi.*ti.*(1-theta) + si.*(1-theta)^2 + vi.*(1-theta)^2).*pn(in)'.^2+...
%         (ai.*mi + ni.*(1-theta) + ji.*(1-theta)).*pn(in)' ...
%         - si.*(theta.*pr(in+1)' + (1-theta).*pn(in+1)').^2 ...
%         - vi.*(theta.*pr(in-1)' + (1-theta).*pn(in-1)').^2 ...
%         - ji.*(theta.*pr(in+1)' + (1-theta).*pn(in+1)') ...
%         - ni.*(theta.*pr(in-1)' + (1-theta).*pn(in-1)');
%
%     F = [(pr(1)-BHCP), F, (pr(end)-pr(end-1)) ]'; % Taking BC into
%                                                   % account
%     

end
