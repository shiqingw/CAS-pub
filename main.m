%% Code for wellbore observer design (annulus)
% Description:
% This code simulates the dynamics in annular region with the boundaries
% given by the bit flow at x=0 and choke at x=1;
% Observer design in based on linearized model but implemented in the 
% non linear friction (Kelda) and choke models.


clear 
clc

%% Simulation

run_time = 20;      % total run time, sec
time_initial = 60;  % initial solution stabilization, sec

% choke ramps
choke_position_initial = 0.59; % initial choke position
choke_position_final = 0.6;   % final choke position  
start_chokeclosure = 0;  % begining of closure, sec
duration_chokeclosure = inf; % duration of the closing, sec

q_pump_initial = 0.03;   % pump flow, kg/s


%% Steady states

time_ss = 60; % steady state calulation stabilization, sec
choke_position_ss = 0.6; % steady state choke opening
q_pump_ss = 0.03; % flow at steady state

%% constants

p_atm = 1e5;    % Atmospheric pressure, bar
p0 = 1e5;       % reference pressure, bar
rho0 = 780;  % reference density, kg/m3
g = 9.81;       % acceleration due to gravity, m/s2

max_itr = 10;   % max. iteration for newton solver
tol = 1e-3;     % tolerance limit for newton solver

Kc = 2.85e-3;   % choke constant

%% wellbore and mud

CFL = 0.9;      % CFL condition for solver

MD = 2000;                      % Measured depth, m
TVD = 2000;                     % True vertical depth, m
theta = asin(TVD/MD) * 180/pi;  % pipe inclination, deg
N_cell = 200;                   % no. of discretization cells

dia_dp_in = 4.5 *0.0254;        % drill pipe inner dia, (inch -> m)
dia_dp_out = 5.5 *0.0254;       % drill pipe out dia (inch -> m)
dia_bit = 9 *0.0254;            % bit diameter (inch -> m)
area_nozzle = 1.1562 *0.0254^2; % bit nozzle diameter (in2 -> m2)
Cd = 0.8;                       % bit nozzle constant

compress_mud = 1.45e-9;         % compressibility of mud, bar-1
mu_mud = 40e-2;                 % viscosity of mud

% function that prepares various parameters of the annular region
[dx,dt,timesteps,well_mesh,in_mesh,bulk_mud,dia_hyd,area_ann,area_dp,Cl,A_plus,...
    A_minus,L_eigv,sig_pp,sig_pm,sig_mp,sig_mm] = prepare_well(MD,N_cell,CFL,...
    compress_mud,mu_mud,rho0,run_time,dia_bit,dia_dp_out,dia_dp_in,g);


%% steady state

% This function calculates the steady state values of the transport PDEs
[q1_ann_ss,q2_ann_ss] = CalculateSteadyState(q_pump_ss,area_ann,Cl,Kc,choke_position_ss,...
    well_mesh,in_mesh,dx,theta,dia_bit,g,mu_mud,rho0,p0,p_atm,A_plus,A_minus,L_eigv,...
dt,time_ss,max_itr,tol);

Qvec_ann_ss = [q1_ann_ss;
               q2_ann_ss];

% Riemann invariants
u_ann_ss = (Cl/2).*q1_ann_ss + (1/2).*q2_ann_ss;
v_ann_ss = -(Cl/2).*q1_ann_ss + (1/2).*q2_ann_ss;

pchoke_ss = (q1_ann_ss(end)-rho0)*Cl^2 + p0;    % steady state choke press
BHCP_ss = (q1_ann_ss(1)-rho0)*Cl^2 + p0;        % steady state BHCP
rho_bar = q1_ann_ss(1);                         % rho in choke eqn (approx)

%% Initialization

dq1_ann_plant = zeros(timesteps+1,length(well_mesh));
dq2_ann_plant = zeros(timesteps+1,length(well_mesh));

dpchoke_plant = zeros(1,timesteps+1);
pchoke_plant = zeros(1,timesteps+1);

BHCP_plant = zeros(1,timesteps+1);
dBHCP_plant = zeros(1,timesteps+1);

% This function initializes the states for the annulus and calculates
% the deviation from the steady state profile

[dq1_ann_plant(1,:),dq2_ann_plant(1,:),Qvec_initial] = ...
    initialization(q_pump_initial,area_ann,Cl,Kc,choke_position_initial,...
    well_mesh,in_mesh,dx,theta,dia_bit,g,mu_mud,rho0,p0,p_atm,A_plus,A_minus,L_eigv,...
dt,time_initial,q1_ann_ss,q2_ann_ss,max_itr,tol);


%% Plant
% This section of the code calculates the states of the plant. 

BHCP_plant(1) = (Qvec_initial(1,1)-rho0)*Cl^2 + p0;
dBHCP_plant(1) = dq1_ann_plant(1,1)*Cl^2;

dpchoke_plant(1) = dq1_ann_plant(1,end)*Cl^2;
pchoke_plant(1) = dpchoke_plant(1) + pchoke_ss;

choke_position = zeros(1,timesteps+1);
choke_position(1) = choke_position_initial;

CTE = 2*(rho0 * Cl^2 + (p_atm - p0));   % constant for simulation


%% Mapping for epsilon
% 
% dc = 0.2;
% choke_position_map = dc:dc:1;
% 
% epsilon_map = cell(length(choke_position_map),1);
% value_map = cell(length(choke_position_map),1);
% 
% H_bar = zeros(1,length(choke_position_map));
% b =2*( -rho0*Cl^2 + p0 - p_atm);
% 
% for i=1:length(choke_position_map)
% 
% [q1_ann_ss_map,q2_ann_ss_map] = CalculateSteadyState(q_pump_ss,area_ann,Cl,Kc,choke_position_map(i),...
%     well_mesh,in_mesh,dx,theta,dia_bit,g,mu_mud,rho0,p0,p_atm,A_plus,A_minus,L_eigv,...
% dt,time_ss,max_itr,tol);
% 
% a_t = Kc*choke_position_map(i);
% 
% q1 = q1_ann_ss_map(end);
% q2 = q2_ann_ss_map(end);
% 
% h1 = -q2/(q1^2)*area_ann + (a_t/sqrt( 2*Cl^2 + b/q1 ))* (b/(q1^2));
% h2 = area_ann/q1;
% 
% H_bar(i) = (h1+Cl*h2)/(h1-Cl*h2);
% % epsilon_map{i} = -1/Cl:1e-8:(H_bar(i)+1)/2/Cl;
% epsilon_map{i} =  -10/Cl:1e-8:(H_bar(i)+1)/2/Cl;
% value_map{i} = (H_bar(i)-epsilon_map{i}*Cl)./(1+epsilon_map{i}*(-Cl));
% 
% end
% 
% figure()
% 
% for i=1:length(choke_position_map)
%     
%    plot(epsilon_map{i,1},value_map{i,1},'linewidth',1.5);
%    hold all
%     
% end
% 
% grid on;



%% Wellbore observer - Kernel gains
% section that calculates the kernels required 

epsilon =8e-4; % design parameter in choke boundary

% constant at the bit boundary
q = (q_pump_initial/(area_ann*Cl) + 1)/(q_pump_initial/(area_ann*Cl) - 1);

% function provides the observer gains for the conservative variable
[L1,L2,L3,P_plus,P_minus] = WellboreKernels(well_mesh,sig_pp,sig_pm,sig_mp,sig_mm,Cl,q,...
    Kc,area_ann,rho0,p0,p_atm,u_ann_ss(end),v_ann_ss(end),q1_ann_ss(end),q2_ann_ss(end),choke_position_ss,epsilon);

%% Controler - Kernels
[K_vu,K_vv] = WellboreControlKernels(well_mesh,sig_pp,sig_pm,...
    sig_mp,sig_mm,Cl,q);

cst = -rho0*Cl^2 + p0 - p_atm; 
h1 = -q2_ann_ss(end)/(q1_ann_ss(end))^2*area_ann +...
    Kc*choke_position_ss*2*cst/sqrt(2*Cl^2+2*cst/q1_ann_ss(end))/(q1_ann_ss(end))^2;
h2 = area_ann/q1_ann_ss(end);
h3 = -Kc*sqrt(2*Cl^2+2*cst/q1_ann_ss(end));

C1 = -h1/(h3*Cl)-h2/h3;
C2 = h1/(h3*Cl)-h2/h3;

%% Observer 
% This section calculated the states by converging to the values in the
% plant by using the output error injection

dq1_ann_obs = zeros(timesteps+1,length(well_mesh));
dq2_ann_obs = zeros(timesteps+1,length(well_mesh));

% L1 = L1*0;
% L2 = L2*0;
% L3 = L3*0;

pchoke_obs = zeros(1,timesteps+1);
dpchoke_obs = zeros(1,timesteps+1);

BHCP_obs = zeros(1,timesteps+1);
dBHCP_obs = zeros(1,timesteps+1);

% initialize to wrong inital states
dq1_ann_obs(1,:) = .8*dq1_ann_plant(1,:);
dq2_ann_obs(1,:) = .8*dq2_ann_plant(1,:);

pchoke_obs(1) = (dq1_ann_obs(1,end)+q1_ann_ss(end) - rho0)*Cl^2 + p0;
dpchoke_obs(1) = dq1_ann_obs(1,end)*Cl^2;

BHCP_obs(1,1)  = (dq1_ann_obs(1,1)+q1_ann_ss(1)-rho0)*Cl^2+p0;
dBHCP_obs(1) = dq1_ann_obs(1,end)*Cl^2;



%% Time loop
for t=2:timesteps+1
    
    % Plant
    % choke ramp gives the corresponding opening of the choke
    choke_position(t) = choke_control(choke_position_initial,choke_position_final,choke_position_ss,...
        start_chokeclosure,duration_chokeclosure,t,dt,...
    K_vu, K_vv, dq1_ann_obs, dq2_ann_obs, C1, C2, Cl, well_mesh);
    % choke_position(t) = choke_position_ss;
    % q = dq + q_bar
    Qvec_ann_plant = [dq1_ann_plant(t-1,:) + Qvec_ann_ss(1,:);
                    dq2_ann_plant(t-1,:) + Qvec_ann_ss(2,:)];
    
    % function to calculate states at x=1, NL choke equation            
    Qvec_ann_plant(:,end) = choke_boundary_plant(Qvec_ann_plant,L_eigv,...
        choke_position(t),Kc,area_ann,Cl,CTE,max_itr,tol);
    
    % function to calculate states at x=0, bit
    Qvec_ann_plant(:,1) = pump_boundary(Qvec_ann_plant,L_eigv,...
        q_pump_initial,area_ann,max_itr,tol);
    
    % function that calcuates flux
    [F_plus_ann, F_minus_ann] = flux_vec(Qvec_ann_plant, A_plus , A_minus);
    
    % update using upwind method
    Qvec_ann_plant(:,in_mesh) = Qvec_ann_plant(:,in_mesh) - ...
        (dt/dx)*(F_plus_ann - F_minus_ann)+...
        dt*[zeros(1,length(in_mesh));
            source_terms(Qvec_ann_plant(:,in_mesh),theta,...
            area_ann,dia_bit,dx,g,mu_mud,1)];
    
    % update the variables    
    dq1_ann_plant(t,:) = Qvec_ann_plant(1,:) - Qvec_ann_ss(1,:);
    dq2_ann_plant(t,:) = Qvec_ann_plant(2,:) - Qvec_ann_ss(2,:);
    
    dpchoke_plant(t) = dq1_ann_plant(t,end)*Cl^2;
    pchoke_plant(t) = dpchoke_plant(t)+ pchoke_ss;
    
    dBHCP_plant(t) = dq1_ann_plant(t,1)*Cl^2;
    BHCP_plant(t) = dBHCP_plant(t) + BHCP_ss;
    
    
    
    
    %Observer
    
    Qvec_ann_obs = [dq1_ann_obs(t-1,:) + Qvec_ann_ss(1,:);
                dq2_ann_obs(t-1,:) + Qvec_ann_ss(2,:)];
            
    
    Qvec_ann_obs(:,end) = choke_boundary_obs(Qvec_ann_obs,L_eigv,...
        choke_position(t),Kc,area_ann,Cl,CTE,max_itr,tol,...
        pchoke_plant(t-1),L3,rho0,p0);
    
    
    Qvec_ann_obs(:,1) = pump_boundary(Qvec_ann_obs,L_eigv,...
        q_pump_initial,area_ann,max_itr,tol);
    
    [F_plus_ann, F_minus_ann] = flux_vec(Qvec_ann_obs, A_plus , A_minus);
    
    Qvec_ann_obs(:,in_mesh) = Qvec_ann_obs(:,in_mesh) - ...
        (dt/dx)*(F_plus_ann - F_minus_ann)+...
        dt*[( (pchoke_plant(t-1)-pchoke_obs(t-1) ).*L1) ;
            (source_terms(Qvec_ann_obs(:,in_mesh),theta,area_ann,dia_bit,...
            dx,g,mu_mud,1) + ( pchoke_plant(t-1)-pchoke_obs(t-1) ).*L2) ];     
    
    dq1_ann_obs(t,:) = Qvec_ann_obs(1,:) - Qvec_ann_ss(1,:);
    dq2_ann_obs(t,:) = Qvec_ann_obs(2,:) - Qvec_ann_ss(2,:);
  
    pchoke_obs(t) = (Qvec_ann_obs(1,end)-rho0)*Cl^2 + p0;
    
    dBHCP_obs(t) = dq1_ann_obs(t,1)*Cl^2;
    BHCP_obs(t) = dBHCP_obs(t) + BHCP_ss;
    
    
   
    
end





%% Figures

time = 0:dt:timesteps*dt;

figure(1)
plot(time,BHCP_plant)
hold on
plot(time,BHCP_obs)
hold off
xlabel('Time, s')
ylabel('Bottom hole pressure, bar')
legend('Plant','Observer')

figure(2)
plot(time,choke_position)
legend('Choke openning')