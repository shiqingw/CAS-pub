function ch_op = choke_control(initial_pos,final_pos,choke_position_ss, start_time,duration,t,dt,...
    Kvu, Kvv, dq1_ann_obs, dq2_ann_obs, C1, C2, Cl, well_mesh)


if t*dt<start_time
    ch_op = initial_pos;


elseif t*dt <= start_time+duration
   n = length(well_mesh)-2;
   du_obs = (Cl/2).*dq1_ann_obs(t-1,:) + (1/2).*dq2_ann_obs(t-1,:);
   dv_obs = -(Cl/2).*dq1_ann_obs(t-1,:) + (1/2).*dq2_ann_obs(t-1,:);
   
   integartion_part = trapz(well_mesh(1:n),transpose(Kvu(1:n,n)).*du_obs(2:n+1)) +...
       trapz(well_mesh(1:n),transpose(Kvv(1:n,n)).*dv_obs(2:n+1));
   dZc = C1*du_obs(n) + C2*integartion_part ;
   ch_op = max(0,min(1,choke_position_ss + dZc));
    

else %t*dt>start_time+duration
    ch_op = final_pos;
end

end