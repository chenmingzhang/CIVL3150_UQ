% the video explanation of this script is accessible from https://www.youtube.com/watch?v=-bPIko39J_E&feature=youtu.be
% in case if narration is slow, increase the playback speed to 1.25
clear

dz_m=0.1;
number_cells=8;
dx_m=1;
dy_m=1;
time_step_s = 30*60;
max_time_steps=24*10;



thermal_cond_wPmPc=2.0;
thermal_conductance_wPm2Pc=100;  %original 100
vol_heat_capacity_jPm3Pc=2959200;

T_c_ay=zeros(number_cells,max_time_steps);
heat_flux_out_jPc=zeros(number_cells,max_time_steps);
heat_flux_in=zeros(number_cells,max_time_steps);

T_c_ay(:,1)=15;  %original 15

T_surface_c=30; %original 30
T_bottom_c=15;%original 15

for time_idx=1:max_time_steps-1
    
    for space_idx=1:number_cells-1
        heat_flux_out_jPs(space_idx,time_idx+1)=thermal_cond_wPmPc*(T_c_ay(space_idx,time_idx)-T_c_ay(space_idx+1,time_idx)) *dx_m*dy_m/dz_m;
    end
    heat_flux_out_jPs(number_cells,time_idx+1)=thermal_conductance_wPm2Pc*(T_c_ay(number_cells,time_idx)-T_bottom_c)*dx_m*dy_m;
    for space_idx=2:number_cells
        heat_flux_in_jPs(space_idx,time_idx+1)=thermal_cond_wPmPc*(T_c_ay(space_idx-1,time_idx)-T_c_ay(space_idx,time_idx)) *dx_m*dy_m/dz_m;
    end
    heat_flux_in_jPs(1,time_idx+1)=thermal_conductance_wPm2Pc*(T_surface_c-T_c_ay(1,time_idx))*dx_m*dy_m;
    
    for space_idx=1:number_cells
        T_c_ay(space_idx,time_idx+1)= (heat_flux_in_jPs(space_idx,time_idx+1)-heat_flux_out_jPs(space_idx,time_idx+1))*time_step_s/vol_heat_capacity_jPm3Pc/dz_m/dx_m/dy_m +T_c_ay(space_idx,time_idx);
    end
    
end

time_s_ay=0:time_step_s:time_step_s*(max_time_steps-1);
x_m_ay=0:dz_m:dz_m*number_cells-dz_m;
fig=figure;
set(gcf,'color','w')
subplot(2,1,1)
plot(x_m_ay,T_c_ay(:,1),'r-','displayname','time = 0','linewidth',2);hold on
plot(x_m_ay,T_c_ay(:,24*4),'g-','displayname','time = 2 days','linewidth',2);hold on
plot(x_m_ay,T_c_ay(:,24*10),'b-','displayname','time = 10 days','linewidth',2);hold on
xlabel('DEPTH (m)')
ylabel('TEMP. {\circ}C')
legend show

dayPs=1/86400;

subplot(2,1,2)
plot(time_s_ay*dayPs,T_c_ay(1,:),'r-','displayname','time = 0','linewidth',2);hold on
plot(time_s_ay*dayPs,T_c_ay(3,:),'g-','displayname','time = 0.6 m','linewidth',2);hold on
plot(time_s_ay*dayPs,T_c_ay(end,:),'b-','displayname','time = 0.8 m','linewidth',2);hold on
xlabel('TIME (days)')
ylabel('TEMP. {\circ}C')
legend show

print(fig,'result.jpg','-djpeg','-r200')




