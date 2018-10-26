%%% the script below draws the hydrograph using a gammar function
clear all
close all
%%% constants
g_per_kg       = 1000.;
m_per_km       = 1000.;
second_per_day = 86400.;

%%% Input parameters
simulation_time_sec          = 50 * second_per_day;


% numerical parameters
dt_second                    = 10000;  % time step


%%% parameters for unit hydrograph Q1
t_day        = 0:0.1:36; % time line in days
hydro_alpha1 = 2.5;  % larger shift higher
hydro_beta1  = 1.8;
flow_amp     = 2000;

Q1_hydrograph_inflow_m3Ps        = flow_amp*(t_day/hydro_beta1).^(hydro_alpha1-1) .* exp(-t_day/hydro_beta1) / hydro_beta1 / gamma(hydro_alpha1);

number_of_time_steps            = ceil(simulation_time_sec / dt_second);   % number of time steps

% time array are stored in a row
time_array_second(1:number_of_time_steps) = ((1:number_of_time_steps)-1)*dt_second;


% enable the below line to allow flow input from hydrograph 
% 1/second_per_day is equivalent to delta in the project description
Q2_hydrograph_inflow_m3Ps      = flow_amp*(time_array_second/second_per_day/hydro_beta1).^(hydro_alpha1-1) .* exp(-time_array_second/second_per_day/hydro_beta1) ...
    / hydro_beta1 / gamma(hydro_alpha1);

%%% enable below to plot hydrograph
figure
set(gcf,'color','w');
plot(t_day,Q1_hydrograph_inflow_m3Ps)
hold on    % make sure two lines can be plotted at the same time
plot(time_array_second/second_per_day,Q2_hydrograph_inflow_m3Ps,'linewidth',2)
grid on    % show the grid
xlabel('Time (days)');  
ylabel('Volumetric Inflow(m3/s)');
set(gca,'linewidth',1)   % thicken the bounding box
img = getframe(gcf);
imwrite(img.cdata, ['hydrograph', '.png']);  % save the figure into png file



