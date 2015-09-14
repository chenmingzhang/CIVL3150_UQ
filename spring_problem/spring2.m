%% This is a script to solve mass damping problem described at
%   Q5 Tutorial 1, CIVL3150. This script serves an example 
%   to solve ordinary differential equation using MATLAB
%  Some good habits to follow when developing your MATLAB script:
%  1. write scripts in three sections: (a) initialization
%    (b) main calculation loop and (c) result analysis
%  2. Declare the meaning and unit of the variable as
%    comment next to it.
%
%% ------- initialize simulation -----
% clear all declared variables
clear

% damper coefficient (kg/sec)
B = 1;

% spring constant (N/m) or (kg m /s2 )
k = 3;

% mass of the object (kg)
M = 5;

% time step (s), remember this value is a numerical parameter
%  Simulation result should always be independent from this parameter.
%  one way to check the dependency is to run the simulation with
%  different dt values and check the difference of the results
dt = 0.1;

% duration of simulation  (s)
tsim = 200;

% number of time steps estimated
nt = tsim/dt+1;

% t_ftd is x-axis (i.e., time) of the key points in the force-time 
%   diagram (s). notice that the last x-axis is the duration of 
%   simulation
t_ftd = [0,10,15,tsim];

% F_ftd is y-axis (i.e., force) of the key points in the force-time 
%   diagram (N or kgms-2). notice that the last y-axis is the force
%   at end of the simulation, which is zero
F_ftd = [5,5,0,0];

% defining arrays. Anything that changes with time can be defined
%   as an array (e.g., flow, converter or stock), the advantage of 
%   using array to store the values is that later one can plot and
%   analyze the change of value over time.
t = zeros(nt,1);
x = zeros(nt,1);
p = zeros(nt,1);
F = zeros(nt,1);
% an array for velocity (m/s), the inital value is 0
% note that the v array here is a converter, not a stock
% more information:
%  converter is a property written as a function of stock, and 
%   it changes over time. The common mistake is that converter
%   is written as function of flow.
%  For example: for the ice melting project, the surface area of
%   ice is a convert, as it changes overtime and it varies with
%   the mass of ice (which is a stoke). However, we never write
%   the surface area of ice as a function of melting rate, which 
%   a flow.
v = zeros(nt,1);

% an array for time (s), the initial time is 0. 
t(1) = 0;

% an array for mass displacement (m), the initial value is 2.
% mass displacement is the first stoke in the system
x(1) = 0;
 
% an array the momentum of the mass
p(1) = 0;

%% ---- main loop, use EULER method to solve the ordinary ---
%  ---- differential equation---
% equation

% defining loop for time
for i = 2:nt
    % the time point (s) of the current time step
    t(i) = t(i-1)+dt;
    % this line reveals that v is a converter as it is a function
    %  of p, which is a stock.
    v(i) = p(i-1)/M;


    % this  subsection gets the force at a time step based on
    % the force-time diagram given in the question brief
    % How it works:
    % cut the diagram into segments which can be represented by a line
    % function. for example, the line made by point (10,5) and (15,0) is the 
    % second segment. 
    % This segment can be expressed by force= ( 5-0)/(15-10) * (time-15) +0
    % The forward time index (F_fIdx) and backward time index (F_bIdx) for 
    % any time between 10 and 15s are 2 and 3. These indices are used to 
    % identify the segments that the current time locates
    % the forward index of time
    % think about how to use similar way to represent flow associated
    % with rainfall in project 2, especially when the simulation needs to 
    % run up to the second year? (function mod would be useful)
    F_fIdx = find((t_ftd-t(i))>0,1);
    % the backward index of time
    F_bIdx = F_fIdx-1;
    % the slope of the line segment
    k_f    = (F_ftd(F_fIdx)-F_ftd(F_bIdx)) / (t_ftd(F_fIdx)-t_ftd(F_bIdx))  ;
    % the force of time step i-1
    F(i-1) = k_f*(t(i-1)-t_ftd(F_fIdx))+F_ftd(F_fIdx) ;



    % net flows associated with the momentum balance
    % as it is using EULER method, all the flow are based on 
    % previous time step; if runge-kutta method is used, this
    % script has to be modified.
    f_n  = F(i-1)-k*x(i-1)-B*v(i-1);

    % calculate the momentum of the next time step (n+1)
    p(i) = p(i-1)+f_n*dt;
    % calculate the displacement of the next time step (n+1)
    x(i) = x(i-1)+v(i-1)*dt;
end


%% --------------- result analysis -------------
% create a figure handle for drawing
a.fig     = figure;
% fontsize of the figure
fontsize  = 12;
% linewidth of the figure
linewidth = 2;

% divide the plot region in 3 rows in a column, plot on the 
% first row
subplot(3,1,1)
% plot x as a function of t
% 'r-o' -- line in red (r), solid line (-) with circles (o)
% specify line width as the variable 'linewidth'
plot(t,x,'r-o','linewidth',linewidth)
% put label in x axis
xlabel('Time (s)','FontSize',fontsize,'FontWeight','bold')
% put label in y axis
ylabel('displacement (m)','FontSize',fontsize,'FontWeight','bold')
% give a title to the plot
title('Displacement of mass over time');


subplot(3,1,2)
plot(t,p,'r-o','linewidth',linewidth)
xlabel('Time (s)','FontSize',fontsize,'FontWeight','bold')
ylabel('displacement (m)','FontSize',fontsize,'FontWeight','bold')
title('Momentum of mass over time');


subplot(3,1,3)
plot(t,v,'r-o','linewidth',linewidth)
xlabel('Time (s)','FontSize',fontsize,'FontWeight','bold')
ylabel('Velocity (m/s)','FontSize',fontsize,'FontWeight','bold')
title('Velocity of mass over time');

print(a.fig,'mass.png','-dpng')   % png figure output
print(a.fig,'mass.tif','-dtiff','-r70')   % tif figure output
print(a.fig,'mass.eps','-depsc')  % vectorized figure output (publishing standard)
                               % it can be viewed using ghostscript, evince
			       % or by insterting it into word
saveas(a.fig,'mass.fig','fig') % save the matlab figure
