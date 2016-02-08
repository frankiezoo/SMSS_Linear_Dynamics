function sys = lin_matlab(sim)

%% Exact linearization of the Simulink model SMSS_nonlinear
%
% This MATLAB script is the command line equivalent of the exact
% linearization tab in linear analysis tool with current settings.
% It produces the exact same linearization results as hitting the Linearize button.

% MATLAB(R) file generated by MATLAB(R) 8.6 and Simulink Control Design (TM) 4.2.1.
%
% Generated on: 25-Jan-2016 17:30:18

%% Specify the model name
model = sim;

%% Specify the analysis I/Os
% Get the analysis I/Os from the model
io = getlinio(model);

%% Specify the operating point
% Use the model initial condition
op = operpoint(model);


%% Linearize the model
sys = linearize(model,io,op);

%% Plot the resulting linearization
% step(sys)

