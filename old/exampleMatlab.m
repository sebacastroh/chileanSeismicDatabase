clear all
close all
clc

EarthquakeName = '20171013_4Mw_20.29S_69.14W_91KM.mat';
load(EarthquakeName);

variablesInCurrentWorkspace = who;
numVariables = length(variablesInCurrentWorkspace);

disp(strcat(['Event: ', EarthquakeName]));
disp(strcat(['Number of stations:', ' ', num2str(numVariables)]));

% Let's plot the acceleration, velocity and displacement of the first component for the first station
attributes = fieldnames(st00);

disp('Attributes of station st00:')
for i=1:length(attributes)
    varname = char(attributes(i));
    value = getfield(st00, varname);
    if varname(1:2) == 'ac'
        disp(strcat(['st00.', varname, ' = ', 'vector of ', num2str(length(value)), ' elements']))
    elseif isstring(value)
        disp(strcat(['st00.', varname, ' = ', value]))
    else
        disp(strcat(['st00.', varname, ' = ', num2str(value)]))
    end
end

t = linspace(0., (length(st00.acc_1)-1)*st00.dt, length(st00.acc_1));

figure

subplot(311)
plot(t, st00.acc_1/9.81)
ylabel('a(t) [g]')
title(strcat(['Station ', st00.name, newline 'Component ', st00.component_1]))

subplot(312)
v = cumtrapz(t, st00.acc_1);
plot(t, v)
ylabel('v(t) [m/s]')

subplot(313)
d = cumtrapz(t, v);
plot(t, d)
ylabel('d(t) [m]')
xlabel('Time [s]')
