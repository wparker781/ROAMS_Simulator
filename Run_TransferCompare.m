%% Script Description
% This script is for generating figures to compare the high thrust transfer
% options. Edit the propulsion system parameters in TransferCompare.m to
% get different results in this plotting script. 
%
% William Parker - March 2021

clc; clear all; close all; 

trans_type = [1,2,4];
thrust_level = 1;
delta_RAAN = [0.1:0.1:10];%in degrees

%preallocate storage for variables
t_trans = zeros(length(trans_type), length(delta_RAAN));
dv = zeros(length(trans_type), length(delta_RAAN));
num_trans = zeros(length(trans_type), length(delta_RAAN));

for i = 1:length(trans_type)
    for j = 1:length(delta_RAAN)
        [t_trans(i,j), dv(i,j),num_trans(i,j)] = TransferCompare(trans_type(i), thrust_level, delta_RAAN(j));
    end
end

subplot(3,1,1)
semilogy(delta_RAAN, t_trans(1,:)./3600, 'b');% Circular GOM
hold on
semilogy(delta_RAAN, t_trans(2,:)./3600, 'r');% Elliptical trans
semilogy(delta_RAAN, t_trans(3,:)./3600, 'g');% Plane change
xlabel('Required Change in RAAN [deg]');
ylabel('Total Time for Transition [hours]');
legend('Circular GOM','Elliptical Transfer Orbit', 'Plane Change Maneuver');

subplot(3,1,2)
semilogy(delta_RAAN, dv(1,:).*1000, 'b');
hold on
semilogy(delta_RAAN, dv(2,:).*1000,'r');
semilogy(delta_RAAN, dv(3,:).*1000, 'g');
xlabel('Required Change in RAAN [deg]')
ylabel('Delta-V Required for Maneuver [m/s]');
% legend('Circular GOM','Elliptical Transfer Orbit', 'Plane Change Maneuver');

subplot(3,1,3)
plot(delta_RAAN, num_trans(1,:), 'b');
hold on
plot(delta_RAAN, num_trans(2,:), 'r');
plot(delta_RAAN, num_trans(3,:), 'g');
xlabel('Required Change in RAAN [deg]');
ylabel('Number of Transfers Possible');
% legend('Circular GOM','Elliptical Transfer Orbit', 'Plane Change Maneuver');


