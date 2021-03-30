%% Script Description
% This script is for generating figures to compare the high thrust transfer
% options. Edit the propulsion system parameters in TransferCompare.m to
% get different results in this plotting script. 
%
% Developed by William Parker - March 2021

clc; clear all; close all; 

trans_type = [1,2,3,4,5]; %simulate over all of the possible transfer types
thrust_level = [1,1,2,1,2];%the only transfer that uses low thrust is the spiral transfer (type 3)and low thrust plane change
delta_RAAN = [0.1:0.1:5];%in degrees

%preallocate storage for variables
t_trans = zeros(length(trans_type), length(delta_RAAN));
dv = zeros(length(trans_type), length(delta_RAAN));
num_trans = zeros(length(trans_type), length(delta_RAAN));

for i = 1:length(trans_type)
    for j = 1:length(delta_RAAN)
        [t_trans(i,j), dv(i,j),num_trans(i,j)] = TransferCompare(trans_type(i), thrust_level(i), delta_RAAN(j), i+j);
    end
end

subplot(3,1,1)
semilogy(delta_RAAN, t_trans(1,:)./(86164.1), 'b-', 'LineWidth', 2);% Circular GOM
hold on
semilogy(delta_RAAN, t_trans(2,:)./86164.1, 'b--','LineWidth', 2);% Elliptical trans
semilogy(delta_RAAN, t_trans(4,:)./86164.1, 'b:','LineWidth', 2);% Plane change
semilogy(delta_RAAN, t_trans(3,:)./86164.1, 'r-','LineWidth', 2);% Spiral Trans
semilogy(delta_RAAN, t_trans(5,:)./86164.1, 'r--','LineWidth', 2);% low thrust Plane change

xlabel('Required Change in RAAN [deg]');
ylabel('Total Time for Transition [days]');
legend('Circular GOM (HT)','Elliptical Transfer Orbit (HT)', 'Plane Change (HT)', 'Spiral Transfer (LT)', 'Plane Change (LT)');

subplot(3,1,2)
semilogy(delta_RAAN, dv(1,:).*1000, 'b-','LineWidth', 2);
hold on
semilogy(delta_RAAN, dv(2,:).*1000,'b--','LineWidth', 2);
semilogy(delta_RAAN, dv(4,:).*1000, 'b:','LineWidth', 2);
semilogy(delta_RAAN, dv(3,:).*1000, 'r-','LineWidth', 2);
semilogy(delta_RAAN, dv(5,:).*1000, 'r--','LineWidth', 2);

xlabel('Required Change in RAAN [deg]')
ylabel('Delta-V Required for Maneuver [m/s]');

subplot(3,1,3)
semilogy(delta_RAAN, num_trans(1,:), 'b-','LineWidth', 2);
hold on
semilogy(delta_RAAN, num_trans(2,:), 'b--','LineWidth', 2);
semilogy(delta_RAAN, num_trans(4,:), 'b:','LineWidth', 2);
semilogy(delta_RAAN, num_trans(3,:), 'r-','LineWidth', 2);
semilogy(delta_RAAN, num_trans(5,:), 'r--','LineWidth', 2);

xlabel('Required Change in RAAN [deg]');
ylabel('Number of Transfers Possible');


