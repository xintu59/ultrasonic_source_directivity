clear
close all

%% define material properties
% Test case 1
% need to specify the elastic constants, C, density, rho, and layer thickness
% for each layer
rotation1 = 0;
rotation2 = 0.5*pi;
rotation3 = 0;

% c11 = 126.83e9/6; % Pa
% c66 = c11/2.01;
% c12 = c11-2*c66;
% c13 = c12*1.1;
% c44 = c66*1.4;
% c33 = c11*1.9;
% rho0 = 1570;
% C0 = [c11, c12, c13, 0, 0, 0
%         c12, c11, c13, 0, 0, 0
%         c13, c13, c33, 0, 0, 0
%         0, 0, 0, c44, 0, 0
%         0, 0, 0, 0, c44, 0
%         0, 0, 0, 0, 0, c66];
% C0 = fn_rotate_stiffness_matrix(C0, rotation1, rotation2, rotation3); % layer 0
% layer1_thickness = 1e-6;

% self-consistency check
% layer_num = 3; % including the two top and bottom half-spaces, layer 0 and layer n
% the layer thickness for layer 0 and layer n are unused (for placeholder purpose)
% for i = 1:layer_num
%     C(:,:,i) = C0;
%     rho(:,i) = rho0;
%     layer_thickness(:,i) = layer1_thickness;
% end

%% Test case 2
% c11 = 70.857e9;
% c12 = 7.2979e9;
% c13 = 10.186e9;
% c22 = 14e9;
% c33 = 63.772e9;
% c44 = 4.4063e9;
% c55 = 7.6455e9;
% c66 = 4.5140e9;

% c11 = 76.015e9;
% c12 = 7.4255e9;
% c13 = 14.475e9;
% c22 = 14e9;
% c33 = 50.036e9;
% c44 = 4.2704e9;
% c55 = 11.937e9;
% c66 = 4.6661e9;

c11 = 68.712e9;
c12 = 7.2979e9;
c13 = 12.331e9;
c22 = 14e9;
c33 = 61.626e9;
c44 = 4.4063e9;
c55 = 9.791e9;
c66 = 4.514e9;
C0 = [c11, c12, c13, 0, 0, 0
        c12, c22, c13, 0, 0, 0
        c13, c13, c33, 0, 0, 0
        0, 0, 0, c44, 0, 0
        0, 0, 0, 0, c55, 0
        0, 0, 0, 0, 0, c66];
rho0 = 1570; % kg/m3

E = 5e9;
rho1 = 1140;
v = 0.389;
c11 = E*(1-v)/((1+v)*(1-2*v));
c66 = E/(2*(1+v));
c12 = E*v/((1+v)*(1-2*v));
c13 = c12*0.99;
c33 = c11*0.99;
c44 = E/(2.01*(1+v));

C1 = [c11, c12, c13, 0, 0, 0
        c12, c11, c13, 0, 0, 0
        c13, c13, c33, 0, 0, 0
        0, 0, 0, c44, 0, 0
        0, 0, 0, 0, c44, 0
        0, 0, 0, 0, 0, c66];

material_properties(1).('density') = rho0;
material_properties(1).('elastic_constants') = C0;
material_properties(1).('layer_thickness') = 100e-6;

material_properties(2).('density') = rho1;
material_properties(2).('elastic_constants') = C1;
material_properties(2).('layer_thickness') = 100e-6;

% material_properties(3).('density') = rho0;
% material_properties(3).('elastic_constants') = C0;
% material_properties(3).('layer_thickness') = 100e-6;

%% incidence wave specification
u_n_depth_sweep = [];
source_depth_list = [0,20e-6,40e-6,60e-6,80e-6,100e-6];
for ii = 1:length(source_depth_list)
    source_depth = source_depth_list(ii);
frequency = 2.5e6; % Hz
in_wave = 3; % 1: qS1, 2:qS2, 3:qL
load_direction = 2; % 1. x; 2. y; 3. z
%source_depth = 100e-6; % unit: metres, greater than or equal to 0
source_type = 'centre_expansion'; % choose between 'centre_expansion', 'dipole', 'point_force'
boundary_type = 'free_surface'; % choose between 'free_surface' - no transmitted wave in layer n, 'bounded' - layer n is semi-infinite half-space
d = 5e-3; % unit: metres, observation distance
resolution = 361; % number of angle increments in the interval [-0.5*pi, 0.5*pi]
source_width = 0; % unit: metres
%%%%%%%%%%%%%%%%%%%%%%%%%% End of inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for i = 1:length(layer_thickness)
%     material_properties(i).('density') = rho(:,i);
%     material_properties(i).('elastic_constants') = C(:,:,i);
%     material_properties(i).('layer_thickness') = layer_thickness(:,i);
% end

material_properties = orderfields(material_properties);

tic
[ray_directions, u_n] = fn_directivity_general_case(frequency, in_wave, load_direction, source_depth, source_type, boundary_type, d, material_properties, resolution, source_width);
toc

figure(101);
scatter(ray_directions/pi*180,abs(u_n))
xlim([-95,95])
%ylim([0,0.2e-4])
crop_fig;
xlabel('Angles (\circ)')
ylabel('Amplitude')
hold on

u_n_depth_sweep = [u_n_depth_sweep, u_n];
end
