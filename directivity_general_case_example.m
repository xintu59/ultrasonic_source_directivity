clear
close all

%% define material properties
% Test case 1
% need to specify the elastic constants, C, density, rho, and layer thickness
% for each layer
rotation1 = 0*pi;
rotation2 = 0.5*pi;
rotation3 = 0;

c11 = 126.83e9/6; % Pa
c66 = c11/2.01;
c12 = c11-2*c66;
c13 = c12*1.1;
c44 = c66*1.4;
c33 = c11*1.9;

% c11 = 51.2e9; % Pa qL velocity in in-plane direction
% c66 = c11*0.3; % qS velocity in in-plane direction
% c13 = c11*0.2; % irrelevant if c13 < c11*0.25- qS1
% c44 = c11*0.18; % qS velocity in through thickness direction
% c33 = c11*0.25; % qL velocity in through thickness direction
% c12 = c11-2*c66; % 
rho0 = 1570;
C0 = [c11, c12, c13, 0, 0, 0
        c12, c11, c13, 0, 0, 0
        c13, c13, c33, 0, 0, 0
        0, 0, 0, c44, 0, 0
        0, 0, 0, 0, c44, 0
        0, 0, 0, 0, 0, c66];
C0 = fn_rotate_stiffness_matrix(C0, rotation1, rotation2, rotation3); % layer 0
% C0 =[63056054687.5000	7250000000.00000	8043945312.50000	2.55996691544737e-09	8.27740445865492e-06	-4.11301351081880e-08
% 7250000000.00000	12800000000.0000	7250000000.00000	-4.36901020236353e-08	4.59242549680257e-08	4.36901020236353e-08
% 8043945312.50000	7250000000.00000	63056054687.5000	4.11301351081880e-08	-5.19741775879932e-06	-2.55996691544738e-09
% 2.55996691544737e-09	-4.36901020236353e-08	4.11301351081879e-08	4459459459.45946	-1.04184953162254e-24	5.16605936090282e-08
% 8.27740445865492e-06	4.59242549680257e-08	-5.19741775879932e-06	-7.24197771741850e-24	5500000000.00000	-4.56673053753518e-24
% -4.11301351081879e-08	4.36901020236353e-08	-2.55996691544738e-09	5.16605936090282e-08	-1.07668587233312e-23	4459459459.45946];
% C1 = C0;
% rho1 = rho0;
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

% c11 = 68.712e9;
% c12 = 7.2979e9;
% c13 = 12.331e9;
% c22 = 14e9;
% c33 = 61.626e9;
% c44 = 4.4063e9;
% c55 = 9.791e9;
% c66 = 4.514e9;
% C0 = [c11, c12, c13, 0, 0, 0
%         c12, c22, c13, 0, 0, 0
%         c13, c13, c33, 0, 0, 0
%         0, 0, 0, c44, 0, 0
%         0, 0, 0, 0, c55, 0
%         0, 0, 0, 0, 0, c66];
% C0 = fn_rotate_stiffness_matrix(C0, rotation1, rotation2, rotation3);
% C0=[56854910714.2857	7250000000.00000	20645089285.7143	4.58621808058780e-24	-2.21842386112878e-06	-9.82634728111692e-24
% 7250000000.00000	14000000000.0000	7250000000.00000	2.77807850389841e-23	-2.56101298433377e-08	3.75967295650474e-23
% 20645089285.7143	7250000000.00000	56854910714.2857	-3.92590019838996e-24	1.00629136492365e-06	1.96288060276560e-23
% 0	0	0	4459459459.45946	6.82207133393163e-24	-3.84273507049481e-08
% 0	0	0	4.39251381950910e-24	18104910714.2857	2.49350234291421e-24
% 0	0	-1.88257350626803e-23	-3.84273507049481e-08	5.69664614810139e-24	4459459459.45946];
% C0= [69459821428.5714	7250000000.00000	8040178571.42857	7.97954352803207e-24	2.93177163411950e-06	1.68629870049420e-24
% 7250000000.00000	14000000000.0000	7250000000.00000	3.14775757425585e-23	-4.59242549680257e-08	3.14775757425585e-23
% 8040178571.42857	7250000000.00000	69459821428.5714	1.68629870049420e-24	-6.33016650175341e-06	7.97954352803207e-24
% 8.33724321278399e-24	2.48005149605226e-23	8.00565979777813e-24	4459459459.45946	6.20012818579596e-24	-5.16605936090282e-08
% 2.93177163411950e-06	-4.59242549680257e-08	-6.33016650175341e-06	-3.10201199115528e-39	5500000000.00000	-3.88292940420090e-40
% 8.00565979777813e-24	2.48005149605226e-23	8.33724321278400e-24	-5.16605936090282e-08	6.20012818579597e-24	4459459459.45946];
% rho0 = 1570; % kg/m3
% C1 = C0;
% rho1 = rho0;

E = 50e9;
rho1 = 1140;
v = 0.389;
c11 = E*(1-v)/((1+v)*(1-2*v));
c66 = E/(2*(1+v));
c12 = E*v/((1+v)*(1-2*v));
c13 = c12*0.999;
c33 = c11*0.999;
c44 = E/(2.001*(1+v));

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

material_properties(3).('density') = rho1;
material_properties(3).('elastic_constants') = C1;
material_properties(3).('layer_thickness') = 50e-6;
% 
% material_properties(4).('density') = rho1;
% material_properties(4).('elastic_constants') = C1;
% material_properties(4).('layer_thickness') = 100e-6;

%% incidence wave specification
frequency = 10e6; % Hz
in_wave = 3; % 1: qS1, 2:qS2, 3:qL
load_direction = 2; % 1. x; 2. y; 3. z
source_depth = 0; % unit: metres, greater than or equal to 0
source_type = 'point_force'; % choose between 'centre_expansion', 'dipole', 'point_force'
boundary_type = 'free_surface'; % choose between 'free_surface' - no transmitted wave in layer n, 'bounded' - layer n is semi-infinite half-space
d = 5e-3; % unit: metres, observation distance
resolution = 361; % number of angle increments in the interval [-0.5*pi, 0.5*pi]
source_width = 0; % unit: metres
dimensions = '2D';
%%%%%%%%%%%%%%%%%%%%%%%%%% End of inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for i = 1:length(layer_thickness)
%     material_properties(i).('density') = rho(:,i);
%     material_properties(i).('elastic_constants') = C(:,:,i);
%     material_properties(i).('layer_thickness') = layer_thickness(:,i);
% end

material_properties = orderfields(material_properties);

tic
[ray_directions, u_n] = fn_directivity_general_case(frequency, in_wave, load_direction, source_depth, source_type, boundary_type, d, material_properties, resolution, source_width, dimensions);
toc

figure(101);
scatter(ray_directions/pi*180,abs(u_n)) %1/3669*2*pi*frequency %3669 %5495
xlim([-95,95])
%ylim([0,0.2e-4])
crop_fig;
xlabel('Angle (\circ)')
ylabel('Amplitude')

