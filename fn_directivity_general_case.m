function [ray_directions, u_n] = fn_directivity_general_case(frequency, in_wave, load_direction, source_depth, source_type, boundary_type, d, material_properties, resolution, source_width,dimension)
%SUMMARY
%   Find the directivity of ultrasonic waves with different source types
%   under different boundary conditions
%USAGE
%   [ray_directions, u_n] = fn_directivity_general_case(frequency, in_wave, load_direction, source_depth, source_type, boundary_type, d, layer_num, C0, C1, C2, Cn)
%INPUTS
%   frequency - incident wave frequency
%   in_wave - incident wave mode
%   load_direction - 1. x; 2. y; 3. z
%   source_depth - unit: metres, greater than or equal to 0
%   source_type - choose between 'centre_expansion', 'dipole', 'point_force'
%   boundary_type - choose between 'free_surface' - no transmitted wave in layer n, 'bounded' - layer n is semi-infinite half-space
%   d - unit: metres, observation distance
%   material_properties - a structure containing stiffness matrix and densities for each layer
%   resolution
%   source_width
%OUTPUTS
%   ray_directions - group velocity angles
%   u_n - beam amplitudes
%AUTHOR
%   Lily Tu (2023)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = linspace(0, pi, resolution)';
size_a = length(a);

%% Layer 0 - bulk CFRP
C0 = material_properties(1).elastic_constants;
rho0 = material_properties(1).density;
switch dimension
    case '2D'
[Gd_in, p0_minus, ray_directions, slowness_total, n_total, ~] = fn_2D_anisotropic_greens(a, C0, rho0, frequency, in_wave, d);
    case '3D'
[Gd_in, p0_minus, ray_directions, slowness_total, n_total, ~] = fn_3D_anisotropic_greens(a, C0, rho0, frequency, in_wave, d);
end

slowness_i = cat(3,slowness_total.*n_total(:,1),slowness_total.*n_total(:,2),slowness_total.*n_total(:,3));

ki_x = squeeze(2*pi*frequency.*slowness_i(:,:,1));
ki_y = squeeze(2*pi*frequency.*slowness_i(:,:,2));
ki_z = squeeze(2*pi*frequency.*slowness_i(:,:,3));
ki = cat(3,ki_x,ki_y,ki_z); % wave vector (angles, wave modes, directions)
[p0_plus, slowness_r, ~, ~] = fn_reflected_polarisation_vector(C0, rho0, in_wave, n_total, slowness_total);

kr_x = squeeze(2*pi*frequency.*slowness_r(:,1,:));
kr_y = squeeze(2*pi*frequency.*slowness_r(:,2,:));
kr_z = squeeze(2*pi*frequency.*slowness_r(:,3,:));
kr = cat(3,kr_x,kr_y,kr_z);
% p_reflect is p0_plus
% p0_minus is the polarising vector for the incident wave, p_total(:,l,in_wave)

% calculate D0_plus
C0 = fn_voigt_to_tensor(C0);
D0_plus = zeros(size_a,3,3);
for i = 1:3 % stress direction
    for n = 1:3
        for l = 1:3
            for j = 1:3 % wave modes
                D0_plus(:,i,j) = D0_plus(:,i,j)+C0(2,i,n,l)*p0_plus(:,l,j).*kr(:,j,n);
            end
        end
    end
end

% calculate D0_minus
D0_minus = zeros(size_a,3,3);
for i = 1:3 % stress direction
    for n = 1:3
        for l = 1:3
            for j = 1:3 % wave modes
                D0_minus(:,i,j) = D0_minus(:,i,j)+C0(2,i,n,l)*p0_minus(:,l,j).*ki(:,j,n);
            end
        end
    end
end

%% a couple of edge points (at 0 and 90degs) throw the warning but does not affect the overall directivity calculation
MSGID='MATLAB:nearlySingularMatrix'; 
warning('off', MSGID)
% to turn the warning back on
%warning on verbose

%% calculate global stiffness matrix, K
K = zeros(size_a,6,6);
% K_11 = K(:,1:3,1:3); K_21 = K(:,1:3,4:6); K_12 = K(:,4:6,1:3); K_22 =
% K(:,4:6,4:6);
if length(material_properties)>2
    Ka = fn_layer_stiffness_matrix(material_properties(2).elastic_constants,...
        material_properties(2).density, material_properties(2).layer_thickness, in_wave, n_total, slowness_total, frequency, size_a);
    if length(material_properties)>3
        for j = 3:length(material_properties)-1
            Kb = fn_layer_stiffness_matrix(material_properties(j).elastic_constants,... 
                material_properties(j).density, material_properties(j).layer_thickness, in_wave, n_total, slowness_total, frequency, size_a);

            for i = 1:size_a
                K(i,:,:) = [squeeze(Ka(i,1:3,1:3))+squeeze(Ka(i,1:3,4:6))/(squeeze(Kb(i,1:3,1:3))-squeeze(Ka(i,4:6,4:6)))*squeeze(Ka(i,4:6,1:3)), ...
                    squeeze(-Ka(i,1:3,4:6))/(squeeze(Kb(i,1:3,1:3))-squeeze(Ka(i,4:6,4:6)))*squeeze(Kb(i,1:3,4:6));...
                    squeeze(Kb(i,4:6,1:3))/(squeeze(Kb(i,1:3,1:3))-squeeze(Ka(i,4:6,4:6)))*squeeze(Ka(i,4:6,1:3)),...
                    squeeze(Kb(i,4:6,4:6))-squeeze(Kb(i,4:6,1:3))/(squeeze(Kb(i,1:3,1:3))-squeeze(Ka(i,4:6,4:6)))*squeeze(Kb(i,1:3,4:6))];            
            end
            Ka = K;
        end
    else
        K = Ka;
    end
else
end

%% Layer n - calculate transmitted waves, Dn1_minus, pn1_minus
% calculate -ve polarisation vector - transmitted waves
Cn = material_properties(end).elastic_constants;
rhon = material_properties(end).density;
[pn1_minus, slowness_minus_xyz, ~, ~] = fn_transmitted_polarisation_vector(Cn, rhon, in_wave, n_total, slowness_total);
Cn = fn_voigt_to_tensor(Cn);

% calculate D_minus
Dn1_minus = zeros(size_a,3,3);
for i = 1:3 % stress direction
    for n = 1:3
        for l = 1:3
            for j = 1:3 % wave modes
                Dn1_minus(:,i,j) = Dn1_minus(:,i,j)+Cn(2,i,n,l)*pn1_minus(:,l,j)*2*pi*frequency.*slowness_minus_xyz(:,n,j);
            end
        end
    end
end

%% Solve for reflection and transmission coefficients, xr and xt
A_in = zeros(size_a,3);
A_in (:,in_wave) = A_in (:,in_wave) + 1;
xr = zeros(size(A_in));
xt = zeros(size(A_in));
switch boundary_type
    case 'free_surface'
        % calculate K_s, no Dn1_minus and pn1_minus term
        K_s = zeros(size_a,3,3);
        if (-squeeze(K(i,4:6,4:6))) == 0
        
        else
            for i = 1:size_a
            K_s(i,:,:) = squeeze(K(i,1:3,1:3)) + squeeze(K(i,1:3,4:6)) /(squeeze(Dn1_minus(i,:,:)*0)/squeeze(pn1_minus(i,:,:))-squeeze(K(i,4:6,4:6))) * squeeze(K(i,4:6,1:3));
            end
        end
        for i = 1:size_a
            B1 = (-squeeze(K_s(i,:,:))*squeeze(p0_minus(i,:,:))+squeeze(D0_minus(i,:,:)))*A_in(i,:)';
            A1 = squeeze(K_s(i,:,:))*squeeze(p0_plus(i,:,:))-squeeze(D0_plus(i,:,:));
            xr(i,:) = A1\B1;
        end
    case 'bounded'
        % calculate K_s
        K_s = zeros(size_a,3,3);
        for i = 1:size_a
            K_s(i,:,:) = squeeze(K(i,1:3,1:3)) + squeeze(K(i,1:3,4:6)) / (squeeze(Dn1_minus(i,:,:))/squeeze(pn1_minus(i,:,:))-squeeze(K(i,4:6,4:6))) * squeeze(K(i,4:6,1:3));
        end
        for i = 1:size_a
            B1 = (-squeeze(K_s(i,:,:))*squeeze(p0_minus(i,:,:))+squeeze(D0_minus(i,:,:)))*A_in(i,:)';
            A1 = squeeze(K_s(i,:,:))*squeeze(p0_plus(i,:,:))-squeeze(D0_plus(i,:,:));
            xr(i,:) = A1\B1;

            B2 = squeeze(K(i,4:6,1:3))*xr(i,:)';
            A2 = squeeze(Dn1_minus(i,:,:))-squeeze(K(i,4:6,4:6))*squeeze(pn1_minus(i,:,:));
            xt(i,:) = A2\B2;
        end
end

switch source_type
    case 'centre_expansion'
        u_n = xr(:,1).*Gd_in.*exp(1i.*kr(:,1,2).*source_depth).*1i.*dot(squeeze(kr(:,1,:)),p0_plus(:,:,1),2) + ...
           xr(:,2).*Gd_in.*exp(1i.*kr(:,2,2).*source_depth).*1i.*dot(squeeze(kr(:,2,:)),p0_plus(:,:,2),2) + ...
           xr(:,3).*Gd_in.*exp(1i.*kr(:,3,2).*source_depth).*1i.*dot(conj(squeeze(kr(:,3,:))),p0_plus(:,:,3),2) + ...
           Gd_in.*exp(1i.*ki(:,in_wave,2).*source_depth).*1i.*dot(squeeze(ki(:,in_wave,:)),squeeze(p0_minus(:,:,in_wave)),2);

    case 'dipole'
        u_n = xr(:,1).*Gd_in.* exp(1i.*kr(:,1,2).*source_depth).*1i.*kr(:,1,load_direction).*p0_plus(:,load_direction,1) + ...
           xr(:,2).*Gd_in.* exp(1i.*kr(:,2,2).*source_depth).*1i.*kr(:,2,load_direction).*p0_plus(:,load_direction,2) + ...
           xr(:,3).*Gd_in.* exp(1i.*kr(:,3,2).*source_depth).*1i.*kr(:,3,load_direction).*p0_plus(:,load_direction,3) + ...
           Gd_in.*exp(1i.*ki(:,in_wave,2).*source_depth).*1i.*ki(:,in_wave,load_direction).*p0_minus(:,load_direction,in_wave);


    case 'point_force'
        u_n = xr(:,1).*p0_plus(:,load_direction,1).*Gd_in.* exp(1i.*kr(:,1,2).*source_depth) + ...
           xr(:,2).*p0_plus(:,load_direction,2).*Gd_in.* exp(1i.*kr(:,2,2).*source_depth) + ...
           xr(:,3).*p0_plus(:,load_direction,3).*Gd_in.* exp(1i.*kr(:,3,2).*source_depth) +...
           Gd_in.*p0_minus(:,load_direction,in_wave).*exp(1i.*ki(:,in_wave,2).*source_depth);
end
u_n(1) = 0;
u_n(size_a)=0;

%% source width term
u_n = u_n .* sinc(sqrt(ki(:,in_wave,1).^2+ki(:,in_wave,2).^2+ki(:,in_wave,3).^2)*source_width.*sin(ray_directions));