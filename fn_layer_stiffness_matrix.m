function K = fn_layer_stiffness_matrix(C, rho, layer_thickness, in_wave, n_total, slowness_total, frequency, size_1)
%SUMMARY
%   Find the directivity of ultrasonic waves with different source types
%   under different boundary conditions
%USAGE
%   [K] = fn_layer_stiffness_matrix(C, rho, in_wave, n_total, slowness_total)
%INPUTS
%   C - material stiffness matrix
%   rho - density
%   layer_thickness
%   in_wave
%   n_total
%   slowness_total
%   frequency
%   size_1
%OUTPUTS
%   K - layer stiffness matrix
%AUTHOR
%   Lily Tu (2023)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate -ve polarisation vector - transmitted waves
[p_minus, slowness_minus_xyz, ~, ~] = fn_transmitted_polarisation_vector(C, rho, in_wave, n_total, slowness_total);

C = fn_voigt_to_tensor(C);
% calculate D_minus
D_minus = zeros(size_1,3,3);
for i = 1:3 % stress direction
    for n = 1:3
        for l = 1:3
            for j = 1:3 % wave modes
                D_minus(:,i,j) = D_minus(:,i,j)+C(2,i,n,l)*p_minus(:,l,j)*2*pi*frequency.*slowness_minus_xyz(:,n,j);
            end
        end
    end
end

% calculate +ve polarisation vector
%[p_plus, slowness_plus_xyz, n_plus, ~] = fn_reflected_polarisation_vector(C_old, rho, 1, n_plus, slowness_plus_total);
slowness_plus_xyz = slowness_minus_xyz; % angles, wave modes, stress direction
slowness_plus_xyz(:,2,:) = -1*slowness_plus_xyz(:,2,:);
p_plus = p_minus;
p_plus(:,2,:) = -1*p_plus(:,2,:);

% calculate D_plus
D_plus = zeros(size_1,3,3);
for i = 1:3 % stress direction
    for n = 1:3
        for l = 1:3
            for j = 1:3 % wave modes
                D_plus(:,i,j) = D_plus(:,i,j)+C(2,i,n,l)*p_plus(:,l,j)*2*pi*frequency.*slowness_plus_xyz(:,n,j);
            end
        end
    end
end

% calculate H-matrix
H1 = zeros(size_1,3,3);
for i = 1:size_1
    H1(i,:,:) = diag([exp(-1i.*2*pi*frequency.*slowness_minus_xyz(i,2,1)*layer_thickness),exp(-1i.*2*pi*frequency.*slowness_minus_xyz(i,2,2)*layer_thickness),exp(-1i.*2*pi*frequency.*slowness_minus_xyz(i,2,3)*layer_thickness)]);
end

% calculate K
K = zeros(size_1,6,6);
for i = 1:size_1
K(i,:,:) = [squeeze(D_minus(i,:,:)), squeeze(D_plus(i,:,:))*squeeze(H1(i,:,:)); squeeze(D_minus(i,:,:))*squeeze(H1(i,:,:)), squeeze(D_plus(i,:,:))]...
        /([squeeze(p_minus(i,:,:)), squeeze(p_plus(i,:,:))*squeeze(H1(i,:,:)); squeeze(p_minus(i,:,:))*squeeze(H1(i,:,:)), squeeze(p_plus(i,:,:))]);
end