function [p_transmit, slowness_t, n_all, s_in] = fn_transmitted_polarisation_vector(C, rho, in_wave, n_total, slowness_total)
%SUMMARY
%   Find the reflected polarisation vector given a specific incident wave
%USAGE
%   [p_transmit, slowness_t, n_all, s_in] = fn_transmitted_polarisation_vector(C, rho, in_wave, n_total, slowness_total)
%INPUTS
%   C - 6x6 stiffness tensor of the transmitted media
%   rho - density of the transmitted media
%   in_wave - incident wave mode
%   n_total
%   slowness_total
%OUTPUTS
%   p_transmit - transmitted polarisation vector
%   slowness_t
%   n_all
%   s_in
%AUTHOR
%   Lily Tu (2023)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
size_1 = size(slowness_total,1);

% 3D slowness curve
% rotated around y-axis
slowness_x = slowness_total.*n_total(:,1); 
slowness_y = slowness_total.*n_total(:,2);
%slowness_z = slowness_total.*n_total(:,3);

%% find transmitted slowness vector at the interface
s_in = [slowness_x(:,in_wave),slowness_y(:,in_wave),...
   zeros(size_1,1)]';
% S1 known, S3=0, find S2

for i = 1:size_1
% p_1*x^n+p_2*x^(n-1)+...+p_n+1=0
    % D, B, A
A = [C(1,1)*s_in(1,i)^2-rho, C(1,6)*s_in(1,i)^2, C(1,5)*s_in(1,i)^2;...
     C(1,6)*s_in(1,i)^2, C(6,6)*s_in(1,i)^2-rho, C(6,5)*s_in(1,i)^2;...
     C(1,5)*s_in(1,i)^2, C(6,5)*s_in(1,i)^2, C(5,5)*s_in(1,i)^2-rho];

B = [2*C(1,6)*s_in(1,i), (C(1,2)+C(6,6))*s_in(1,i), (C(1,4)+C(6,5))*s_in(1,i);...
     (C(1,2)+C(6,6))*s_in(1,i), 2*C(6,2)*s_in(1,i), (C(6,4)+C(2,5))*s_in(1,i);...
     (C(1,4)+C(6,5))*s_in(1,i), (C(6,4)+C(2,5))*s_in(1,i), 2*C(5,4)*s_in(1,i)];

D = [C(6,6),C(6,2),C(6,4);...
     C(6,2),C(2,2),C(2,4);...
     C(6,4),C(2,4),C(4,4)];

p = [det([D(:,1),D(:,2),D(:,3)])... %beta^6

det([B(:,1),D(:,2),D(:,3)])+det([D(:,1),B(:,2),D(:,3)])+det([D(:,1),D(:,2),B(:,3)]) ... %beta^5

det([A(:,1),D(:,2),D(:,3)])+det([B(:,1),B(:,2),D(:,3)])+det([B(:,1),D(:,2),B(:,3)])+...
    det([D(:,1),A(:,2),D(:,3)])+det([D(:,1),B(:,2),B(:,3)])+det([D(:,1),D(:,2),A(:,3)]) ... %beta^4

det([A(:,1),B(:,2),D(:,3)])+det([A(:,1),D(:,2),B(:,3)])+det([B(:,1),A(:,2),D(:,3)])+...
    det([B(:,1),B(:,2),B(:,3)])+det([B(:,1),D(:,2),A(:,3)])+det([D(:,1),A(:,2),B(:,3)])+...
    det([D(:,1),B(:,2),A(:,3)]) ... %beta^3

det([A(:,1),A(:,2),D(:,3)])+det([A(:,1),B(:,2),B(:,3)])+det([A(:,1),D(:,2),A(:,3)])+...
    det([B(:,1),A(:,2),B(:,3)])+det([B(:,1),B(:,2),A(:,3)])+det([D(:,1),A(:,2),A(:,3)]) ... %beta^2

det([A(:,1),A(:,2),B(:,3)])+det([A(:,1),B(:,2),A(:,3)])+det([B(:,1),A(:,2),A(:,3)]) ... % beta^1

det([A(:,1),A(:,2),A(:,3)])]; %beta^0

t(:,i) = roots(p);

transmitted_slowness=zeros(6,1);
for j = 1:6
if isreal(t(j,i))
transmitted_slowness(j)=(t(j,i)>0)*t(j,i);
else 
    transmitted_slowness(j)=(imag(t(j,i))>0)*t(j,i);
end
end
transmitted_slowness(transmitted_slowness==0)=[];


if length(transmitted_slowness)<3
    transmitted_slowness = [transmitted_slowness;zeros((3-length(transmitted_slowness)),1)];
end
transmitted_slowness = sort(transmitted_slowness,'ComparisonMethod','real');
s_tr(:,i) = transmitted_slowness(1:3);
if ~isreal(s_tr(2,i)) && ~isreal(s_tr(3,i))
    s_tr(2:3,i) = sort(s_tr(2:3,i));
end
end

% for x-z plane
slowness_t = zeros(size_1,3,3); 
% no of planes, angles in plane, 3 axes, 3 modes 
% in the order 1. slow qS, 2. fast qS, 3. qL
% 1. x, 2. y, 3. z
slowness_t(:,:,3) = [slowness_x(:,in_wave),...
    s_tr(1,:)', zeros(size_1,1)];
slowness_t(:,:,2) = [slowness_x(:,in_wave), ...
    s_tr(2,:)', zeros(size_1,1)];
slowness_t(:,:,1) = [slowness_x(:,in_wave), ...
    s_tr(3,:)', zeros(size_1,1)];

%% transmitted polarisation vector from slowness vector at the interface
n_all = slowness_t(:,:,:);
p_transmit = zeros(size(n_all));
n_all = n_all ./ sqrt(sum(n_all .^ 2, 2));
C = fn_voigt_to_tensor(C);
for mm = 1:3
    n = n_all(:,:,mm);
for nn = 1:size_1 %loop over phase velocity directions
    %Prepare Christoffel equation matrix
    Cnn = zeros(3);
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    Cnn(i,k) = Cnn(i,k) + C(i,j,k,l) * n(nn,j) * n(nn,l);
                end
            end
        end
    end
    %Solve Eigenvalue problem to get phase velocity magnitude
    [V, D] = eig(Cnn);
    D = diag(D);
    
    %sort in ascending order to keep mode order consistent
    [~, i] = sort(D);
    D = D(i);
    V = V(:, i);
    
    %Create polarisation vector from Eigenvectors
    %%% new edits - make the polarisation vector pointing away from the
    %%% origin
    if dot(n(nn,:), V(:,mm)) >= 0
        p_transmit(nn, :, mm) = V(:, mm);
    else
        p_transmit(nn, :, mm) = -V(:, mm);
    end
    %%%
end
end