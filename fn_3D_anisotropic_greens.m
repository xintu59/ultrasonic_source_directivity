function [Gd_in, p_total, ray_directions, slowness_total, n_total, wave_number] = fn_3D_anisotropic_greens(a, C, rho, f, in_wave, d)
%SUMMARY
%   Find the wave displacement generated by an arbitrary force using 3D
%   Green's function in an anisotropic medium
%USAGE
%   [Gd_in, p_in, ray_directions, slowness_in, n_in, wave_number] = fn_3D_anisotropic_greens(a, C, rho, f, in_wave, d)
%INPUTS
%   a - angles
%   C - 6x6 stiffness tensor
%   rho - density
%   f - frequency
%   in_wave - wave mode
%   d - observation distance
%OUTPUTS
%   Gd_in - wave amplitude
%   p_in - polarisation vector
%   ray_directions - group velocity angles
%   slowness_in
%   n_in
%   wave_number
%AUTHOR
%   Lily Tu (2023)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
size_a = length(a);
da = a(2)-a(1);
b = (0:da:2*pi);
slice_number = floor(length(b)/2)+1;
% a is in the range [0, pi]
%% compute phase and group velocities
slowness_total = zeros(5,size_a,1,3);
n_total = zeros(5,size_a,1,3);
for i = 1:5
    idx = slice_number -3 + i;
    n = [cos(a)*cos(b(idx)), sin(a), ones(size_a,1).*cos(a)*sin(b(idx))];
    [vph, vgr, p] = fn_anisotropic_vel_profile(C, rho, n);

    phase_velocity = sqrt(vph(:,1,:).^2+vph(:,2,:).^2+vph(:,3,:).^2);
    slowness = 1./ phase_velocity;
    group_velocity = sqrt(vgr(:,1,:).^2+vgr(:,2,:).^2+vgr(:,3,:).^2);

    % matrix dimensions
    % no. of curves x no. of points on the curve x 3 axis x 3 wave modes
    slowness_total(i,:,:,:)=slowness;
    n_total(i,:,:,:)=n;
    if i == 3
        vgr_total=vgr;
        p_total=p;
        group_velocity_total=squeeze(group_velocity);
    else
    end
end
% 3D slowness curve
% rotated around y-axis
slowness_x = slowness_total.*n_total(:,:,:,1); 
slowness_y = slowness_total.*n_total(:,:,:,2);
slowness_z = slowness_total.*n_total(:,:,:,3);

%% incident wave amplitude - Media 1 (substrate)
% gaussian curvature of slowness surface

K_total=zeros(5,size_a,3);
H_total=zeros(5,size_a,3);
Pmax_total=zeros(5,size_a,3);
Pmin_total=zeros(5,size_a,3);

for m =1:3
[K,H,Pmax,Pmin] = surfature(squeeze(slowness_x(:,:,:,m)),squeeze(slowness_y(:,:,:,m)),squeeze(slowness_z(:,:,:,m))); % gaussian curvature of 3D surface
K(:,floor(size_a/2)+1)=K(:,floor(size_a/2)); % at these points K is erratic
H(:,floor(size_a/2)+1)=H(:,floor(size_a/2)); % at these points K is erratic
K(:,1)=K(:,3);
K(:,2)=K(:,3);
K(:,end-1)=K(:,end-2);
K(:,end)=K(:,end-2);
H(:,1)=H(:,3);
H(:,2)=H(:,3);
H(:,end-1)=H(:,end-2);
H(:,end)=H(:,end-2);
K_total(:,:,m)=K;
H_total(:,:,m)=H;
Pmax_total(:,:,m)=Pmax;
Pmin_total(:,:,m)=Pmin;
end

% wavenumber
wave_number = 2*pi*f.*squeeze(slowness_total(3,:,1,:)); % angular frequency / phase velocity
sig_0 = 1 + 0.5*sign(Pmin_total(1,:,in_wave)) + 0.5*sign(Pmax_total(1,:,in_wave));
p_in = dot(p_total(:,:,in_wave),p_total(:,:,in_wave),2);
Gd_in = p_in'./(4*pi*rho*sqrt(K_total(1,:,in_wave)).*...
group_velocity_total(:,in_wave)'*d*100).*exp(1i*(wave_number(:,in_wave)'*d+0.5*pi*sig_0)); % -w.*t %
Gd_in = Gd_in';

% p_in_vector = zeros(3,3,size_a);
% for i = 1:3
%     for n = 1:3
%         p_in_vector(i,n,:) =  p_total(:,i,in_wave).*p_total(:,n,in_wave);
%         Gd_in_vector(:,i,n) = squeeze(p_in_vector(i,n,:))'./(4*pi*rho*sqrt(abs(K_total(3,:,in_wave))).*...
%     group_velocity_total(:,in_wave)'*d).*exp(1i*(wave_number(:,in_wave)'*d+0.5*pi*sig_0)); % -w.*t %
%     end
% end
% 
% Gd_in=0;
% for i=1:3
%     for n=1:3
%     Gd_in = Gd_in + (abs(Gd_in_vector(:,i,n))).^2;
%     end
% end
% Gd_in = sqrt(Gd_in);

d1 = sqrt(vgr_total(:,1,in_wave).^2+vgr_total(:,3,in_wave).^2); %%%
ray_directions = atan2(abs(d1),abs(vgr_total(:,2,in_wave)))*(-1).*sign(vgr_total(:,1,in_wave));

slowness_total = squeeze(slowness_total(3,:,1,:));
n_total = squeeze(n_total(3,:,:,:));
return