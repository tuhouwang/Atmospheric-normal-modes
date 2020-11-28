clear;
close all;
clc;
% edit 'input_CAA.txt';
tic
[casename,N,cpmax,dr,zs,rmax,freq,H,tlmin,tlmax,...
          alpha,dep,c] = Read_input('input_CAA.txt');

w  = 2 * pi * freq;
nr = rmax / dr;
r  = dr : dr : rmax;

x  = - cos((0 : N) * pi / N)'; %[-1, 1]
z  = (1 + x) * H / 2.0;        %[0 , H]

c      = interp1(dep, c   ,z,'linear');
alpha  = interp1(dep,alpha,z,'linear');

k  = w ./ c .* (1.0 + 1i * alpha / (40 * pi * log10(exp(1.0))));

[kr,eigvector] = EigenValueVectorCollocation(N,H,k,x);

% PrintWavenumbers(kr,casename);

[nmodes,kr,eigvector] = NumofModes(freq,kr,eigvector,cpmax);

[psi,psizs] = NormCollocation(eigvector,x,z,zs);

tl = SynthesizeSoundField(r,kr,psizs,psi);

ShowSoundField(r,z,tl,tlmin,tlmax,casename);

toc
%--------------------------------------------------------------------------

function [casename,N,cpmax,dr,zs,rmax,freq,H,...
    tlmin,tlmax,alpha,dep,c]=Read_input(env_file)
fid      = fopen(env_file);
casename = fgetl(fid);
N        = fscanf(fid,'%d',1);
cpmax    = fscanf(fid,'%f',1);
freq     = fscanf(fid,'%f',1);
zs       = fscanf(fid,'%f',1);
rmax     = fscanf(fid,'%f',1);
dr       = fscanf(fid,'%f',1);
H        = fscanf(fid,'%f',1);
tlmin    = fscanf(fid,'%f',1);
tlmax    = fscanf(fid,'%f',1);
n        = fscanf(fid,'%d',1);
if(H >0.0  &&  N>2)
    Profile    = fscanf(fid,'%f %f',[3,n]);
    dep(1:n)   = Profile(1,1:n);
    c(1:n)     = Profile(2,1:n);
    alpha(1:n) = Profile(3,1:n);
else
    error('Error! H must greater than 0!');
end
% Check the input sound profile
if(dep(1) ~=0.0 || dep(n)~=H)
    error('Error! input sound profile is unsuitable!');
end

if((rmax / dr - floor(rmax / dr))~=0)
    error('Please reinput the dr and rmax!');
end

if(tlmin >= tlmax)
    error('tlmin must less than tlmax!');
end

fclose(fid);
end

function [kr,eigvector] = EigenValueVectorCollocation(N,H,k,x)

D          = Collocation(N,x);
A          = 4.0 / H /H * D * D + diag(k.^2);
%Set boundary conditions
A(1,:)     = 2.0 / H * D(1,:);
A(1,1)     = A(1,1) + 1i * k(1) / (12.97 + 12.38i);
A(N+1,:)   = 2.0 / H * D(N+1,:);
A(N+1,N+1) = A(N+1,N+1) - 1i * k(N+1);
% Row transformations and column transformations
A          = [A(2:N,:);A(1,:);A(N+1,:)];
A          = [A(:,2:N),A(:,1),A(:,N+1)];

L11        = A(1:N-1,1:N-1);
L12        = A(1:N-1,N:N+1);
L21        = A(N:N+1,1:N-1);
L22        = A(N:N+1,N:N+1);
L          = L11 - L12 * (L22 \ L21);

[v, kr]    = eig(L);
v2         = - (L22 \ L21) * v;

eigvector  = [v2(1,:); v ; v2(2,:)];
kr         = sqrt(diag(kr));
[~,ind]    = sort(real(kr),'descend');
kr         = kr(ind);
eigvector  = eigvector(:,ind);

end

function D = Collocation(N,x)

c  = [2; ones(N-1,1); 2] .* (-1).^(0 : N)';

X  = repmat(x,1,N+1);
dX = X - X';

D = (c * (1 ./ c)') ./ (dX + eye(N + 1));
D = D - diag(sum(D, 2));

end

function [nmodes,kr,v] = NumofModes(freq,kr,v,cpmax)
% cp = 2 * pi * freq ./real(kr);
% nmodes = 0;
% for i = 1 : length(kr)
%     if(cp(i) <= cpmax )
%         nmodes = i;
%     end
% end
% 
% if(cp(length(kr)) < cpmax)
%     nmodes = length(kr);
% end
% 
% if(nmodes == 0)
%     error('Incorrect maximum phase speed input!');
% end
% kr = kr(1 : nmodes);
% v  = v(:,1 : nmodes);

%--------example1----------------------------------------------------
ind    = find(real(kr)<1.85 & real(kr)>1.6 & imag(kr)>0 & imag(kr)<0.005);
nmodes = length(ind);
kr     = kr(ind);
v      = v(:,ind);

end

function [psi,psizs] = NormCollocation(eigvector,x,z,zs)

N    = length(z) - 1;
w    = [pi / 2 / N; pi / N * ones(N-1,1); pi/ 2 / N];
norm = zeros(size(eigvector, 2), 1);

for i = 1 : size(eigvector, 2)

    f = eigvector(:,i);
    f = f.^2;
    
    for j = 1 : length(z)
         norm(i) = norm(i) + w(j) * f(j) * (sqrt(1 - x(j) ^2));
    end
     
end

    psi   = eigvector * diag(1.0 ./ sqrt(z(length(z)) * norm));
    psizs = interp1(z,psi,zs,'linear');

end

function tl = SynthesizeSoundField(r,kr,psizs,psi)

range = kr * r;
range = exp(1i * range) ./ sqrt(range);

p     = psi * diag(psizs) * range * sqrt(2 * pi);
tl    = -20 * log10(abs(p));

end

function ShowSoundField(r,z,tl,tlmin,tlmax,casename)

figure;
disp('plot the transmission loss field');
pcolor( r./1000, z, tl );
caxis( [tlmin tlmax] ); colormap( gray );
shading flat; view( 0, 90 );
xlabel( 'Range (km)'); ylabel( 'Depth (m)');
colorbar;title(casename);
set(gca,'FontSize',16,'FontName','Times New Roman');

end

function PrintWavenumbers(kr,casename)

disp('plot the modal wavenumbers!');
figure;
plot(real(kr),imag(kr),'r*');grid on;
axis( [1.6 1.85 0 0.035] );
xlabel('Real Wave Number (1/m)');
ylabel('Imaginary Wave Number (1/m)');
title(casename);
set(gca,'FontSize',14,'FontName','Times New Roman');

end