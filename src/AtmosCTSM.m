
clear;
% close all;
clc;
% edit 'input_CAA.txt';
tic
[casename, N, cpmax, dr, zs, zr, dz, rmax, freq, H, Zground, tlmin, ...
         tlmax, alpha, dep, c] = ReadEnvParameter('input_CAA.txt');

w  = 2 * pi * freq;
nr = rmax / dr;
r  = dr : dr : rmax;

x  = cos((0 : N) * pi / N)';
z  = (1 - x) .* H / 2.0;

zl = 0 : dz : H;
xl = 1 - 2 * zl ./ H;

c      = interp1(dep, c, z, 'linear');
alpha  = interp1(dep, alpha, z, 'linear');

k  = w ./ c .*(1.0 + 1i * alpha / (40 * pi * log10(exp(1.0))));

[kr, eigvector] = EigenValueVectorTau(N, H, Zground, k);

% ShowWavenumbers(kr, casename);

[nmodes, kr, eigvector] = NumOfModes(w, kr, eigvector, cpmax);

[psi, psizs] = Normalization(eigvector, nmodes, H, x, xl, zl, zs);

% ShowMode(psi, zl);

[tl,  tl_zr] = SynthesizeSoundField(r, kr, zl, zr, psizs, psi);

% ShowSoundField(r, zl, tl, tlmin, tlmax, casename);

ShowTLcurve(r, zr, tl_zr);

toc
%-------------------------------------------------------------------------- 

function [kr, eigvector] = EigenValueVectorTau(N, H, Zground, k)
    D  = DerivationMatrix(N+1);
    %----------------------------------------------------
    A =4.0 / H ^2 * D * D + ConvolutionMatrix( ChebTransFFT(N, k.^2) );

    Pg=-2 / H * ones(1, N+1) * D + 1i * k(1) / Zground * ones(1, N+1);
    Pa=-2 / H * ((-1.0) .^ (0 : N)) * D - 1i * k(N + 1) * (-1.0).^(0 : N);

    A(N,   :) = Pg;
    A(N+1, :) = Pa;
    %blocking
    L11  = A(1:N-1, 1:N-1);
    L12  = A(1:N-1, N:N+1);
    L21  = A(N:N+1, 1:N-1);
    L22  = A(N:N+1, N:N+1);
    L    = L11 - L12 * (L22 \ L21);

    [v, kr] = eig(L);
    v2 = - (L22 \ L21) * v;

    eigvector = [v; v2];

    kr = sqrt(diag(kr));
    [~, ind] = sort(real(kr), 'descend');
    kr = kr(ind);
    eigvector = eigvector(:, ind);

end

function [psi, psizs] = Normalization(eigvector, nmodes, H, x, xr, zr, zs)
   
   psix  = InvChebTrans(eigvector, x); 
   norm  = zeros(nmodes,1);   

     for j=1 : nmodes
          norm(j) = ChebGaussQuadrature( psix(:, j) .^2 ) * H / 2;  
     end
  
   psi   = InvChebTrans(eigvector, xr);      
   psi   = psi * diag(1.0 ./ sqrt( norm ) ); 
   psizs = interp1(zr, psi, zs, 'linear');

%    psi   = InvChebTrans(eigvector, xr); 
%    N     = size(eigvector, 1) - 1;
%    f     = zeros(N+1, nmodes);
%    
%    P = zeros(1, N+1);
%    k = 0 : 2 : N;
%    P(k+1) = -2 ./(k.^2-1);
%    
%    parfor j=1 : nmodes
%       f(:,j) = ConvolutionMatrix(eigvector(:, j)) * eigvector(:, j);  
%    end
%    
%    psi   = psi * diag(1.0 ./sqrt(P * f * H / 2));  
%    psizs = interp1(zr, psi, zs, 'linear');

end
