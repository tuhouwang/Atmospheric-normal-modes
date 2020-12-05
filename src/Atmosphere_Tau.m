%点取得太少就容易出问题

clear;
% close all;
clc;
% edit 'input_CAA.txt';
tic
[casename,N,cpmax,dr,zs,zr,dz,rmax,freq,H,...
    tlmin,tlmax,alpha,dep,c] = ReadEnvParameter('input_CAA.txt');

w  = 2 * pi * freq;
nr = rmax / dr;
r  = dr : dr : rmax;

x  = cos((0 : N) * pi / N)';
z  = (1 - x) .* H / 2.0;

zl = 0 : dz : H;
xl = 1 - 2 * zl ./ H;

c      = interp1(dep,c,z,'linear');
alpha  = interp1(dep,alpha,z,'linear');

k  = w ./ c .*(1.0 + 1i * alpha / (40 * pi * log10(exp(1.0))));

[kr,eigvector] = EigenValueVectorTau(N,H,k);

% ShowWavenumbers(kr,casename);

[nmodes,kr,eigvector] = NumofModes(freq,kr,eigvector,cpmax);

[psi,psizs] = Normalization(eigvector,nmodes,H,xl,zl,zs);

[tl, tl_zr] = SynthesizeSoundField(r,kr,zl,zr,psizs,psi);

ShowSoundField(r,zl,tl,tlmin,tlmax,casename);

ShowTLcurve(r,zr,tl_zr);

toc
%-------------------------------------------------------------------------- 

function D  = DerivationMatrix(n)

    D = zeros(n, n);
    for k = 1 : n
        j = k+1 : 2 : n;
        D(k, j) = 2 * j - 2;
    end
    D(1, :) = D(1, :) / 2;

end

function C  = ConvolutionMatrix(v)

    n = length(v);
    C = zeros(n, n);
    
    for i = 1 : n
        for k = 1 : n
            
            j = k - i + 1;       
            if (j >= 1 && j <= n)
                C(k,i) = C(k,i) + v(j);
            end
            
            j = i - k + 1;
            if (j <= n && j >= 1)
                C(k,i) = C(k,i) + v(j);
            end
            
            if (k > 1)
                j = i + k - 1 ;
                if (j <= n && j >= 1)
                    C(k,i) = C(k,i) + v(j);
                end  
            end
            C(k,i) = C(k,i) * 0.5; 
        end
    end

end

function fk = ChebTransFFT(N, fx)

%  The function computes the Chebyshev tranforms by FFT in the nodes:
%  x_j = cos(j pi/N) from Physical space to Spectral space.
%  
%  Input:
%  N:           Degree of polynomials, must be even
%  fx:          vector to be transformed

%  Output:
%  fk:          Chebyshev coefficients of fx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fk = fft([fx; flipud(fx(2:N))]);              % Extend and compute fft
    fk = fk(1:N+1).*[0.5; ones(N-1,1); 0.5]/N;    % fk contains Chebyshev
                                                  % coefficients of fx
end

function fx = InvChebTrans(fk,z)

    n = size(fk, 1);
    m = length(z);
    T = zeros(m, n);
    
    for k = 0 : n-1
        T(:, k+1) = cos( k * acos(z) );
    end
    
    fx = T * fk;

end

function [kr,eigvector] = EigenValueVectorTau(N,H,k)
    D  = DerivationMatrix(N+1);
    %----------------------------------------------------
    A =4.0 / H ^2 * D * D + ConvolutionMatrix( ChebTransFFT(N, k.^2) );

    Pg=-2/H * ones(1, N+1) * D + 1i * k(1) / (12.97 + 12.38i) * ones(1, N+1);
    Pa=-2/H * ((-1.0).^(0 : N)) * D - 1i * k(N+1) * (-1.0).^(0 : N);

    A(N,:)   = Pg;
    A(N+1,:) = Pa;
    %矩阵分块
    L11  = A(1:N-1,1:N-1);
    L12  = A(1:N-1,N:N+1);
    L21  = A(N:N+1,1:N-1);
    L22  = A(N:N+1,N:N+1);
    L    = L11 - L12 * (L22 \ L21);

    [v, kr] = eig(L);
    v2 = - (L22 \ L21) * v;

    eigvector = [v; v2];

    kr = sqrt(diag(kr));
    [~,ind] = sort(real(kr),'descend');
    kr = kr(ind);
    eigvector = eigvector(:,ind);

end

function [psi,psizs] = Normalization(eigvector,nmodes,H,xr,zr,zs)
   
   psi   = InvChebTrans(eigvector, xr); 
   N     = size(eigvector, 1) - 1;
   f     = zeros(N+1, nmodes);
   
   P = zeros(1, N+1);
   k = 0 : 2 : N;
   P(k+1) = -2 ./(k.^2-1);
   
   for j=1 : nmodes
      f(:,j) = ConvolutionMatrix(eigvector(:, j)) * eigvector(:, j);  
   end
   
   psi   = psi * diag(1.0 ./sqrt(P * f * H / 2));  
   psizs = interp1(zr, psi, zs, 'linear');
   
end
