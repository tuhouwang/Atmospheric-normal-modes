%点取得太少就容易出问题

clear;
close all;
clc;
% edit 'input_CAA.txt';
tic
[casename,N,cpmax,dr,zs,rmax,freq,H,...
    tlmin,tlmax,alpha,dep,c]=Read_input('input_CAA.txt');

w  = 2 * pi * freq;
nr = rmax / dr;
r  = dr : dr : rmax;

x  = cos((0 : N) * pi / N)';
z  = (1 - x) .* H / 2.0;

zr = 0 : 2 : H;
xr = 1 - 2 * zr ./ H;

c      = interp1(dep,c,z,'linear');
alpha  = interp1(dep,alpha,z,'linear');

k  = w ./ c .*(1.0 + 1i * alpha / (40 * pi * log10(exp(1.0))));

[kr,eigvector] = EigenValueVectorTau(N,H,k);

% ShowWavenumbers(kr,casename);

[nmodes,kr,eigvector] = NumofModes(freq,kr,eigvector,cpmax);

[psi,psizs] = Normalization(eigvector,nmodes,H,xr,zr,zs);

tl = SynthesizeSoundField(r,kr,psizs,psi);

ShowSoundField(r,zr,tl,tlmin,tlmax,casename);

toc
%-------------------------------------------------------------------------- 


function [casename,N,cpmax,dr,zs,rmax,freq,H,...
    tlmin,tlmax,alpha,dep,c]=Read_input(env_file)
    fid      = fopen(env_file);
    casename = fgetl(fid);
    N = fscanf(fid,'%d',1);
    cpmax = fscanf(fid,'%f',1);
    freq  = fscanf(fid,'%f',1);
    zs    = fscanf(fid,'%f',1);
    rmax  = fscanf(fid,'%f',1);
    dr    = fscanf(fid,'%f',1);
    H     = fscanf(fid,'%f',1);
    tlmin = fscanf(fid,'%f',1);
    tlmax = fscanf(fid,'%f',1);
    n   = fscanf(fid,'%d',1);
    if(H >0.0  && N>2)
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

    if((rmax/dr-floor(rmax/dr))~=0)
        error('Please reinput the dr and rmax!');
    end

    if(tlmin >= tlmax)
        error('tlmin must less than tlmax!');
    end

    fclose(fid);
end

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

function [nmodes,kr,v]= NumofModes(freq,kr,v,cpmax)
    % cp=2*pi*freq./real(kr);
    % nmodes =0;
    % for i=1:length(kr)
    %     if(cp(i)<=cpmax )
    %         nmodes=i;
    %     end
    % end
    % 
    % if(cp(length(kr))<cpmax)
    %     nmodes=length(kr);
    % end
    % 
    % if(nmodes==0)
    %     error('Incorrect maximum phase speed input!');
    % end
    % kr=kr(1:nmodes);
    % v=v(:,1:nmodes);

    %--------example1----------------------------------------------------
    ind = find(real(kr)<1.85 & real(kr)>1.6 & imag(kr)>0 & imag(kr)<0.005);
    nmodes = length(ind);
    kr = kr(ind);
    v  = v (:, ind);

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

function tl = SynthesizeSoundField(r,kr,psizs,psi)

    range = kr * r;
    range = exp(1i * range) ./ sqrt(range);

    p     = psi * diag(psizs) * range * sqrt(2 * pi);
    tl    = -20 * log10(abs(p));

end

function ShowWavenumbers(kr,casename)

    disp('plot the modal wavenumbers!');
    figure;
    plot(real(kr),imag(kr),'r*');grid on;
    axis( [1.6 1.85 0 0.035] );
    xlabel('Real Wave Number (1/m)');
    ylabel('Imaginary Wave Number (1/m)');
    title(casename);
    set(gca,'FontSize',14,'FontName','Times New Roman');

end

function ShowSoundField(r,z,tl,tlmin,tlmax,casename)

    figure;
    disp('plot the transmission loss field');
    pcolor( r./1000, z, tl );
    caxis( [tlmin tlmax] ); colormap( gray );
    shading flat; view( 0, 90 );
    xlabel( 'Range (km)'); ylabel( 'Depth (m)');
    colorbar;
    title(casename);
    set(gca,'FontSize',16,'FontName','Times New Roman');

end
