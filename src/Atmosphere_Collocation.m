clear;
% close all;
clc;
% edit 'input_CAA.txt';
tic
[casename,N,cpmax,dr,zs,rmax,freq,H,tlmin,tlmax,...
          alpha,dep,c] = ReadEnvParameter('input_CAA.txt');

w  = 2 * pi * freq;
nr = rmax / dr;
r  = dr : dr : rmax;

x  = - cos((0 : N) * pi / N)'; %[-1, 1]
z  = (1 + x) * H / 2.0;        %[0 , H]
zr = 0 : 2 : H;

c      = interp1(dep, c   ,z,'linear');
alpha  = interp1(dep,alpha,z,'linear');

k  = w ./ c .* (1.0 + 1i * alpha / (40 * pi * log10(exp(1.0))));

[kr,eigvector] = EigenValueVectorCollocation(N,H,k,x);

% ShowWavenumbers(kr,casename);

[nmodes,kr,eigvector] = NumofModes(freq,kr,eigvector,cpmax);

[psi,psizs] = Normalization(eigvector,x,z,zs);

tl = SynthesizeSoundField(r,kr,psizs,psi);

tl = interp1(z,tl,zr,'linear');

ShowSoundField(r,zr,tl,tlmin,tlmax,casename);

toc
%--------------------------------------------------------------------------

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

function [psi,psizs] = Normalization(eigvector,x,z,zs)

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

    psi   = eigvector * diag(1.0 ./ sqrt(z(length(z)) / 2 * norm));
    psizs = interp1(z,psi,zs,'linear');

end
