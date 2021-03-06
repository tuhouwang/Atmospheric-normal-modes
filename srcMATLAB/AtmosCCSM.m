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

x  = - cos((0 : N) * pi / N)'; %[-1,  1]
z  = (1 + x) * H / 2.0;        %[0 ,  H]
zl = 0 : dz : H;

c      = interp1(dep,  c   , z, 'linear');
alpha  = interp1(dep, alpha, z, 'linear');

k  = w ./ c .* (1.0 + 1i * alpha / (40 * pi * log10(exp(1.0))));

[kr, eigvector] = EigenValueVectorCollocation(N, H, Zground, k, x);

ShowWavenumbers(kr, casename);

[nmodes, kr, eigvector] = NumOfModes(w, kr, eigvector, cpmax);

[psi, psizs] = Normalization(eigvector, H, z, zs);

[tl, tl_zr] = SynthesizeSoundField(r, kr, z, zr, psizs, psi);

% ShowMode(psi, z);

tl = interp1(z, tl, zl, 'linear');

ShowSoundField(r, zl, tl, tlmin, tlmax, casename);

% ShowTLcurve(r, zr, tl_zr);

toc
%--------------------------------------------------------------------------

function [kr, eigvector] = EigenValueVectorCollocation(N, H, Zground, k, x)

    D          = Collocation(N, x);
    A          = 4.0 / H /H * D * D + diag(k .^ 2);
    %Set boundary conditions
    A(1,     :) = 2.0 / H * D(1, :);
    A(1,     1) = A(1, 1) + 1i * k(1) / Zground;
    A(N+1,   :) = 2.0 / H * D(N + 1, :);
    A(N+1, N+1) = A(N + 1, N + 1) - 1i * k(N + 1);
    % Row transformations and column transformations
    A          = [A(2:N, :); A(1, :); A(N+1, :)];
    A          = [A(:, 2:N), A(:, 1), A(:, N+1)];

    L11        = A(1:N-1, 1:N-1);
    L12        = A(1:N-1, N:N+1);
    L21        = A(N:N+1, 1:N-1);
    L22        = A(N:N+1, N:N+1);
    L          = L11 - L12 * (L22 \ L21);

    [v, kr]    = eig(L);
    v2         = - (L22 \ L21) * v;

    eigvector  = [v2(1, :); v; v2(2, :)];
    kr         = sqrt(diag(kr));
    [~, ind]   = sort(real(kr), 'descend');
    kr         = kr(ind);
    eigvector  = eigvector(:, ind);

end

function D = Collocation(N, x)

    c  = [2; ones(N-1, 1); 2] .* (-1).^(0 : N)';

    X  = repmat(x, 1, N + 1);
    dX = X - X';

    D = (c * (1 ./ c)') ./ (dX + eye(N + 1));
    D = D - diag(sum(D, 2));

end

function [psi, psizs] = Normalization(eigvector, H, z, zs)

    norm = zeros( size( eigvector, 2 ), 1 );

    for i = 1 : size( eigvector, 2 )

        norm(i) = ChebGaussQuadrature( eigvector(:, i) .^ 2 ) * H / 2 ;

    end

    psi   = eigvector * diag( 1.0 ./ sqrt( norm ) );
    psizs = interp1(z, psi, zs, 'linear');

end
