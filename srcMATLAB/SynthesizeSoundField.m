function [tl, tl_zr] = SynthesizeSoundField(r,kr,z,zr,psizs,psi)

    range = kr * r;
    range = exp(1i * range) ./ sqrt(range);
    p     = psi * diag(psizs) * range * sqrt(2 * pi);
    tl    = -20 * log10(abs(p));
    tl_zr = interp1(z,tl,zr,'linear');

end
