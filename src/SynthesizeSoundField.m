function tl = SynthesizeSoundField(r,kr,psizs,psi)

    range = kr * r;
    range = exp(1i * range) ./ sqrt(range);
    p     = psi * diag(psizs) * range * sqrt(2 * pi);
    tl    = -20 * log10(abs(p));

end