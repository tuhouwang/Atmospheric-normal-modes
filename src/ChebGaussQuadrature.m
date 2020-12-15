function integral = ChebGaussQuadrature(f)

    integral = 0;
    n  = length(f);
    w  = [pi / 2 / (n - 1); pi / (n - 1) * ones(n - 2, 1); pi / 2 / (n - 1)];
    x  = cos( (0 : n-1) * pi / (n - 1) )';
    for k = 1 : n
         integral = integral + w(k) * f(k) * ( sqrt( 1 - x(k) ^ 2 ) );
    end

end
