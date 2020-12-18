function [nmodes,kr,v] = NumOfModes(w,kr,v,cpmax)

    cp = w ./ real(kr);
    nmodes = length( find( cp <= cpmax ) );

    if(nmodes == 0)
        error('Incorrect maximum phase speed input!');
    end
    kr = kr(  1 : nmodes);
    v  = v (:,1 : nmodes);

end
