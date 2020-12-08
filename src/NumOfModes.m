function [nmodes,kr,v] = NumOfModes(w,kr,v,cpmax)

%     cp = w ./ real(kr);
%     nmodes = length( find( cp <= cpmax ) );
% 
%     if(nmodes == 0)
%         error('Incorrect maximum phase speed input!');
%     end
%     kr = kr(  1 : nmodes);
%     v  = v (:,1 : nmodes);

    %--------example1----------------------------------------------------
    ind    = find(real(kr)<1.85 & real(kr)>1.6 & imag(kr)>0 & imag(kr)<0.005);
    nmodes = length(ind);
    kr     = kr(ind);
    v      = v(:,ind);

end
