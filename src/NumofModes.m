function [nmodes,kr,v] = NumofModes(freq,kr,v,cpmax)

    % cp = 2 * pi * freq ./real(kr);
    % nmodes = 0;
    % for i = 1 : length(kr)
    %     if(cp(i) <= cpmax )
    %         nmodes = i;
    %     end
    % end
    % 
    % if(cp(length(kr)) < cpmax)
    %     nmodes = length(kr);
    % end
    % 
    % if(nmodes == 0)
    %     error('Incorrect maximum phase speed input!');
    % end
    % kr = kr(1 : nmodes);
    % v  = v(:,1 : nmodes);

    %--------example1----------------------------------------------------
    ind    = find(real(kr)<1.85 & real(kr)>1.6 & imag(kr)>0 & imag(kr)<0.005);
    nmodes = length(ind);
    kr     = kr(ind);
    v      = v(:,ind);

end