function ShowWavenumbers(kr,casename)

%     figure;
    disp('plot the modal wavenumbers!');
    plot(real(kr),imag(kr),'bo');grid on;
    xlabel('Real Wave Number (1/m)');
    ylabel('Imaginary Wave Number (1/m)');
    title(casename);set(gca,'FontSize',14);
%     axis( [1.6 1.85 0 0.035] );
end