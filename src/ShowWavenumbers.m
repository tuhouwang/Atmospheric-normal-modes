function ShowWavenumbers(kr,casename)

    figure;
    disp('plot the modal wavenumbers!');
    plot(real(kr),imag(kr),'r*');grid on;
    axis( [1.6 1.85 0 0.035] );
    xlabel('Real Wave Number (1/m)');
    ylabel('Imaginary Wave Number (1/m)');
    title(casename);
    set(gca,'FontSize',14,'FontName','Times New Roman');

end