function ShowTLcurve(r,zr,tl_zr)

    figure;
    disp('plot the transmission loss curve at zr!');
    plot(r./1000,tl_zr,'r.-','LineWidth',1.5);
    set(gca,'YDir','reverse');
    xlabel( 'Range (km)'); ylabel('TL (dB)');
    title(['Depth=',num2str(zr),'m']);
    set(gca,'FontSize',14);

end