function ShowSoundField(r,z,tl,tlmin,tlmax,casename)

    figure;
    disp('plot the transmission loss field');
    pcolor( r./1000, z, tl );
    caxis( [tlmin tlmax] ); colormap( gray );
    shading flat; view( 0, 90 );
    xlabel( 'Range (km)'); ylabel( 'Height (m)');
    colorbar;title(casename);box off;
    set(gca,'FontSize',14);

end