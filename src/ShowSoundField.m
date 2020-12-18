function ShowSoundField(r,z,tl,tlmin,tlmax,casename)

    figure;
    disp('plot the transmission loss field');
    pcolor( r./1000, z, tl );
    caxis( [tlmin tlmax] ); colormap( hot );
    shading flat; view( 0, 90 );box off;
    xlabel( 'Range (km)'); ylabel( 'Height (m)');
    colorbar('YDir', 'Reverse');title(casename);
    set(gca,'FontSize',14);

end