clc
clear
fid=fopen('tl.bin','rb');
 
nz=fread(fid,1,'int32');
nr=fread(fid,1,'int32');  
tlmin=fread(fid,1,'double');
tlmax=fread(fid,1,'double');
z = fread(fid, [nz, 1], 'double');
r = fread(fid, [1, nr], 'double');
tl= fread(fid, [nz, nr],'double');

fclose(fid);
figure;
pcolor(r./1000, z, tl ); view( 0, 90 );
caxis( [ tlmin tlmax ] ); colormap( hot );
shading flat; colorbar( 'YDir', 'Reverse' ); 
xlabel( 'Range (km)' ); ylabel( 'Depth (m)' );
set(gca,'FontSize',14);axis([0 10 0 900]);
