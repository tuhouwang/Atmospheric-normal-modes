function ShowMode(psi, z)

    mode_num = input('What mode number do you want to plot?:');

    figure;
    plot(real(psi(1:381, mode_num)), z(1:381), 'b-', 'LineWidth', 1) ;
    hold on;
    plot(imag(psi(1:381, mode_num)), z(1:381), 'c--', 'LineWidth', 1);

    ylabel( 'Height (m)');
    set(gca, 'FontSize', 20);

end
