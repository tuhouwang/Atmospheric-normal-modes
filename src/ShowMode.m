function ShowMode(psi, z)

    mode_num = input('What mode number do you want to plot?:');

%     figure;
    plot(real(psi(:, mode_num)), z, 'k-', 'LineWidth', 2) ;
    hold on;
    plot(imag(psi(:, mode_num)), z, 'y--', 'LineWidth', 2);

    ylabel( 'Height (m)');
    set(gca, 'FontSize', 20);

end
