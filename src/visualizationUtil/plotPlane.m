function plotPlane(a, x_mink, color)
[x, y]= meshgrid(-1.5:0.1:1.5);

zPlane = -a(1)/a(3) * x - a(2)/a(3)*y + a(1)/a(3)*x_mink(1) + a(2)/a(3)*x_mink(2) + x_mink(3);

scatter3(x_mink(1), x_mink(2), x_mink(3), 'MarkerFaceColor', color);

surf(x, y, zPlane, 'FaceColor', color, 'FaceAlpha', 0.5);
end