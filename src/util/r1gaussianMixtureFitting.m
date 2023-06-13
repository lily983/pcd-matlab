function r1gaussianMixtureFitting
% create a Java slider bar for mb
figure
jSlider = javax.swing.JSlider;
%hjSlider: handle of jSlider. Callback
%functions is attached to it.
%hContainer: set properties of the jSlider position and size relative to
%the figure.
[hjSlider, hContainer] = javacomponent(jSlider, [10, 70, 400, 45]);
set(jSlider, 'MajorTickSpacing', 50, 'PaintLabels', true);
jSlider.setMaximum(300);

set(hContainer, 'position', [10, 20, 500, 45])
set(hjSlider, 'StateChangedCallback', @gaussianMixtureValue);

end

function  gaussianMixtureValue(hObject, eventdata, handles)
%get value of slider
value = get(hObject, 'value');
%convert integer to floating number
mb = value/100;
ratio = 3;
b1 = mb;
b2 = ratio*mb;
syms a;
eqn =ratio*a*exp(-b1/2) - a*exp(-b2/2) == sqrt(2*pi);
ma = double(solve(eqn, a, 'Real', true));

x = -10:0.01:10;
y = ratio*ma*exp(-b1*x.^2/2)/sqrt(2*pi) - ma*exp(-b2*x.^2/2)/sqrt(2*pi);
z = rectangularPulse(-1, 1, x);
w = 1/exp(-0.5) *  exp(-x.^2/2);
plot(x, y, 'g')
hold on
plot(x, z, 'b')
plot(x, w, 'c')
hold off
mb 
ma

end