function r1gaussianMixtureFittingAB
    S.x = -5:.1:5;
    S.fh = figure('units','normalized',...
    'Position', [0.515 0.025 0.415 0.87],... %%%%
              'name','slider_plot');

     S.ma = 0.5;
     S.mb = 3;
     S.k = 10;
     
     S.maSlider = uicontrol('style','slider',...
                                       'unit','normalized',...
                                       'position',[0.2 0.0 0.7 0.01],...
                                       'min',0,'max',5,...
                                       'sliderstep',[0.1, 0.1],...
                                       'callback', {@gaussianMixtureValue, 'ma'}); 
                                   
     S.mbSlider = uicontrol('style','slide',...
               'unit','normalized',...
               'position',[0.2 0.03 0.7 0.01],...
               'min',0,'max',3,...
               'sliderstep',[0.1, 0.1],...
               'callback', {@gaussianMixtureValue, 'mb'});
           
     S.kSlider = uicontrol('style','slide',...
               'unit','normalized',...
               'position',[0.2 0.06 0.7 0.01],...
               'min',0,'max',10,...
               'sliderstep',[1, 1],...
               'callback', {@gaussianMixtureValue, 'k'});
     
    
     guidata(S.fh, S);  

end

function gaussianMixtureValue(hSlider, eventData, param)
S = guidata(hSlider); 
if(param == 'ma')
%     disp("ma value is ")
%     get(hSlider, 'value')
    S.ma = get(hSlider, 'value');
elseif(param == 'mb')
%     disp("mb value is ")
%     get(hSlider, 'value')
    S.mb = get(hSlider, 'value'); 
elseif(param == 'k')
    S.k = get(hSlider, 'value'); 
end
x = -5:0.01:5;
y = S.k*S.ma*exp(-S.mb*x.^2/2)/sqrt(2*pi) - S.ma*exp(-S.k*S.mb*x.^2/2)/sqrt(2*pi);
z = rectangularPulse(-1, 1, x);
w = 1/exp(-0.5) *  exp(-x.^2/2);
axis equal
plot(x, y, 'g', 'Linewidth', 2)
hold on
plot(x, z, 'b', 'Linewidth', 2)
plot(x, w, 'c', 'Linewidth', 2)
hold off
disp("mb value is ")
S.mb
disp("ma value is ")
S.ma
disp("k value is ")
S.k
guidata(S.fh, S);
end
