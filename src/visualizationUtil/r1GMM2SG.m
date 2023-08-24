function r1GMM2SG
    S.x = -5:.1:5;
    S.fh = figure('units','normalized',...
    'Position', [0.515 0.025 0.415 0.87],... %%%%
              'name','slider_plot');
          
     a=[2.9436   15.4697];
     b=[30.0000    3.6882];
     
     % Declar the member variables of class S
     S.a1=a(1); S.a2=a(2);
     S.b1=b(1); S.b2=b(2);
     
     distance = 0.03;
     
    S.a1Slider = uicontrol('style','slider',...
                   'unit','normalized',...
                   'position',[0.2 3*distance 0.7 0.01],...
                   'min',-10,'max',30,...
                   'sliderstep',[0.1, 0.1],...
                   'callback', {@gaussianMixtureValue, 'a1'}); 

     S.a2Slider = uicontrol('style','slide',...
               'unit','normalized',...
               'position',[0.2 2*distance 0.7 0.01],...
               'min',-10,'max',30,...
               'sliderstep',[0.1, 0.1],...
               'callback', {@gaussianMixtureValue, 'a2'});
           
     S.b1Slider = uicontrol('style','slide',...
               'unit','normalized',...
               'position',[0.2 1*distance 0.7 0.01],...
               'min',-10,'max',30,...
               'sliderstep',[1, 1],...
               'callback', {@gaussianMixtureValue, 'b1'});
           
     S.b2Slider = uicontrol('style','slide',...
               'unit','normalized',...
               'position',[0.2 0.0 0.7 0.01],...
               'min',-10,'max',30,...
               'sliderstep',[1, 1],...
               'callback', {@gaussianMixtureValue, 'b2'});
           
     guidata(S.fh, S);  

end


function gaussianMixtureValue(hSlider, eventData, param)
S = guidata(hSlider); 
if(param == 'a1')
    S.a1 = get(hSlider, 'value');
elseif(param == 'a2')
    S.a2 = get(hSlider, 'value'); 
elseif(param == 'b1')
    S.b1 = get(hSlider, 'value'); 
elseif(param == 'b2')
    S.b2 = get(hSlider, 'value'); 
end
x = -5:0.01:5;

z = rectangularPulse(-1, 1, x);
w = 1/exp(-0.5) *  exp(-x.^2/2);
y = S.a1 * exp(-S.b1 * x.^2/2)/sqrt(2*pi) + S.a2 * exp(-S.b2 * x.^2/2)/sqrt(2*pi);

plot(x, y, 'g')
hold on
plot(x, z, 'b')
plot(x, w, 'c')
hold off

disp("a is ")
[S.a1 S.a2]
disp("b is ")
[S.b1 S.b2]

guidata(S.fh, S);
end
