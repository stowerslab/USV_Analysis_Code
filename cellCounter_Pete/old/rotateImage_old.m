function [ Irotated ] = rotateImage( Ioriginal )
%rotateImage - interactive GUI image rotation function

    %# setup GUI
    hFig = figure('menu','none');
    hAx = axes('Parent', hFig);
    hRotControl = uicontrol('Parent',hFig, 'Style','slider', 'Value',0, 'Min',0,'Max',360, 'SliderStep',[1 10]./360, 'Position',[150 5 300 20], 'Callback',@slider_callback); 
    hDoneControl = uicontrol(hFig, 'Style', 'pushbutton', 'String', 'Done', 'Callback', @doneCallback, 'Position', [10 10 50 20]);    %make button to press when done manipulating
                          
    hTxt = uicontrol('Style','text', 'Position',[290 28 20 15], 'String','0');

    %# show image
    imshow(Ioriginal, 'Parent', hAx);
    waitfor(hDoneControl);
    close(hFig);

    %# Callback function
    function slider_callback(hObj, eventdata)
        angle = round(get(hObj,'Value'));        %# get rotation angle in degrees
        imshow(imrotate(Ioriginal,angle), 'Parent',hAx)  %# rotate image 
        set(hTxt, 'String',num2str(angle))       %# update text
        Irotated = imrotate(Ioriginal,angle,'bilinear','loose');
    end

    function doneCallback(hObject, ~)
        delete(hObject);
    end

end
