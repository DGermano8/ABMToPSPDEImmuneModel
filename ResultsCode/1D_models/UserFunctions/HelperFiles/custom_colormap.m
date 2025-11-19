function cmap = custom_colormap(color1, color2, color3, numPoints)
    % CUSTOM_COLORMAP Generates a colormap transitioning through three colors.
    %
    % Inputs:
    %   color1    - 1x3 vector specifying the RGB values of the first color
    %   color2    - 1x3 vector specifying the RGB values of the second color
    %   color3    - 1x3 vector specifying the RGB values of the third color
    %   numPoints - Total number of colors in the colormap
    %
    % Output:
    %   cmap - numPoints x 3 colormap array
    numPoints = numPoints + 1;
    % Ensure colors are row vectors
    color1 = color1(:)';
    color2 = color2(:)';
    color3 = color3(:)';
    
    % Split numPoints roughly equally between two transitions
    n1 = round(numPoints / 2);
    n2 = numPoints - n1;
    
    % Generate linear transitions between color1->color2 and color2->color3
    cmap1 = [linspace(color1(1), color2(1), n1)', linspace(color1(2), color2(2), n1)', linspace(color1(3), color2(3), n1)'];
    cmap2 = [linspace(color2(1), color3(1), n2)', linspace(color2(2), color3(2), n2)', linspace(color2(3), color3(3), n2)'];
    
    % Combine both segments
    cmap = [cmap1; cmap2(2:end, :)]; % Avoid repeating the middle color
    
    % Ensure values stay within valid RGB range [0,1]
    cmap = min(max(cmap, 0), 1);
end