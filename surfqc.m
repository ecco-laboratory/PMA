function h = surfqc(varargin)
    %SURFC  Combination surf/contour plot.
    %   SURFC(...) is the same as SURF(...) except that a contour plot
    %   is drawn beneath the surface.
    %
    %   See also SURF, SHADING.
    
    %   Clay M. Thompson 4-10-91
    %   Copyright 1984-2017 The MathWorks, Inc.

    if nargin > 0
        [varargin{:}] = convertStringsToChars(varargin{:});
    end
    
    [~, cax, args] = parseplotapi(varargin{:},'-mfilename',mfilename);
    [reg, prop] = parseparams(args);
    nargs = length(reg);
    
    if nargs < 1
        error(message('MATLAB:narginchk:notEnoughInputs'));
    elseif nargs > 4
        error(message('MATLAB:narginchk:tooManyInputs'));
    end
    if rem(length(prop), 2) ~= 0
        error(message('MATLAB:surfc:InvalidPVPair'))
    end
    
    if nargs == 1  % Generate x, y matrices for surface z.
        z = reg{1};
        [m, n] = size(z);
        [x, y] = meshgrid(1 : n, 1 : m);
    elseif nargs == 2
        z = reg{1};
        [m, n] = size(z);
        [x, y] = meshgrid(1 : n, 1 : m);
    elseif nargs == 3
        [x, y, z] = deal(reg{1 : 3});
    elseif nargs == 4
        [x, y, z] = deal(reg{1 : 3});
    end
    x = matlab.graphics.chart.internal.datachk(x,'numeric');
    y = matlab.graphics.chart.internal.datachk(y,'numeric');
    z = matlab.graphics.chart.internal.datachk(z,'numeric');
    
    if min(size(z)) == 1
        error(message('MATLAB:surfc:NonMatrixInput'));
    end
    
    % Determine state of system
    cax = newplot(cax);
    nextPlot = cax.NextPlot;
    
    % Plot surface
    hs = surf(cax, args{:});
    
    % Set NextPlot to 'add' so that the contour object is added to the
    % existing axes. 'surf' calls 'newplot', so the Figure's NextPlot
    % property will already be set to 'add' at this point.
    cax.NextPlot = 'add';
    
    a = get(cax, 'ZLim');
    
    % Always put contour below the plot at the ZMin.

    % Get D contour data
    [~, hh] = contour(cax, x, y, z, 'ZLocation', "ZMin");

    [u,v] = gradient(-1*downsample((downsample(z,5))',5));
    xd = downsample((downsample(x,5))',5);
    yd = downsample((downsample(y,5))',5);
    quiver3(xd,yd,min(hh.ZData(:))*ones(size(xd)),u,v,zeros(size(u)),'k-');

    % Disable data tips on the contour
    d = hggetbehavior(hh, 'DataCursor');
    d.Enable = false;
    
    % Restore the original value for NextPlot.
    cax.NextPlot = nextPlot;
    
    if nargout > 0
        h = [hs; hh(:)];
    end
end
