function [fhandle,exitflag,message] = curvefit(points,varargin)

    options = setOptions(varargin)

end


function [options]=setOptions(input)
    % Here you can add new options if needed
    p = inputParser;
    p.CaseSensitive = false;
    % Helper functions for input parser
    checkEmpetyOrChar = @(x) (isempty(x) || ischar(x));
    checkEmptyOrNumericPositive = @(x) (isempty(x) || (isnumeric(x) && all(x > 0)));
    checkNumericPositive = @(x) ((isnumeric(x) && all(x > 0)));
    % General settings
    p.addParameter('display','off',  @(x) checkEmpetyOrChar(x));
    p.addParameter('plot',false,  @(x) islogical(x));

    % Curve settings
    p.addRequired('nCurves', @(x) checkNumericPositive(x));
    p.addRequired('startpos', @(x) checkNumericPositive(x));
    p.addParameter('endpos',[], @(x) checkEmptyOrNumericPositive(x));
    p.addRequired('curveOrder',  @(x) checkNumericPositive(x));
    p.addParameter('curveContinuity',2,  @(x) checkNumericPositive(x));
    p.addParameter('floating',false,  @(x) islogical(x));

    % Additional constraints
    p.addParameter('c0',[],  @(x) checkEmptyOrNumericPositive(x));
    p.addParameter('c1',[],  @(x) checkEmptyOrNumericPositive(x));
    p.addParameter('pointlb',[],  @(x) checkEmptyOrNumericPositive(x));
    p.addParameter('pointub',[],  @(x) checkEmptyOrNumericPositive(x));
    p.addParameter('kappa',[],  @(x) checkEmptyOrNumericPositive(x));

    % Solver settings
    p.addParameter('method','bound',  @(x) checkEmpetyOrChar(x));
    p.addParameter('solver','fminslp',  @(x) checkEmpetyOrChar(x));

    % pars input
    if nargin < 1 || isempty(input)
        parse(p);
    else
        parse(p,input{:});
    end

    % Output results to options structure
    options = p.Results;
end