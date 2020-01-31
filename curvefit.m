function [fhandle,exitflag,message] = curvefit(points,varargin)
    
    options = setOptions(varargin)

end

function [f,df] = HS(xs,ps,beta)
  % Unit step function / projection filter
  % xs is the normalized design variable associated with either a start or end position of a curve.
  % ps is the normalized position of a point which we try to fit the curve to. 
  % beta is a slope parameter for this particular unit step approximation. 
  % The formulation is based on Wang,  F.,  Lazarov,  B.,  and  Sigmund,  O 2011 and Soerensen, R, and Lund, E 2015
    f = (tanh(beta*xs)+tanh(beta*(ps-xs)))./(tanh(beta*xs)+tanh(beta*(1-xs)));
    if nargout > 1
        df =((1-tanh(beta*xs).^2)-(1-tanh(beta*(ps-xs)).^2)*beta)/(tanh(beta*xs)+tanh(beta*(1-xs))) ...
        - (tanh(beta*ps)+tanh(beta*(xs-ps)))*beta*((1-tanh(beta*xs).^2)-(1-tanh(beta*(1-xs)).^2))./(tanh(beta*ps)+tanh(beta*(1-xs))).^2;
    end

end

function [options]=setOptions(input)
    % Here you can add new options if needed
    p = inputParser;
    p.CaseSensitive = false;
    % Helper functions for input parser
    checkEmpetyOrChar = @(x) (isempty(x) || ischar(x));
    checkEmptyOrNumericPositive = @(x) (isempty(x) || (isnumeric(x) && all(x > 0)));
    checkNumericPositive = @(x) ((isnumeric(x) && all(x > 0)));

    % Curve settings
    p.addRequired('startpos', @(x) isnumeric(x));
    p.addRequired('nCurves', @(x) checkNumericPositive(x));
   
    p.addRequired('curveOrder',  @(x) checkNumericPositive(x));
    p.addParameter('endpos',[], @(x) isnumeric(x));
   
    p.addParameter('curveContinuity',2,  @(x) checkNumericPositive(x));
    p.addParameter('floating',false,  @(x) islogical(x));
    p.addParameter('beta',100,  @(x) checkNumericPositive(x));

    % General settings
    p.addParameter('display','off',  @(x) checkEmpetyOrChar(x));
    p.addParameter('plot',false,  @(x) islogical(x));
    
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