function [Sxint,Spxint] = deval(sol,xint,idx)
%DEVAL  Evaluate the solution of a differential equation problem.
%   SXINT = DEVAL(SOL,XINT) evaluates the solution of a differential equation
%   problem at all the entries of the vector XINT. SOL is a structure returned
%   by an initial value problem solver (ODE45, ODE23, ODE113, ODE15S, ODE23S,
%   ODE23T, ODE23TB, ODE15I), a boundary value problem solver (BVP4C, BVP5C),
%   or a delay differential equations solver (DDE23, DDESD). The elements of
%   XINT must be in the interval [SOL.x(1) SOL.x(end)]. For each I, SXINT(:,I)
%   is the solution corresponding to XINT(I).
%
%   SXINT = DEVAL(SOL,XINT,IDX) evaluates as above but returns only
%   the solution components with indices listed in IDX.
%
%   SXINT = DEVAL(XINT,SOL) and SXINT = DEVAL(XINT,SOL,IDX) are also acceptable.
%
%   [SXINT,SPXINT] = DEVAL(...) evaluates as above but returns also the value
%   of the first derivative of the polynomial interpolating the solution.
%
%   For multipoint boundary value problems or initial value problems extended
%   using ODEXTEND, the solution might be discontinuous at interfaces.
%   For an interface point XC, DEVAL returns the average of the limits from
%   the left and right of XC. To get the limit values, set the XINT argument
%   of DEVAL to be slightly smaller or slightly larger than XC.
%
%   Class support for inputs SOL and XINT:
%     float: double, single
%
%   See also
%       ODE solvers:  ODE45, ODE23, ODE113, ODE15S,
%                     ODE23S, ODE23T, ODE23TB, ODE15I,
%                     RKN1210.
%       DDE solvers:  DDE23, DDESD
%       BVP solvers:  BVP4C, BVP5C


% Please report bugs and inquiries to:
%
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@gmail.com    (personal)
%              oldenhuis@luxspace.lu  (professional)
% Affiliation: LuxSpace s�rl
% Licence    : BSD


% If you find this work useful, please consider a donation:
% https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=6G3S5UYM7HJ3N
    
    if ~isa(sol, 'struct')
        % Try DEVAL(XINT,SOL)
        temp = sol;
        sol  = xint;
        xint = temp;
    end
    
    if isstruct(sol)
        if isfield(sol, 'solver')
            solver = sol.solver;
        else
            msg = ['Missing ''solver'' field in the structure ''' inputname(1) '''.'];
            if isfield(sol,'yp')
                warning('MATLAB:deval:MissingSolverField',['%s\n         DEVAL will ' ...
                    'treat %s as an output of the MATLAB R12 version of BVP4C.'], ...
                    msg, inputname(1));
                solver = 'bvp4c';
            else
                error('MATLAB:deval:NoSolverInStruct', msg);
            end
        end
    else
        error('MATLAB:deval:NoSolverStruct', ...
            ['The first or second input to DEVAL must be a structure returned by a ',...
            'differential equation solver.']);        
    end
    
    % Delegate if solver is not RKN1210
    prevPath = pwd;
    cd([matlabroot filesep 'toolbox' filesep 'matlab' filesep 'funfun']);
    deval_original = @deval;
    cd(prevPath);      
    switch solver
        
        %case 'rkn1210'
        %    % TODO: make magic happen
            
        otherwise
            % call builtin DEVAL()
            [Sxint, Spxint] = deval_original(sol,xint,idx);
            return;
    end
    
    % TODO: make MORE magic happen here
    
end