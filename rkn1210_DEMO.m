% RKN1210_DEMO         small demonstration of RKN1210-integrator
%
% Run this function for a small demonstration of the Runge-Kutta-Nystrom
% 12th/10th order integrator for second-order ODE's. 
%
% See also RKN1210, RKN86, ODE86, ODE113, DEVAL, ODESET, ODEGET. 
function rkn1210_DEMO

% Please report bugs and inquiries to: 
%
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@gmail.com    (personal)
%              oldenhuis@luxspace.lu  (professional)
% Affiliation: LuxSpace sàrl
% Licence    : BSD
%
% Please report any bugs or suggestions to oldenhuis@gmail.com


% If you find this work useful and want to show your appreciation:
% https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=6G3S5UYM7HJ3N

    %% initialize
    
    % optimize figure window sizes
    scz      = get(0, 'screensize');
    position = @(width,height) [scz(3:4)/2-[width, height]/2, width, height];
  
    oldInterpreter = get(0, 'defaulttextinterpreter');
    set(0, 'defaulttextinterpreter', 'none')
    
    %% Elementary example: Circular orbit
    
    % define function
    f = @(t,y) [-y(1) * (y(1)^2+y(2)^2)^(-3/2);
                -y(2) * (y(1)^2+y(2)^2)^(-3/2)];

    % show message
    uiwait(msgbox({
        'The class of Runge-Kutta-Nystrom (RKN) integrators are very efficient for ODE''s of the type'
        ''
        '        y'''' = F(t,y)'
        ''
        'i.e., where the first derivative of the function [y] does not explicitly appear in the ODE.'
        ''
        ['One example of a problem of this type is a body in freefall, in an ',...
        'environment without air drag or any other resistence; a body in orbit, for example.'];
        ''
        ['As a first example, we will integrate the trajectory of a massless body in a ',...
        'simple circular orbit around a massive body, for 25 revolutions:']}, ...
        'First demonstration','modal'));

    % start integration
    options = odeset('abstol', 5e-16, 'reltol', 5e-16);
    [t, y, dummy,dummy, output] = rkn1210(f, [0, 50*pi], [1; 0], [0; 1],options);%#ok

    % show results
    
    % The orbit should be *exactly* circular
    figure(1), clf, hold on
    set(1, 'position', position(800, 600));
    subplot(2,2,1)
    plot(y(:,1), y(:,2), '.k', 'markersize', 1);
    title('Integrated circular orbit'), xlabel('X'), ylabel('Y')
    axis equal, axis tight
    
    % We can easily plot the x- and y-error
    subplot(2,2,2), hold on
    Xerr = abs(y(:,1) - cos(t));
    Yerr = abs(y(:,2) - sin(t));
    plot(sqrt(Xerr.^2+Yerr.^2), 'r');
    title('Absolute error')
    xlabel('step'), ylabel('error')
    axis tight
    
    % Also plot the step magnitude and error estimate
    subplot(2,2,3), hold on
    plot(output.h, 'b')    
    title('Step size [s]') 
    xlabel('step'), ylabel('magnitude')
    axis tight
    
    subplot(2,2,4), hold on
    plot(output.delta, 'r')
    title('Error estimate') 
    xlabel('step'), ylabel('magnitude')
    axis tight

    %% Elliptic orbit
    
    % show message
    pause(1)
    uiwait(msgbox({[
        'Note how the error remains below 1.5e-13, and how the step size ',...
        'starts conservatively but automatically adapts itself to the largest possible value that ',...
        'still gives the desired accuracy.']
        ''
        'That''s pretty good, but of course it''s an exceedingly simple example.';
        ' '
        [' Now let''s try a more realistic example; a 200×10000 km elliptical ',...
        'orbit around the Earth, including the Earth''s oblateness around the equator, for 5 ',...
        'revolutions:']}, ...
        'Second demonstration','modal')); 
    
    % define differential equation
    muE = 398600.4418;
    RE  = 6738;
    RE2 = RE*RE;
    J2  = 0.00108263;
    function d2ydt2 = f2(t,y)%#ok
        % some renaming
        rmag = norm(y);
        % main term
        d2ydt2 = -muE*y/rmag^3;
        % J2 perturbation
        pp = 3/2*muE*J2*RE2/rmag^5;
        z2r2 = 5*(y(3)/rmag)^2;
        d2ydt2 = d2ydt2 - pp*y.*[1-z2r2; 1-z2r2; 3-z2r2];
    end      
    
    % start integration
    options = odeset('abstol',5e-14, 'reltol',5e-14); % loosen constraints a bit (for speed)
    ra  = RE+10000;  % apogee
    rp  = RE+200;    % perigee
    a   = (rp+ra)/2; % semi-major axis
    T   = 2*pi*sqrt(a^3/muE);   % orbital period
    Vp  = sqrt(muE*(2/rp-1/a)); % speed at perigee
    y0  = [rp 0 0]; % initial position
    yp0 = [0 Vp Vp]*sqrt(2)/2; % initial velocity
    [t, y, dummy,dummy, output] = rkn1210(@f2, [0 5*T], y0, yp0, options);%#ok
        
    % show results
    figure(1), clf, hold on
    set(1, 'position', position(800, 600));
    subplot(1,2,1), hold on
    [xyz{1:3}] = sphere(20);
    xyz = cellfun(@(x)x*RE,xyz,'uniformoutput',false);
    mesh(xyz{:}, 'edgecolor', [.5 .5 .5]); hold on
    plot3(y(:,1), y(:,2), y(:,3), 'r', 'markersize', 1);
    title('Integrated elliptic orbit')
    xlabel('X [km]'), ylabel('Y [km]'), ylabel('Z [km]')
    axis equal, axis tight, view(-102,44)
    
    % Also plot the step magnitude and error estimate
    subplot(2,2,2), hold on
    plot(output.h, 'b')    
    title('Time step per iteration') 
    xlabel('step #'), ylabel('Step size [s]')
    axis tight
    
    subplot(2,2,4), hold on
    plot(output.delta, 'r')
    title('Error estimate') 
    xlabel('step #'), ylabel('magnitude [km]')
    axis tight
    
    %% Using output functions
    
    % show message
    pause(1)
    uiwait(msgbox({[
        'Note how the accumulated roundoff error remains very small for this more ',...
        'complicated case. Also note how the integrator automatically decreases its stepsize each ',...
        'time the satellite approaches the Earth more closely. This is because its speed increases ',...
        'at those points, so that the error will accumulate more rapidly if large stepsizes are ',...
        'maintained. These are also the regions where the ODE45 integrator often fails to produce ',...
        'the correct trajectory.']
        ' '
        ['This RKN1210-function also supports ''output functions'', functions that are called ',...
        'after each successful integration step. This can be used to produce a progress bar, for example:']},...
        'Third demonstration','modal'));
                
    % initialize waitbar
    wait = waitbar(0, 'integrating circular orbit...');

    % start integration (circular orbit again)
    options = odeset('abstol', 1e-16, 'outputfcn', @OutputFcn1);
    [t, y, dummy,dummy, output] = rkn1210(f, [0, 50*pi], [1; 0], [0; 1],options);%#ok
    % kill waitbar
    close(wait)

    % the output function
    function stop = OutputFcn1(t,y,dy,flag)%#ok
        % don't stop
        stop = false;
        % only after sucessfull steps
        if ~isempty(flag), return
        else
            wait = waitbar(t/50/pi, wait);
        end
    end
    
    % Plot the circular orbit
    figure(1), clf, hold on
    set(1, 'position', position(800, 600));
    plot(y(:,1), y(:,2), '.k', 'markersize', 1);
    title('Integrated circular orbit'), xlabel('X'), ylabel('Y')
    axis equal, axis tight
    
    %% Using event functions
  
    pause(1)
    uiwait(msgbox({[
        'Aside from output functions, this integrator also supports ''event functions''. ',...
        'Such functions can be used to detect where user-defined events occur, possibly halting the ',...
        'integration after these events have occured.']
        ' '
        ['As a final example, we will use an output function ',...
        'to display the instantaneous integration step, and an event function that detects when the ',...
        'Y-coordinate becomes negative. When this happens, the integration is terminated:']}, ...
        'Fourth demonstration','modal'));
    
    % initialize figure
    figure(1), clf, hold on
                
    % start integration (circular orbit again)    
    options = odeset('abstol', 1e-16, 'outputfcn', @OutputFcn2, 'events', @EventFcn);
    [t, y, dy, TE,YE] = rkn1210(f, [0, 50*pi], [1; 0], [0; 1],options);%#ok
        
    % plot final point
    plot(YE(1), YE(2), 'rx')
    text(-.9, 0 ,'<--  Event: y < 0; integration stopped', 'interpreter', 'none')
            
    % the output function
    function stop = OutputFcn2(t,y,dy,flag)%#ok
        % don't stop
        stop = false;
        % only after sucessfull steps
        if ~isempty(flag), return
        else
            plot(y(1), y(2), 'b.')
            axis equal, axis([-1.2 1.2 -.2 1.2])            
            pause(0.1), drawnow % flush plotting commands
        end
    end
    
    % the event function
    function [value, isterminal, direction] = EventFcn(t,y,dy)%#ok
        % stop upon event? Yup:
        isterminal = true;
        % direction should be DEcreasing (pos to neg)
        direction = -1;
        % value is simply Y-coordinate
        value = y(2);
    end
      
    %% thanks & enjoy
    
    pause(1)
    uiwait(msgbox(['That''s about it. Please look inside this Demo to learn how to use these features.',...
    'Please send any bugs you find to oldenhuis@gmail.com'],...
    'That''s all Folks!','modal'));
    close(1)
    
    set(0, 'defaulttextinterpreter', oldInterpreter);
    
end
 
