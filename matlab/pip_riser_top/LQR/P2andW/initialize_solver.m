FUN = @v4_new_4d_model;
if (odepar.pair==0),
   k_ = x(:,1)*zeros(1,7);  % k_ needs to be initialized as an Nx7 matrix where N=number of rows in x
                       % (just for speed so octave doesn't need to allocate more memory at each stage value)
   % Compute the first stage prior to the main loop.  This is part of the Dormand-Prince pair caveat.
   % Normally, during the main loop the last stage for x(k) is the first stage for computing x(k+1).
   % So, the very first integration step requires 7 function evaluations, then each subsequent step
   % 6 function evaluations because the first stage is simply assigned from the last step's last stage.
   % note: you can see this by the last element in c_ is 1.0, thus t+c_(7)*h = t+h, ergo, the next step.
   k_(:,1)=feval(FUN,t,x(:,1),u(:,1),'derivatives',par); % first stage
    if odepar.count==1,
    global rhs_counter
    if ~exist('rhs_counter'),rhs_counter=0; end
    end % if count
elseif (odepar.pair==1)
    k_ = x(:,1)*zeros(1,6);  % k_ needs to be initialized as an Nx6 matrix where N=number of rows in x
                       % (just for speed so octave doesn't need to allocate
                       % more memory at each stage value)
end

hi = odepar.hmin; % initial guess at a step size