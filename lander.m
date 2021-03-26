close all; clear classes; clear all; clc;
format longG;
format compact;

%
% BODY CONSTANTS
%

body.mu = 4.90486959e+12;          % m^3/s^2; moon
body.r  = 1.7371e+6;               % m; moon radius
body.w  = 2.6616665e-6 * [0,0,1]'; % rad/s; sidereal angular

%
% VEHICLE
%

stages(1).m0     = 10334 + 4700;  % kg
stages(1).thrust = 45040;         % N
stages(1).isp    = 311;           % sec
stages(1).bt     = 555;           % sec

% we shouldn't use the ascent stage to land but this demonstrates multistage
stages(2).m0     = 4700;      % kg
stages(2).thrust = 16 * 1000; % N
stages(2).isp    = 311;       % sec
stages(2).bt     = 463.250;   % sec

r0 = [ 1.7371e+6*1.1 0 0]';
v0 = [-100 100 0]' * 5;

ct = land(stages, body, r0, v0)

function ct = land(stages, body, r0, v0)
  phases = stages;

  % insert our coast phase
  phases = [ phases(1) phases ];
  phases(1).thrust = 0;
  phases(1).bt = 120;

  ct = solver(phases, body, r0, v0);
end

function ct = solver(phases, body, r0, v0)
  global indexes r_scale v_scale r0_bar v0_bar g_bar Nphases rT_bar

  %
  % INTERNAL CONSTANTS
  %

  g0     = 9.80665;                 % m/s; standard gravity
  Nphases = length(phases);
  indexes.r = 1:3;
  indexes.v = 4:6;
  indexes.m = 7;

  %
  % PROBLEM SCALING / SANITY
  %

  g_bar = body.mu / norm(r0)^2;
  r_scale = norm(r0);
  v_scale = sqrt( norm(r0) * g_bar );
  t_scale = sqrt( norm(r0) / g_bar );

  r0_bar = r0 / r_scale;
  v0_bar = v0 / v_scale;

  body.r_bar = body.r / r_scale;
  body.w_bar = body.w * t_scale;

  for p = 1:Nphases;
    phases(p).ve     = phases(p).isp * g0;
    phases(p).a0     = phases(p).thrust / phases(p).m0;
    phases(p).tau    = phases(p).ve / phases(p).a0;
    phases(p).mdot   = phases(p).thrust / phases(p).ve; % informational: not used
    phases(p).dv     = -phases(p).ve * log(1 - phases(p).bt / phases(p).tau); % informational: not used
    phases(p).c      = g0 * phases(p).isp / t_scale;
    phases(p).bt_bar = phases(p).bt / t_scale;
    phases(p)
  end

  %
  % INITIAL STATE AND COAST GUESS
  %

  x0(indexes.r)  = r0_bar;
  x0(indexes.v)  = v0_bar;
  x0(indexes.m)  = phases(1).m0;

  %
  % MAIN SOLVER
  %

  ct_bar = fzero(@(ct) residualFunction(ct, x0, phases, body), 0);

  phases(1).bt_bar = ct_bar;

  [xf] = multipleShooting(x0, phases, body);

  rf = xf(indexes.r)'
  vf = xf(indexes.v)'

  dot(rf, vf)

  ct = ct_bar * t_scale;
end

function z = residualFunction(ct, x0, phases, body)
  global indexes

  phases(1).bt_bar = ct;

  [xf] = multipleShooting(x0, phases, body);

  rf = xf(indexes.r)';
  vf = xf(indexes.v)';

  % FIXME: this should look up the actual height of the terrain at lat,lng of rf
  z = norm(rf) - body.r_bar;
end


function [z,isterminal,direction] = zeroHorizontalEvent(t, X)
  global indexes

  r  = X(indexes.r)';
  v  = X(indexes.v)';

  z = dot(r,v);
  isterminal = 1;  % halt integration
  direction = 1;   % only when increasing from negative to positive
end

function [xf] = multipleShooting(x0, phases, body)
  global indexes Nphases
  ode45options = odeset('RelTol',1e-8,'AbsTol',1e-10,'Events',@zeroHorizontalEvent);

  x = x0;

  for p = 1:Nphases
    bt = phases(p).bt_bar;
    if bt ~= 0
      [ts,xs,te,xe,ie] =  ode45(@(t,x) EOM(t, x, p, phases, body), [0 bt], x, ode45options);
      x = xs(end,:);
      if ~isempty(te)
        break % do not integrate any further upper stages (this may be weird if we coast through periapsis)
      end
    end
  end

  xf = x;
end

%
% EQUATIONS OF MOTION
%

function dX_dt = EOM(t, X, p, phases, body)
  global indexes g_bar

  thrust   = phases(p).thrust;
  c        = phases(p).c;

  r  = X(indexes.r);
  v  = X(indexes.v);
  m  = X(indexes.m);

  vhT = cross(body.w_bar, r);
  vv = dot(v, r) / norm(r);
  vh = v - vv;

  % the 0.0001 here is a hack to avoid terminal guidance issues and keep the rocket pointing up
  u = (0.0001-vv) * r/norm(r) + ( vhT - vh );
  u = u/norm(u);
  %u = r/norm(r);
  T = thrust / (m * g_bar);

  r2 = dot(r,r);
  r3 = r2^(3/2);
  r5 = r2 * r3;

  rdot  = v;
  vdot  = - r / r3 + T * u;

  if thrust == 0
    mdot = 0;
  else
    mdot = - thrust / c;
  end

  dX_dt = [ rdot' vdot' mdot ]';
end

