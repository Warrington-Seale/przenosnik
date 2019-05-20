clc;
clear all;

t0 = 0;
tk = 100;
ts = [t0, tk];

y0 = [0;0;0;0;0;0;0;1;0;0];


[t, y] = ode45('modelv1',ts,y0);