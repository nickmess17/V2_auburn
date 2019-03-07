%script to run GA on MyLake
clc,clear,close all
fun=@MyLake_Auburn_run;
nvars=5;%number of parameters to calibrate
lb=[0.5 0.1 10 0.000008 .000008];%lower bounds of parameters
ub=[3 2 25 0.008 0.008];%upper bounds of parameters
population_size = 50;  % Populations size for each generation of the genetic algorithm
max_generations = 20;  % How many generations to run the genetic algorithm for
parallelize     = false; % runs function in series
options=optimoptions('ga', 'MaxGenerations', max_generations, 'PopulationSize', population_size, 'UseParallel', parallelize);
[par,fval]=ga(fun,nvars,[],[],[],[],lb,ub,[],[],options);