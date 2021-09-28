% The root directory, used in the problem definition script for saving the
% data. 

root = pwd;
root = [root '/../'];

% add working directory from fastfins
root    = [root, '/fastfins'];
addpath(root);
addpath(genpath([root '/optimizer']));
addpath(genpath([root '/samplers']));
addpath(genpath([root '/solvers']));
addpath(genpath([root '/library']));

% GOMOS model
addpath([pwd '/model']);
addpath([pwd '/utility']);

% set default for figures, you may not need this
figure_default;
