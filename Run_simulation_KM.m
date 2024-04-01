
close all
clear all
load('SimulationParameters.mat')


%%
%%%%%% #1 Generate the cell structures and build a cell mask   %%%%%%
% Each simulation step takes ~ 30 min and 20 MB of storage
DiffSim_ToolBox.Simulate_cell_structure(Sim.Cell)

%%
%%%%%% #2 Simulate the diffusion water distribution from the cell mask   %%%%%%
% Each simulation step takes ~ 30 s and 500 MB of storage
DiffSim_ToolBox.Simulate_diffusion_distribution(Sim.Diff)

%%
%%%%%% #3 Simulate the MR phase distribution from the water distribution   %%%%%%
% Each simulation step takes ~ 20 s and 500 MB of storage
DiffSim_ToolBox.Simulate_MR_distribution(Sim.MR)