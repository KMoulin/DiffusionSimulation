
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This code generates structures that can be used to generated cells masks
%   and then simulate a diffusion distribution and a MR simulation. 
%
%   Each Nb cell structure is generated once   
%   Each Nb Diffusion distribution is generated per cell structure
%   Each Nb MR distribution is generated per diffusion distribution
%
%   The duration and storage of the simulation is roughtly :
% 
%   Nb_cell * 200Mb                         Nb_cell * 30 min
%   Nb_diff * Nb_cell * 500Mb               Nb_diff * Nb_cell * 30s
%   Nb_diff * Nb_cell * Nb_MR * 500Mb       Nb_diff * Nb_cell * Nb_MR *20s
%
%   45 cells structures * 5 Diff distribution * 1 MR distribution =
%   234Gb/26hours 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all
clear all

Sim=[];
Sim.dir= uigetdir;
cd(Sim.dir)
Sim.Cell.dir=[Sim.dir '\CellStructure'];
Sim.Diff.dir=[Sim.dir '\Diffusion'];
Sim.MR.dir  =[Sim.dir '\MRseq'];

mkdir(Sim.Cell.dir);
mkdir(Sim.Diff.dir);
mkdir(Sim.MR.dir);


%%%%%% #1 Setup the structures for the Cell mask simulation  %%%%%%

ECV=[0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1]; 

%ECV=[0.50 0.60 0.70 0.80 0.90 1]; 

cpt_ECV=1;
for cpt_cell=1:1:2
    
    
Sim.Cell.struct{cpt_cell}.ID=cpt_cell; 
%%%%% Voxel parameters %%%%%%
Sim.Cell.struct{cpt_cell}.Voxel=[0.5 0.5 0.5].*1e-3;                                    % mm -> m
Sim.Cell.struct{cpt_cell}.Voxel_vol=Sim.Cell.struct{cpt_cell}.Voxel(1)*Sim.Cell.struct{cpt_cell}.Voxel(2)*Sim.Cell.struct{cpt_cell}.Voxel(3);              % Volume in m3
Sim.Cell.struct{cpt_cell}.Voxel_sur=Sim.Cell.struct{cpt_cell}.Voxel(1)*Sim.Cell.struct{cpt_cell}.Voxel(2);                             % Surface in m2
                                          
Sim.Cell.struct{cpt_cell}.Voxel_pixel=[500 500 500];
Sim.Cell.struct{cpt_cell}.MaskRes=Sim.Cell.struct{cpt_cell}.Voxel./Sim.Cell.struct{cpt_cell}.Voxel_pixel;  % m 

%%%%% Cells parameters %%%%%%

Sim.Cell.struct{cpt_cell}.Cell=[9 20 100 100].*1e-6;                                   % 100um long and 10-25um of diameters
Sim.Cell.struct{cpt_cell}.Cell_vol=pi*Sim.Cell.struct{cpt_cell}.Cell(2).*(Sim.Cell.struct{cpt_cell}.Cell(1)/2)^2;                      % Volume in m3 
Sim.Cell.struct{cpt_cell}.Cell_sur=pi*(Sim.Cell.struct{cpt_cell}.Cell(1)/2)^2;                    % Surface in m2 
Sim.Cell.struct{cpt_cell}.ECV=0.60;%ECV(cpt_ECV);                  % Percentage of Extra-cellular volume (%)
Sim.Cell.struct{cpt_cell}.PV=0.15;                                                          % Percentage of perfusion (%)
Sim.Cell.struct{cpt_cell}.ICV=1-Sim.Cell.struct{cpt_cell}.ECV;                              % Percentage of Intra-cellular volume (%)

Sim.Cell.struct{cpt_cell}.Bridge=cpt_cell-1;

Sim.Cell.struct{cpt_cell}.folder_out=[Sim.Cell.dir '\Cell_struct' ]; 

    cpt_ECV=cpt_ECV+1;
    if cpt_ECV>length(ECV)
        cpt_ECV=1;
    end

end


%%%%%% #3 Setup the structures for the MR simulation  %%%%%%

MaxTime=0;
for cpt_MR=1:1
%%%%% Waveforms parameters %%%%%%
Sim.MR.struct{cpt_MR}.ID=cpt_MR;
Sim.MR.struct{cpt_MR}.Gamma= 267.513*10^6;                                                   % rad/(s . Tesla )
Sim.MR.struct{cpt_MR}.TE= 54e-3;                                                             % ms -> s   
Sim.MR.struct{cpt_MR}.Bval= [300];                                                           % s/mm² 
Sim.MR.struct{cpt_MR}.Waveform=1;                                                            % How many waveform we simulate
Sim.MR.struct{cpt_MR}.Refoct=4300;                                                           % Duration of the refocussing pulse in us
Sim.MR.struct{cpt_MR}.MixingTime=1e6;                                                       % Duration of the Mixing Time in us
Sim.MR.struct{cpt_MR}.RampTime=1500;                                                         % Ramp time in us
Sim.MR.struct{cpt_MR}.dT= 10e-6;  
%%% Sequence Timings and Amplitudes
Sim.MR.struct{cpt_MR}.Dn=[
                                7180  12620   (27880-7180-12620)    Sim.MR.struct{cpt_MR}.Refoct      12620  7180   0;
                                8690  8690    (25470-2*8690)        Sim.MR.struct{cpt_MR}.Refoct     8690   8690   0;
                                8300     0    (16380-8300)          Sim.MR.struct{cpt_MR}.Refoct      8300      0   0;

                                8890  9790    1790                 Sim.MR.struct{cpt_MR}.Refoct      8190   8400   0;
                                8990  10500   0                    Sim.MR.struct{cpt_MR}.Refoct     7090   8500   0;
                                9300  2700    0                    Sim.MR.struct{cpt_MR}.Refoct     8100      0   0;];
   
Sim.MR.struct{cpt_MR}.Gn=[   
                                   1    -1    0                         0      1     -1   0; %AMC
                                   1    -1    0                         0     -1      1   0; %Bip
                                   1     0    0                         0     -1      0   0; %Mono

                                   1    -1    1                         0     -1      1   0; %Code M2
                                   1    -1    0                         0     -1      1   0; %Code M1
                                   1    -1    0                         0     -1      0   0;]; %Code M0


Sim.MR.struct{cpt_MR}.Dir=[0.831862183551238,-0.978020640898505,-0.746594402068549,0.618503898022940,0.0357105488849347,0.0690933229952567,0.285311347832409,0.541890310423585,0.696013916160111,-0.586381611851074,-0.361852167119390,-0.107233806522707;
         -0.474849677041022,-0.117331253822020,-0.566355417065775,-0.0500144660490100,-0.296253369933249,-0.451075167295507,-0.874813688131557,-0.803757865561315,-0.264300268733346,-0.487477487497046,-0.922081582695813,-0.849061587106308;
         -0.287268170935018,-0.172362813165683,0.349053404091507,-0.784188429834895,0.954441552711641,-0.889807334332884,-0.391533602085168,0.245617655503534,0.667615050863187,-0.646933288596948,-0.137217684961672,0.517296200740132]';

MaxTime=max(max(sum(Sim.MR.struct{cpt_MR}.Dn(:,:),2)*1e-6),MaxTime);                               % Save the longest sequence duration                     

Sim.MR.struct{cpt_MR}.folder_in=  [Sim.Diff.dir '\Diff_sim'  ] ;  
Sim.MR.struct{cpt_MR}.folder_out= [Sim.MR.dir   '\MR_sim' ];  

end


%%%%%% #2 Setup the structures for the diffusion simulation  %%%%%%

for cpt_diff=1:3
    
Sim.Diff.struct{cpt_diff}.ID=cpt_diff;
%%%%% Diffusion parameters %%%%%%
Sim.Diff.struct{cpt_diff}.N=4000;                                                               % Number of particles
Sim.Diff.struct{cpt_diff}.dim=3;                                                                % Diff parameters                                                             
Sim.Diff.struct{cpt_diff}.D=3*(10^-3)*(10^-6);                                                  % Coefficient of diffusion extravascular mm²/s -> m²/s
Sim.Diff.struct{cpt_diff}.Din=Sim.Diff.struct{cpt_diff}.D*2/3;                                               % Coefficient of diffusion intravascular mm²/s -> m²/s
Sim.Diff.struct{cpt_diff}.V=[0 0 0];                                                            % Linear Velocity along m/s 
Sim.Diff.struct{cpt_diff}.VGrad=[0 0 0];                                                        % Velocity Grad along x m/s 
Sim.Diff.struct{cpt_diff}.Perma=0;                                                              % Cell permability (%)
Sim.Diff.struct{cpt_diff}.dT= 10e-6;                                                            % Diffusion dT in s
Sim.Diff.struct{cpt_diff}.dur= MaxTime/Sim.Diff.struct{cpt_diff}.dT;                       % Duration in us / dT = nb points
Sim.Diff.struct{cpt_diff}.Mt= Sim.MR.struct{cpt_MR}.Mt/Sim.Diff.struct{cpt_diff}.dT;                       % Duration in us / dT = nb points
Sim.Diff.struct{cpt_diff}.Iterative= 0;
Sim.Diff.struct{cpt_diff}.Compartement= cpt_diff-1;
if cpt_diff>2
    Sim.Diff.struct{cpt_diff}.Compartement=[];
end
Sim.Diff.struct{cpt_diff}.folder_in=  [Sim.Cell.dir '\Cell_struct'  ] ;  
Sim.Diff.struct{cpt_diff}.folder_out= [Sim.Diff.dir '\Diff_sim'  ] ; 

end




save([Sim.dir '/SimulationParameters.mat'],'Sim');
       