warning off;
clear all;

World.Gamma= 267.513*10^6;                                                 % rad/(s . Tesla )
World.dT= 1e-6;    

%%%%% Voxel parameters %%%%%%
World.Voxel=[0.5 0.5 0.5].*1e-3;                                           % mm -> m
World.Voxel_vol=World.Voxel(1)*World.Voxel(2)*World.Voxel(3);              % Volume in m3
World.Voxel_sur=World.Voxel(1)*World.Voxel(2);                             % Surface in m2
                                                                           % us -> s
%%%%% Cells parameters %%%%%%
% 100um long and 10-25um of diameters
World.Cell=[9 20 100 100].*1e-6;                                           % um -> m
World.Cell_vol=pi*World.Cell(2).*(World.Cell(1)/2)^2;                      % Volume in m3 
World.Cell_sur=pi*(World.Cell(1)/2)^2;                                     % Surface in m2 
World.ECV=0.25;                                                            % Purcentage of Extra-cellular volume (%)
World.PV=0.15;                                                             % Purcentage of perfusion (%)
World.ICV=1-World.ECV;                                                     % Purcentage of Intra-cellular volume (%)

%%%%% Diffusion parameters %%%%%%
Diff.N=1000;                                                               % Number of particles
Diff.dim=3;                                                                % Diff parameters                                                             
Diff.D=2.2*(10^-3)*(10^-6);                                                % Coefficient of diffusion extravascular mm²/s -> m²/s
Diff.Din=Diff.D*2/3;                                                       % Coefficient of diffusion intravascular mm²/s -> m²/s
Diff.V=[0 0 0];                                                            % Linear Velocity along m/s 
Diff.VGrad=[0 0 0];                                                        % Velocity Grad along x m/s 
Diff.Perma=0;
Diff.dT= 25e-6; 
                                                     
%%%%% Waveforms parameters %%%%%%
Seq.TE= 54e-3;                                                             % ms -> s   
Seq.Bval= [300];                                                           % s/mm² %[0 5 10 15 20 25 30 35 45 50 75 100 200 300]
Seq.Waveform=6;
Seq.Refoct=4300;
Seq.RampTime=1500;

Seq.Dn=[
        7180  12620   (27880-7180-12620)    Seq.Refoct      12620  7180   0;
        8690  8690    (25470-2*8690)        Seq.Refoct      8690   8690   0;
        8300     0    (16380-8300)          Seq.Refoct      8300      0   0;
        
        8890  9790    1790                  Seq.Refoct      8190   8400   0;
        8990  10500   0                     Seq.Refoct      7090   8500   0;
        9300  2700    0                     Seq.Refoct      8100      0   0;];
        
        %3300     0    0                     1011380-3300    3300      0   0;];
   
Seq.Gn=[   
           1    -1    0                         0      1     -1   0; %AMC
           1    -1    0                         0     -1      1   0; %Bip
           1     0    0                         0     -1      0   0; %Mono
           
           1    -1    1                         0     -1      1   0; %Code M2
           1    -1    0                         0     -1      1   0; %Code M1
           1    -1    0                         0     -1      0   0;]; %Code M0
           
           %1     0    0                         0     -1      0   0;]; %STEAM
       
Seq.Dir=[0.831862183551238,-0.978020640898505,-0.746594402068549,0.618503898022940,0.0357105488849347,0.0690933229952567,0.285311347832409,0.541890310423585,0.696013916160111,-0.586381611851074,-0.361852167119390,-0.107233806522707;
         -0.474849677041022,-0.117331253822020,-0.566355417065775,-0.0500144660490100,-0.296253369933249,-0.451075167295507,-0.874813688131557,-0.803757865561315,-0.264300268733346,-0.487477487497046,-0.922081582695813,-0.849061587106308;
         -0.287268170935018,-0.172362813165683,0.349053404091507,-0.784188429834895,0.954441552711641,-0.889807334332884,-0.391533602085168,0.245617655503534,0.667615050863187,-0.646933288596948,-0.137217684961672,0.517296200740132]';
%Seq.Dir=[1,0,0; 0,1,0; 0,0,1];

Diff.dur= ceil((max(sum(Seq.Dn(:,:),2))*World.dT)/(Diff.dT));

%% Simulate a 2D cell structure from the parameters
[Cells] = cell_structure_2D(World.Cell,World.Voxel,World.ECV,1000,5);


%% ADD a Z structure to our cells
Cells=Cellule_ToolBox.Add_Cell_Z(Cells,World.Cell,World.Voxel);

%% Calculate ECV as a checkup
tic
disp(['Generated ECV ' num2str(Cellule_ToolBox.Calculate_ECV(Cells,World.Voxel)) ]);
toc


%% Simulation diffusion inside and outside the cells
[Water,In_out_List] = simulate_diffusion_distribution(Cells,World.Voxel,Diff);
    
%% Display some diffusion molecules inside the cells
figure
DisplayStruct_ToolBox.Cells_3D(Cells(1:500,:,:,:,:,:,:,:,:),World.Voxel);
DisplayStruct_ToolBox.Water_3D(Water(1:25:end,:,:));
view([20 50])

figure
if size(Cells,2)>7
    DisplayStruct_ToolBox.Cells_Poly(Cells,World.Voxel)
else
    DisplayStruct_ToolBox.Cells(Cells,World.Voxel)
end
DisplayStruct_ToolBox.Water(Water(1:25:end,:,:));
axis([0 (World.Voxel(1)) 0 (World.Voxel(2))])


%% Simulate the MR diffusion experiments

MD=[];
FA=[];
ADC=[];
for cpt_g=1:1:Seq.Waveform
    [grad m1]= diff_grad_gen(Seq.Dn(cpt_g,:),Seq.Gn(cpt_g,:), Seq.Bval(1), Seq.RampTime,1) ;              % mT/m (us) & T.us² /m 
%    [grad m1]= diff_grad_gen([8000 0 0 4000 8000 0 0],[1 0 0 0 -1 0 0], Seq.Bval, 500,1) ;              % mT/m (us) & T.us² /m 
    grad=grad.*10^3;
    Att_approx(cpt_g,:) = DiffSim_ToolBox.Simulate_MR(Water,Seq,grad,World); % [Waveform Direction] % Grad in [T/m], dimension : NPoints x Waveform with dT corresponding to World.dT  

    % Calculate the MD/FA from the diffusion tensor 

    ADC(cpt_g,:)= DiffSim_ToolBox.ADC_dir(abs(squeeze(Att_approx(cpt_g,:))),Seq.Bval);
    ADC_mean(cpt_g)= DiffSim_ToolBox.ADC_mean(abs(squeeze(Att_approx(cpt_g,:))),Seq.Bval);
    
    [~, EigValue(cpt_g,:),~,MD(cpt_g),FA(cpt_g),~]= DiffSim_ToolBox.Tensor(abs(squeeze(Att_approx(cpt_g,:))),Seq.Dir,Seq.Bval);
    
    ADC_ref(cpt_g,:) = DiffSim_ToolBox.ADC_Water_dir(Water(:,:,1:size(grad,1)),Seq.Dir,World.dT);
    ADC_mean_ref(cpt_g) = DiffSim_ToolBox.ADC_Water_mean(Water(:,:,1:size(grad,1)),World.dT);
    EigValue_ref(cpt_g,:) = DiffSim_ToolBox.EigenValue_Water(Water(:,:,1:size(grad,1)),World.dT);
    FA_ref(cpt_g) = DiffSim_ToolBox.FA_Water(Water(:,:,1:size(grad,1)),World.dT);
end

