warning off;
clear all;

World.Gamma= 267.513*10^6;                                                 % rad/(s . Tesla )
World.dT= 1e-6;    

%%%%% Voxel parameters %%%%%%
World.Voxel=[0.2 0.2 0.2].*1e-3;                                           % mm -> m
World.Voxel_vol=World.Voxel(1)*World.Voxel(2)*World.Voxel(3);              % Volume in m3
World.Voxel_sur=World.Voxel(1)*World.Voxel(2);                             % Surface in m2
                                                                           % us -> s
%%%%% Diffusion parameters %%%%%%
Diff.N=1000;                                                               % Number of particles
Diff.dim=3;                                                                % Diff parameters                                                             
Diff.D=2.2*(10^-3)*(10^-6);                                                % Coefficient of diffusion extravascular mm²/s -> m²/s
Diff.Din=Diff.D*2/3;                                                       % Coefficient of diffusion intravascular mm²/s -> m²/s
Diff.V=[0 0 0];                                                            % Linear Velocity along m/s 
Diff.VGrad=[0 0 0];                                                        % Velocity Grad along x m/s 
Diff.Perma=0;
Diff.dT= 1e-6; 
                                                     
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

% Take the waveform with the maximum duration to simulate diffusion
Diff.dur= ceil((max(sum(Seq.Dn(:,:),2))*World.dT)/(Diff.dT));

%% Load a cell structure and create a mask from a H5 file 
[File.name File.path]=uigetfile('.h5');
cd(File.path);
Cells_mask = squeeze(h5read(File.name,'/exported_data'));
Cells_mask(Cells_mask==2)=0; % 2 Extracellular
Cells_mask(Cells_mask==1)=1; % 1 Intracellular
Cells_mask=logical(Cells_mask); % faster computing and memory storage
World.Resolution=World.Voxel./size(Cells_mask);

%% Calculate and estimate ECV as a checkup
World.ECV_estimate=Cellule_ToolBox.Calculate_ECV_Mask(Cells_mask,World.Voxel);
World.ECV_calculate=1-sum(Cells_mask(:))/(size(Cells_mask,1)*size(Cells_mask,2)*size(Cells_mask,3));
disp(['Calculated ECV ' num2str(World.ECV_calculate) ' Estimated ECV ' num2str(World.ECV_estimate)]);

%% Simulation diffusion inside and outside the cells
[Water,In_out_List] = simulate_diffusion_distribution_mask(Cells_mask,World.Voxel,Diff);
    
%% Display some diffusion molecules inside the cells
figure
DisplayStruct_ToolBox.Cells_Mask(Cells_mask(:,:,1));
%DisplayStruct_ToolBox.Water_Mask(Water(squeeze(~logical(In_out_List(:,1))),:,:),World.Resolution);

%% Simulate the MR diffusion experiments
% Init
Att_approx=[];
MD=[];FA=[];EigValue=[];ADC=[];ADC_mean=[];
ADC_ref=[];ADC_mean_ref=[];EigValue_ref=[];FA_ref=[];

% Per waveform
for cpt_g=1:1:Seq.Waveform
    tic
    
    % Generate the waveform based on the timing furnished before
    [grad m1]= diff_grad_gen(Seq.Dn(cpt_g,:),Seq.Gn(cpt_g,:), Seq.Bval, Seq.RampTime,1) ;              % mT/m (us) & T.us² /m    
    grad=grad.*10^3;
    
    % Simuate the MR experiement 
    Att_approx(cpt_g,:,:) = DiffSim_ToolBox.Simulate_MR(Water,Seq,grad,World) ; % [Waveform Direction] % Grad in [T/m], dimension : NPoints x Waveform with dT corresponding to World.dT  

    % Calculate the ADC and tensor paramaters from simulation for each b-values  
    for cpt_b=1:1:length(Seq.Bval)
        ADC(cpt_g,cpt_b,:)      = DiffSim_ToolBox.ADC_dir(abs(squeeze(Att_approx(cpt_g,cpt_b,:))),Seq.Bval(cpt_b));
        ADC_mean(cpt_g,cpt_b)   = DiffSim_ToolBox.ADC_mean(abs(squeeze(Att_approx(cpt_g,cpt_b,:))),Seq.Bval(cpt_b));
        [~, EigValue(cpt_g,cpt_b,:),~,MD(cpt_g,cpt_b),FA(cpt_g,cpt_b),~]= DiffSim_ToolBox.Tensor(abs(squeeze(Att_approx(cpt_g,cpt_b,:))),Seq.Dir,Seq.Bval(cpt_b));     
    end
    
    % Calculate the ground truth ADC and tensor paramaters from water displacement at TE 
    ADC_ref(cpt_g,:)        = DiffSim_ToolBox.ADC_Water_dir(Water(:,:,1:size(grad,1)),Seq.Dir,World.dT);
    ADC_mean_ref(cpt_g)     = DiffSim_ToolBox.ADC_Water_mean(Water(:,:,1:size(grad,1)),World.dT);
    EigValue_ref(cpt_g,:)   = DiffSim_ToolBox.EigenValue_Water(Water(:,:,1:size(grad,1)),World.dT);
    FA_ref(cpt_g)           = DiffSim_ToolBox.FA_Water(Water(:,:,1:size(grad,1)),World.dT);
    
    disp(['Simulate MR diffusion sequence ' num2str(cpt_g) ' in ' num2str(toc) 'ms']);
end

