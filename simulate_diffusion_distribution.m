
 function [txq,In_out_List] = simulate_diffusion_distribution(Cells,Voxel_dim,Diff)

 disp(['Simulate diffusion distribution ']);
 tic
 
%% Generate Diffusion Distribution
In_out_List=[];

for cpt_dim=1:1:Diff.dim
   xt(:,cpt_dim,1)= rand(Diff.N,1)*Voxel_dim(cpt_dim);  % Uniform distribution along the voxel
end

h = waitbar(0,['Diffusion distribution ']);
[In_out_List(:,1)] = Collision_ToolBox.Collision_Detection(Cells, squeeze(xt(:,:,1)));

for cpt_t=2:1:(Diff.dur+1)
    
    % Diffusion-Collision solver - final positions
    % Reroll Diffusion until we don't have collision issues anymore
    List_coll=[1:1:Diff.N];
    
    while List_coll
        % Extracellular
        List_Extra = find(In_out_List(:,cpt_t-1));
        List_Extra=List_Extra(ismember(List_Extra, List_coll));
        tmp_water=[];
        tmp_water = DiffSim_ToolBox.diffusion_step(length(List_Extra),Diff.D,Diff.dim,Diff.dT);
        xt(List_Extra,:,cpt_t)=xt(List_Extra,:,cpt_t-1)+tmp_water;

        % Intracellular
        List_Intra = find(~In_out_List(:,cpt_t-1));
        List_Intra=List_Intra(ismember(List_Intra, List_coll));
        tmp_water=[];
        tmp_water= DiffSim_ToolBox.diffusion_step(length(List_Intra),Diff.Din,Diff.dim,Diff.dT);
        xt(List_Intra,:,cpt_t)=xt(List_Intra,:,cpt_t-1)+tmp_water;

        % Collision detection
        [tmp_List] = Collision_ToolBox.Collision_Detection(Cells, squeeze(xt(:,:,cpt_t)));
        In_out_List(:,cpt_t)=tmp_List; % New list of locations
        
        % Diffusion-Collision solver with or without Permability
        [List_coll] = Collision_ToolBox.Permeability(In_out_List(:,cpt_t-1),In_out_List(:,cpt_t),Diff.Perma);
    end
    

    waitbar(cpt_t/ Diff.dur ,h);
end

% reinterpolate to a dt of 1e-6
%txq=xt;
tt=[0:Diff.dT:(Diff.dur*Diff.dT)];
tq=[0:1e-6:(Diff.dur*Diff.dT)-1e-6];

% 
tic
for cpt_mol=1:1:size(xt,1)  
    txq(cpt_mol,1,:) = interp1(tt,squeeze(xt(cpt_mol,1,:)),tq);
    txq(cpt_mol,2,:) = interp1(tt,squeeze(xt(cpt_mol,2,:)),tq);
    txq(cpt_mol,3,:) = interp1(tt,squeeze(xt(cpt_mol,3,:)),tq);
end
toc

% tic
% for cpt_mol=1:1:size(xt,1) 
%     
%     
%     for cpt_dim=1:1:3
%         ttxq=[];
%         for cpt_time=2:1:Diff.dur
%             ttxq=[ttxq linspace(xt(cpt_mol,cpt_dim,cpt_time-1),xt(cpt_mol,cpt_dim,cpt_time),round(Diff.dT/1e-6))];
%         end
%         txq2(cpt_mol,cpt_dim,:)=ttxq;
%     end
% end
% toc

close(h)

 end