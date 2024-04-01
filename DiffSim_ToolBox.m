classdef DiffSim_ToolBox

    methods(Static)
        
        %% General function to calculate the Tensor out of signal attenuation
        function [T6, EigValue,EigVector,MD,FA,Trace_DTI]= Tensor(dwi,dir,bval)

            Tensor=[];
            EigValue=[];
            EigVector=[];
            MD=[];
            FA=[];
            Trace_DTI=[];
            tmpInput=[];  

              tmpInput(1,1,1)=1;
              tmpInput(1,1,(2:size(dir,1)+1))=dwi;
              [T6,EigVector,EigValue, Mat] = DiffSim_ToolBox.Tensor_local(tmpInput,dir,bval);

              [MD, FA, Trace] =  DiffSim_ToolBox.Maps_local(EigValue); 
              
        end

        function [T6, EigValue,EigVector,MD,FA,Trace_DTI]= Tensor_Water(Water,dT)

            Tensor=[];
            EigValue=[];
            EigVector=[];
            MD=[];
            FA=[];
            Trace_DTI=[];
            tmpInput=[];  
            
            % b-value and dir are arbitrary to give us enought precision.
            bval=1000;
            dir=[0.831862183551238,-0.978020640898505,-0.746594402068549,0.618503898022940,0.0357105488849347,0.0690933229952567,0.285311347832409,0.541890310423585,0.696013916160111,-0.586381611851074,-0.361852167119390,-0.107233806522707;
                 -0.474849677041022,-0.117331253822020,-0.566355417065775,-0.0500144660490100,-0.296253369933249,-0.451075167295507,-0.874813688131557,-0.803757865561315,-0.264300268733346,-0.487477487497046,-0.922081582695813,-0.849061587106308;
                 -0.287268170935018,-0.172362813165683,0.349053404091507,-0.784188429834895,0.954441552711641,-0.889807334332884,-0.391533602085168,0.245617655503534,0.667615050863187,-0.646933288596948,-0.137217684961672,0.517296200740132]';
            
             % create synthetic dwi values.
             for cpt_dwi=1:1:size(dir,1)
                   [tmp_adc] = DiffSim_ToolBox.ADC_Water_dir(Water,squeeze(dir(cpt_dwi,:)),dT);        
                   dwi(cpt_dwi)=exp(-bval*tmp_adc);
             end
              tmpInput(1,1,1)=1;
              tmpInput(1,1,(2:size(dir,1)+1))=dwi;
              [T6,EigVector,EigValue, Mat] = DiffSim_ToolBox.Tensor_local(tmpInput,dir,bval);

              [MD, FA, Trace] =  DiffSim_ToolBox.Maps_local(EigValue); 
              
        end
        
        %% Local function to calculte the fit the tensor
        function [Tensor,EigVector,EigValue, Mat2] = Tensor_local(slice,dir,bvalue)

        % slice : x y [B0 B1-dir1 B2-dir2...]
        % dir : x1 y1 z1
        %       x2 y2 z2
        %        .  .  .
        %        .  .  .
        %
        %  Tensor : x y [6]
        %  EigVector : x y [3 coor][3 num]
        %  EigValue : x y [3]


        nb_dir=size(dir,1);
        H = zeros(nb_dir,6);

        for i=1:nb_dir
           H(i,:)=[dir(i,1)*dir(i,1),dir(i,2)*dir(i,2),dir(i,3)*dir(i,3),2*dir(i,1)*dir(i,2),2*dir(i,1)*dir(i,3),2*dir(i,2)*dir(i,3)];
        end
        [U,S,V] = svd(H,0);      
        H_inv=V*inv(S)*U';                 %H is a 6*30 matrix



        Tensor=zeros(size(slice,1),size(slice,2),6);
        EigVector=zeros(size(slice,1),size(slice,2),3,3);
        EigValue= zeros(size(slice,1),size(slice,2),3);
        Mat=zeros(size(slice,1),size(slice,2),1,9);
        Mat2=zeros(size(slice,1),size(slice,2),3,3);
        for y=1:1:size(slice,1)
            for x=1:1:size(slice,2)
                Y=[];

                if slice(y,x,1)~=0
                    for z=1:1:nb_dir
                        if double(slice(y,x,z+1))~=0
                            Y=[Y;log(double(slice(y,x,1))/double(slice(y,x,z+1)))/bvalue]; % Y = [log(S0/S1), log(S0/S2), log(S0,S3)....]
                        else
                            Y=[Y;0]; % Y = [log(S0/S1), log(S0/S2), log(S0,S3)....]
                        end
                    end
                    Tensor(y,x,:)=H_inv*Y;
                    Mat=[Tensor(y,x,1),Tensor(y,x,4),Tensor(y,x,5);Tensor(y,x,4),Tensor(y,x,2),Tensor(y,x,6);Tensor(y,x,5),Tensor(y,x,6),Tensor(y,x,3)];
                     if mean(mean(isnan(Mat)))>0 || mean(mean(isinf(Mat)))>0
                        Mat=zeros(3,3); 
                    end
                    Mat2(y,x,:,:)=Mat;

                    [Vect,Diag]=eig(Mat);

                    EigValue(y,x,:)=[abs(Diag(1,1)),abs(Diag(2,2)),abs(Diag(3,3))];

                    [t,index]=sort(EigValue(y,x,:),'descend');
                    EigValue(y,x,:)=EigValue(y,x,index);
                    EigVector(y,x,:,:)=Vect(:,index);
                else
                    Mat2(y,x,:,:)=0;
                    Tensor(y,x,:)=0;
                    EigValue(y,x,:)=0;
                    EigVector(y,x,:,:)=0;
                end

            end
        end
        end

        %% Local function to calculte the invarient of the tensor
        function [MD, FA, TRACE, I1, I2, I3] = Maps_local(EigValue)

        MD=zeros(size(EigValue,1),size(EigValue,2));
        FA=zeros(size(EigValue,1),size(EigValue,2));
        TRACE=zeros(size(EigValue,1),size(EigValue,2));
        I1=zeros(size(EigValue,1),size(EigValue,2));
        I2=zeros(size(EigValue,1),size(EigValue,2));
        I3=zeros(size(EigValue,1),size(EigValue,2));

        for y=1:1:size(EigValue,1)
            for x=1:1:size(EigValue,2)
            TRACE(y,x)=EigValue(y,x,1)+EigValue(y,x,2)+EigValue(y,x,3);
            MD(y,x)=(EigValue(y,x,1)+EigValue(y,x,2)+EigValue(y,x,3))/3;
            FA(y,x)=sqrt((3*((EigValue(y,x,1)-MD(y,x))^2+(EigValue(y,x,2)-MD(y,x))^2+(EigValue(y,x,3)-MD(y,x))^2))/(2*(EigValue(y,x,1)^2+EigValue(y,x,2)^2+EigValue(y,x,3)^2)));

            I1(y,x)=EigValue(y,x,1)+EigValue(y,x,2)+EigValue(y,x,3);
            I2(y,x)=EigValue(y,x,1)*EigValue(y,x,2)+EigValue(y,x,2)*EigValue(y,x,3)+EigValue(y,x,1)*EigValue(y,x,3);
            I3(y,x)=EigValue(y,x,1)*EigValue(y,x,2)*EigValue(y,x,3);

            end
        end
        end
        
        %% Function to calculate the mean ADC from several signal attenuations
        function [ADC] = ADC_mean(dwi,bval)
               ADC=log(mean(dwi))/-bval;                      
        end
         
        %% Function to calculate the ADC per direction from signal attenuations
        function [ADC] = ADC_dir(dwi,bval)
                for cpt_dir=1:1:length(dwi)
                   ADC(cpt_dir)=log(dwi(cpt_dir))/-bval; 
                end
         end
         
        %% Calculate the coefficient of diffusion from a dispacement distribution
        function [ADC] = ADC_Water_mean(Water,dT) 
            
             % sum of the mean square of distance
            dx=abs(Water(:,:,end)-Water(:,:,1));
            % sum of x y z
            ssdx=squeeze(sum(dx.*dx,2));
            % Coefficient of diffusion mean square distance over dT in 3d, unit in mm2/s
            ADC=1e6*mean(ssdx)/(6*size(Water,3)*dT);
         end
        
        %% Calculate the coefficient of diffusion from a dispacement distribution per direction
        function [ADC] = ADC_Water_dir(Water,dir,dT) 
            
             % sum of the mean square of distance
            dx=(Water(:,:,end)-Water(:,:,1));   
            for cpt_dir=1:1:size(dir,1)
               
                % Make the projection on the corresponding axis 
                
               dx_dir= Vector_ToolBox.Projection_vect_n(dx, dir(cpt_dir,:));
              %  dx_dir=dx.*(repmat((dir(cpt_dir,:)),size(Water,1),1));
              %  dx_dir3=dx.*(repmat(abs(dir(cpt_dir,:)),size(Water,1),1));
                % Take the norm of the vectors
                %dx_norm=sqrt(dx_dir(:,1).*dx_dir(:,1)+dx_dir(:,2).*dx_dir(:,2)+dx_dir(:,3).*dx_dir(:,3));
                dx_norm=Vector_ToolBox.Norm_vect_n(dx_dir);
                % sum of x y z
                ssdx2=squeeze(mean(dx_norm.*dx_norm));
                % Coeffiecient of diffusion mean square distance over dT in 3d, unit in mm2/s
                % NB: the true equation to calculate the diffusion from a
                % displacement distribution is D= <X^2>/ ( 2^Dim * t )
                % By calculating my diffusion coeffcient in 3 Dimension, 
                % my ADC is too small by a factor 3
                % I believe that because in MR diffusion we only calculate
                % the ADC based on a projection of diffusion gradient
                % the correct dimension is actually only 1
                % To be verified... 
                ADC(cpt_dir)=1e6*ssdx2/(2*size(Water,3)*dT);
            end
        end
        
        %% Calculate the coefficient of diffusion from a dispacement distribution in X Y Z
        function [ADC] = ADC_Water_xyz(Water,dT)  
            dir=[[1 0 0];[0 1 0];[0 0 1];];
            [ADC] = DiffSim_ToolBox.ADC_Water_dir(Water,dir,dT) ;
        end
        
        %% Calculate the theorical FA considering XYZ direction as principal Eigen vector directions
        function [FA] = FA_Water(Water,dT)
            [ADC_xyz] = DiffSim_ToolBox.ADC_Water_xyz(Water,dT); 
            FA= sqrt(1/2)* sqrt( (ADC_xyz(1)-ADC_xyz(2)).^2 +(ADC_xyz(2)-ADC_xyz(3)).^2 +(ADC_xyz(3)-ADC_xyz(1)).^2)./sqrt( ADC_xyz(1).^2 +ADC_xyz(2).^2 +ADC_xyz(3).^2 ) ;
        end
        
        %% Calculate the theorical EigenValues considering XYZ direction as principal Eigen vector directions
        function [Eigen] = EigenValue_Water(Water,dT)
            [ADC_xyz] = DiffSim_ToolBox.ADC_Water_xyz(Water,dT); 
            Eigen=sort(ADC_xyz,'descend');
        end
        
        %% Simulate a diffusion step  
        function [Dstep] = diffusion_step(N,D,Dim,dT)  
   
        KD = sqrt(D * 3 * Dim * dT);
        Dvar=KD * randn(N,1);
         
        VectorDir=randn(N,3);
        VectorDir=VectorDir./repmat(Vector_ToolBox.Norm_vect_n(VectorDir),1,3);
        
        Dstep=VectorDir.*Dvar;
      
        end
        
        %% Simulate MR from a vector gradient file and water 
        function [Att_approx, Phi_Dist]= Simulate_MR(Water,Seq,Grad)
        %%% Simulation gradients
        Phi_Dist=[];
        Nb_Water=size(Water,1);

            %%% Loop over Direction %%%
            for cpt_dir=1:1:size(Seq.Dir,1)
                %%% Loop over Waveform and or B-values %%%
                for cpt_g=1:1:size(Grad,2)          
                     %%% Intregral Diffusion + gradients  %%%           
                     phi=zeros(Nb_Water,1);   
                    % tmp=[];
                     for cpt=2:1:size(Grad,1)
                        phi(:)=phi(:)+ nansum( Seq.Gamma*Seq.dT*( repmat(Grad(cpt,cpt_g),Nb_Water,3).*repmat(Seq.Dir(cpt_dir,:),Nb_Water,1).*squeeze(Water(:,:,cpt))),2);
                        %tmp(cpt,:)=phi(:);
                     end

                    %%% Phase Distribution
                    Phi_Dist(cpt_g,cpt_dir,:)=phi(:);
                    %%% MR Signal attenuation
                    Att_approx(cpt_g,cpt_dir)=nanmean(cos(phi(:))+i*sin(phi(:)));  
                    %%%
                end
            end
        end
    
        %% Simulate MR from a vector gradient file and water 
        function [Att_approx, Phi_Dist]= Simulate_MR_optim(Water,Seq,Grad,World)
        %%% Simulation gradients
        Phi_Dist=[];
        Nb_Water=size(Water,1);
        L_Grad=size(Grad,1);
        %tmp_w=permute(squeeze(Water(:,:,1:L_Grad)),[2 3 1]);
            %%% Loop over Direction %%%    
                for cpt_g=1:1:size(Grad,2)          
                     %%% Intregral Diffusion + gradients  %%%           
                     phi=zeros(Nb_Water,1);   
%                      for cpt=2:1:size(Grad,1)
%                         phi(:)=phi(:)+ nansum( World.Gamma *World.dT*( repmat(Grad(cpt,cpt_g),Nb_Water,3).*repmat(Seq.Dir(cpt_dir,:),Nb_Water,1).*squeeze(Water(:,:,cpt))),2);
%                      end
                       
%                        tic 
%                       for cpt=1:1:Nb_Water
%                         tmp= cumsum( nansum(World.Gamma *World.dT*( repmat(Grad(:,cpt_g),1,3).*Seq.Dir(cpt_dir,:).*squeeze(Water(cpt,:,1:L_Grad))')));
%                         phi(cpt)=tmp(end);
%                       end
%                        toc
                       
                       tmp_g=squeeze(Grad(:,cpt_g));
                       for cpt=1:1:Nb_Water
                        tmp_x= World.Gamma *World.dT*(tmp_g.*squeeze(Water(cpt,1,1:L_Grad)));
                        tmp_y= World.Gamma *World.dT*(tmp_g.*squeeze(Water(cpt,2,1:L_Grad)));
                        tmp_z= World.Gamma *World.dT*(tmp_g.*squeeze(Water(cpt,3,1:L_Grad)));
                         
                        
                        for cpt_dir=1:1:size(Seq.Dir,1) 
                             tmp_x2=cumsum(tmp_x.*Seq.Dir(cpt_dir,1));
                             tmp_y2=cumsum(tmp_y.*Seq.Dir(cpt_dir,2));
                             tmp_z2=cumsum(tmp_z.*Seq.Dir(cpt_dir,3));
                             Phi_Dist(cpt_g,cpt_dir,cpt)=tmp_x2(end)+tmp_y2(end)+tmp_z2(end);
                             
                        end
                       end
                       
                       for cpt_dir=1:1:size(Seq.Dir,1) 
                         Att_approx(cpt_g,cpt_dir)=nanmean(cos( Phi_Dist(cpt_g,cpt_dir,:))+i*sin( Phi_Dist(cpt_g,cpt_dir,:)));  
                       end
                    %%% Phase Distribution
                    %%% MR Signal attenuation
                    
                    
                    %%%
            end
        end
        
        %% Old simulation
        function [txq,In_out_List] = SimDiff(Cells,Voxel_dim,Diff)

%                 disp(['Simulate diffusion distribution ']);
%                 tic

                %% Generate Diffusion Distribution
                In_out_List=[];

                for cpt_dim=1:1:Diff.dim
                   xt(:,cpt_dim,1)= rand(Diff.N,1)*Voxel_dim(cpt_dim);  % Uniform distribution along the voxel
                end

                h = waitbar(0,['Diffusion distribution ']);
                [In_out_List(:,1)] = Collision_ToolBox2.Collision_Detection(Cells, squeeze(xt(:,:,1)));

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
                        [tmp_List] = Collision_ToolBox2.Collision_Detection(Cells, squeeze(xt(:,:,cpt_t)));
                        In_out_List(:,cpt_t)=tmp_List; % New list of locations

                        % Diffusion-Collision solver with or without Permability
                        [List_coll] = Collision_ToolBox2.Permeability(In_out_List(:,cpt_t-1),In_out_List(:,cpt_t),Diff.Perma);
                    end


                    waitbar(cpt_t/ Diff.dur ,h);
                end

                % reinterpolate to a dt of 1e-6
                %txq=xt;
                tt=[0:Diff.dT:(Diff.dur*Diff.dT)];
                tq=[0:1e-6:(Diff.dur*Diff.dT)-1e-6];

                % 
%                 tic
                for cpt_mol=1:1:size(xt,1)  
                    txq(cpt_mol,1,:) = interp1(tt,squeeze(xt(cpt_mol,1,:)),tq);
                    txq(cpt_mol,2,:) = interp1(tt,squeeze(xt(cpt_mol,2,:)),tq);
                    txq(cpt_mol,3,:) = interp1(tt,squeeze(xt(cpt_mol,3,:)),tq);
                end
               % toc

                 close(h)

        end

        %% Old simulation mask
        function [txq,In_out_List] = SimDiff_mask(Cells_mask,Voxel_dim,Diff)

        % disp(['Simulate diffusion distribution from a mask']);
        % tic

        %% Generate Diffusion Distribution
        In_out_List=[];

        for cpt_dim=1:1:Diff.dim
           xt(:,cpt_dim,1)= rand(Diff.N,1)*Voxel_dim(cpt_dim);  % Uniform distribution along the voxel
        end

        h = waitbar(0,['Diffusion distribution ']);

        Resolution=Voxel_dim./size(Cells_mask);


        [In_out_List(:,1)] = Collision_ToolBox2.Collision_Detection_Mask(Cells_mask, squeeze(xt(:,:,1)),Resolution);

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
                [tmp_List] = Collision_ToolBox2.Collision_Detection_Mask(Cells_mask, squeeze(xt(:,:,cpt_t)),Resolution);
                In_out_List(:,cpt_t)=tmp_List; % New list of locations

                % Diffusion-Collision solver with or without Permability
                [List_coll] = Collision_ToolBox2.Permeability(In_out_List(:,cpt_t-1),In_out_List(:,cpt_t),Diff.Perma);

        %         figure,
        %         hold on
        %         imagesc(Cells_mask(:,:,1))
        %         %scatter(xt(:,1,cpt_t)./Resolution(1),xt(:,2,cpt_t)./Resolution(1),40,'g')
        %         plot(squeeze(xt(List_coll,2,1:cpt_t)./Resolution(1))',squeeze(xt(List_coll,1,1:cpt_t)./Resolution(1))','g')
        %         scatter(xt(List_coll,2,cpt_t)./Resolution(1),xt(List_coll,1,cpt_t)./Resolution(1),40,'r','filled')

            end


            waitbar(cpt_t/ Diff.dur ,h);
        end

        txq=xt;

        close(h)
        end
        
        %% Old simulation mask
        function [txq,In_out_List] = Water_distribution_mask_old(Cells_mask,Voxel_dim,Diff)

       %  disp(['Simulate diffusion distribution from a mask']);
       %  tic

        %% Generate Diffusion Distribution
        In_out_List=[];

        for cpt_dim=1:1:Diff.dim
           xt(:,cpt_dim,1)= rand(Diff.N,1)*Voxel_dim(cpt_dim);  % Uniform distribution along the voxel
        end

        h = waitbar(0,['Diffusion distribution ']);

        Resolution=Voxel_dim./size(Cells_mask);

        

        [In_out_List(:,1)] = Collision_ToolBox2.Collision_Detection_Mask(Cells_mask, squeeze(xt(:,:,1)),Resolution);

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
                [tmp_List] = Collision_ToolBox2.Collision_Detection_Mask(Cells_mask, squeeze(xt(:,:,cpt_t)),Resolution);
                In_out_List(:,cpt_t)=tmp_List; % New list of locations

                % Diffusion-Collision solver with or without Permability
                [List_coll] = Collision_ToolBox2.Permeability(In_out_List(:,cpt_t-1),In_out_List(:,cpt_t),Diff.Perma);

        %         figure,
        %         hold on
        %         imagesc(Cells_mask(:,:,1))
        %         %scatter(xt(:,1,cpt_t)./Resolution(1),xt(:,2,cpt_t)./Resolution(1),40,'g')
        %         plot(squeeze(xt(List_coll,2,1:cpt_t)./Resolution(1))',squeeze(xt(List_coll,1,1:cpt_t)./Resolution(1))','g')
        %         scatter(xt(List_coll,2,cpt_t)./Resolution(1),xt(List_coll,1,cpt_t)./Resolution(1),40,'r','filled')

            end


            waitbar(cpt_t/ Diff.dur ,h);
        end

        % reinterpolate to a dt of 1e-6
        txq=xt;
        % tt=[0:Diff.dT:(Diff.dur*Diff.dT)];
        % tq=[0:1e-6:(Diff.dur*Diff.dT)-1e-6];
        % 
        % 
        % tic
        % for cpt_mol=1:1:size(xt,1)  
        %     txq(cpt_mol,1,:) = interp1(tt,squeeze(xt(cpt_mol,1,:)),tq);
        %     txq(cpt_mol,2,:) = interp1(tt,squeeze(xt(cpt_mol,2,:)),tq);
        %     txq(cpt_mol,3,:) = interp1(tt,squeeze(xt(cpt_mol,3,:)),tq);
        % end
       % toc

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
        
        %% Simulation of water inside the mask 
        function [txq,In_out_List] = Water_distribution_mask(xt,Cells_mask,Voxel_dim,Diff)

         %disp(['Simulate diffusion distribution from a mask']);
        % tic
         Diff.Res=Voxel_dim./size(Cells_mask);
        % Generate Diffusion Distribution
        In_out_List=[];

%         for cpt_dim=1:1:Diff.dim
%            xt(:,cpt_dim,1)= rand(Diff.N,1)*Voxel_dim(cpt_dim);  % Uniform distribution along the voxel %repmat(250*Diff.Res(cpt_dim),Diff.N,1);%
%         end

        h = waitbar(0,['Diffusion distribution ']);


        xt(:,:,1)=DiffSim_ToolBox.Water_distribution_init(Voxel_dim,Cells_mask,Diff);
        [In_out_List(:,1)] = Collision_ToolBox2.Collision_Detection_Mask(Cells_mask, squeeze(xt(:,:,1)),Diff.Res);

        tmp_in_out=In_out_List(:,1);
        tmp_xt= xt(:,:,1);
        cpt_t=2;
        txq(:,:,1)=xt(:,:,1);
        for cpt_t_tmp=2:1:(Diff.dur+1)
             if Diff.Iterative
                       [tmp_xt,tmp_in_out]=DiffSim_ToolBox.Water_collision_iterative(tmp_xt,tmp_in_out,Cells_mask,Diff,1,1);
             else
                       [tmp_xt,tmp_in_out]=DiffSim_ToolBox.Water_collision(tmp_xt,tmp_in_out,Cells_mask,Diff);
             end
             % Save only the relevant diffusion displacement
             if isfield(Diff,'vMt')
                 if ~isempty(Diff.vMt)
                     Idx=find(cpt_t_tmp<Diff.vMt(:,1));
                     if (~isempty(Idx))
                         if Diff.vMt(Idx(1),2)
                             txq(:,:,cpt_t)=tmp_xt;
                             In_out_List(:,cpt_t)=tmp_in_out;
                             cpt_t=cpt_t+1;
                         end
                     else
                            txq(:,:,cpt_t)=tmp_xt;
                            In_out_List(:,cpt_t)=tmp_in_out;
                            cpt_t=cpt_t+1;
                     end
                 else
                    txq(:,:,cpt_t)=tmp_xt;
                    In_out_List(:,cpt_t)=tmp_in_out;
                    cpt_t=cpt_t+1;
                 end
             else
                 txq(:,:,cpt_t)=tmp_xt;
                 In_out_List(:,cpt_t)=tmp_in_out;
                 cpt_t=cpt_t+1;
             end
            waitbar(cpt_t_tmp/ Diff.dur ,h);
        end

        % reinterpolate to a dt of 1e-6
       % txq=xt;
       % toc;

        close(h)

        end

        %% Iterative collision solver
        function [xt2,In_out_List] = Water_collision(xt,In_out_List,Cells_mask,Diff)
            
         % Diffusion-Collision solver - final positions
            % Reroll Diffusion until we don't have collision issues anymore
            List_coll=[1:1:size(xt,1)];

            % Extracellular
            List_Extra = find(In_out_List);
            if List_Extra
                tmp_water=[];
                tmp_water = DiffSim_ToolBox.diffusion_step(length(List_Extra),Diff.D,Diff.dim,Diff.dT);
                xt2(List_Extra,:)=xt(List_Extra,:)+tmp_water;
            end
            
             % Intracellular
             List_Intra = find(~In_out_List); 
             if List_Intra
                tmp_water=[];
                tmp_water= DiffSim_ToolBox.diffusion_step(length(List_Intra),Diff.Din,Diff.dim,Diff.dT);
                xt2(List_Intra,:)=xt(List_Intra,:)+tmp_water;
             end
            % Collision detection
            [tmp_List] = Collision_ToolBox2.Collision_Detection_Mask(Cells_mask, squeeze(xt2),Diff.Res);
            In_out_List_new=tmp_List; % New list of locations

            % Diffusion-Collision solver with or without Permability only
            [List_coll] = Collision_ToolBox2.Permeability(In_out_List,In_out_List_new,Diff.Perma);
            
            % We reroll smaller case until we have no more collision. 
              while List_coll
                    % Extracellular
                    List_Extra = find(In_out_List);
                    List_Extra=List_Extra(ismember(List_Extra, List_coll));
                    tmp_water=[];
                    tmp_water = DiffSim_ToolBox.diffusion_step(length(List_Extra),Diff.D,Diff.dim,Diff.dT);
                    xt2(List_Extra,:)=xt(List_Extra,:)+tmp_water;

                    % Intracellular
                    List_Intra = find(~In_out_List);
                    List_Intra=List_Intra(ismember(List_Intra, List_coll));
                    tmp_water=[];
                    tmp_water= DiffSim_ToolBox.diffusion_step(length(List_Intra),Diff.Din,Diff.dim,Diff.dT);
                    xt2(List_Intra,:)=xt(List_Intra,:)+tmp_water;

                    % Collision detection
                    [tmp_List] = Collision_ToolBox2.Collision_Detection_Mask(Cells_mask, squeeze(xt2),Diff.Res);
                    In_out_List_new=tmp_List; % New list of locations

                    % Diffusion-Collision solver with or without Permability
                    [List_coll] = Collision_ToolBox2.Permeability(In_out_List,In_out_List_new,Diff.Perma);    
              end
                 
            
            
        end
        
         %% Iterative collision solver
        function [xt2,In_out_List] = Water_collision_iterative(xt,In_out_List,Cells_mask,Diff,scale,iter)
            
            % Diffusion-Collision solver - final positions
            % Reroll Diffusion until we don't have collision issues anymore
            List_coll=[1:1:size(xt,1)];

            % Extracellular
            List_Extra = find(In_out_List);
            if List_Extra
                tmp_water=[];
                tmp_water = DiffSim_ToolBox.diffusion_step(length(List_Extra),Diff.D,Diff.dim,Diff.dT/scale);
                xt2(List_Extra,:)=xt(List_Extra,:)+tmp_water;
            end
            
             % Intracellular
             List_Intra = find(~In_out_List); 
             if List_Intra
                tmp_water=[];
                tmp_water= DiffSim_ToolBox.diffusion_step(length(List_Intra),Diff.Din,Diff.dim,Diff.dT/scale);
                xt2(List_Intra,:)=xt(List_Intra,:)+tmp_water;
             end
            % Collision detection
            [tmp_List] = Collision_ToolBox2.Collision_Detection_Mask(Cells_mask, squeeze(xt2),Diff.Res);
            In_out_List_new=tmp_List; % New list of locations

            % Diffusion-Collision solver with or without Permability only
            % in the first iterative test
            if scale==1
                [List_coll] = Collision_ToolBox2.Permeability(In_out_List,In_out_List_new,Diff.Perma);
            else
                [List_coll] = Collision_ToolBox2.Permeability(In_out_List,In_out_List_new,0); 
            end
             
            % We reroll smaller case until we have no more collision. 
            if List_coll
                if scale<=2             
                    [xt2(List_coll,:),~] = DiffSim_ToolBox.Water_collision_iterative(xt(List_coll,:),In_out_List(List_coll,:),Cells_mask,Diff,scale*2,1);
               
                
                else % To avoid ultra long recursive call, we impose a bottom floor
                      while List_coll
                            % Extracellular
                            List_Extra = find(In_out_List);
                            List_Extra=List_Extra(ismember(List_Extra, List_coll));
                            tmp_water=[];
                            tmp_water = DiffSim_ToolBox.diffusion_step(length(List_Extra),Diff.D,Diff.dim,Diff.dT/scale);
                            xt2(List_Extra,:)=xt(List_Extra,:)+tmp_water;

                            % Intracellular
                            List_Intra = find(~In_out_List);
                            List_Intra=List_Intra(ismember(List_Intra, List_coll));
                            tmp_water=[];
                            tmp_water= DiffSim_ToolBox.diffusion_step(length(List_Intra),Diff.Din,Diff.dim,Diff.dT/scale);
                            xt2(List_Intra,:)=xt(List_Intra,:)+tmp_water;

                            % Collision detection
                            [tmp_List] = Collision_ToolBox2.Collision_Detection_Mask(Cells_mask, squeeze(xt2),Diff.Res);
                            In_out_List_new=tmp_List; % New list of locations

                            % Diffusion-Collision solver with or without Permability
                            [List_coll] = Collision_ToolBox2.Permeability(In_out_List,In_out_List_new,0);    
                      end
                end 
            end
            
            % At this point we have run everything for dT/scale so now we
            % have to rerun the same dT/scale to get the other half of the
            % trajectory
            if scale~=1 & iter  
                    [xt2,~] = DiffSim_ToolBox.Water_collision_iterative(xt2,In_out_List,Cells_mask,Diff,scale,0);
            end 
            
        end
        
        %% Old simulation
        function [txq,In_out_List] = Water_distribution(Cells,Voxel_dim,Diff)

         %disp(['Simulate diffusion distribution ']);
        % tic

        %% Generate Diffusion Distribution
        In_out_List=[];

        for cpt_dim=1:1:Diff.dim
           xt(:,cpt_dim,1)= rand(Diff.N,1)*Voxel_dim(cpt_dim);  % Uniform distribution along the voxel
        end

        h = waitbar(0,['Diffusion distribution ']);
        [In_out_List(:,1)] = Collision_ToolBox2.Collision_Detection(Cells, squeeze(xt(:,:,1)));

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
                [tmp_List] = Collision_ToolBox2.Collision_Detection(Cells, squeeze(xt(:,:,cpt_t)));
                In_out_List(:,cpt_t)=tmp_List; % New list of locations

                % Diffusion-Collision solver with or without Permability
                [List_coll] = Collision_ToolBox2.Permeability(In_out_List(:,cpt_t-1),In_out_List(:,cpt_t),Diff.Perma);
            end


            waitbar(cpt_t/ Diff.dur ,h);
        end

        % reinterpolate to a dt of 1e-6
        %txq=xt;
        tt=[0:Diff.dT:(Diff.dur*Diff.dT)];
        tq=[0:1e-6:(Diff.dur*Diff.dT)-1e-6];

        % 
       % tic
        for cpt_mol=1:1:size(xt,1)  
            txq(cpt_mol,1,:) = interp1(tt,squeeze(xt(cpt_mol,1,:)),tq);
            txq(cpt_mol,2,:) = interp1(tt,squeeze(xt(cpt_mol,2,:)),tq);
            txq(cpt_mol,3,:) = interp1(tt,squeeze(xt(cpt_mol,3,:)),tq);
        end
       % toc

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
        
        
       
            
           %% Simulation of isotropique water 
        function [xt]= Water_distribution_iso(Voxel_dim,Diff)
           
            
            
            for cpt_dim=1:1:Diff.dim
                   xt(:,cpt_dim,1)= rand(Diff.N,1)*Voxel_dim(cpt_dim);  % Uniform distribution along the voxel %repmat(250*Diff.Res(cpt_dim),Diff.N,1);%
                 %  xt(:,cpt_dim,1)= zeros(Diff.N,1)*Voxel_dim(cpt_dim);
            end
            
            h = waitbar(0,['Diffusion distribution ']);
            for cpt_t=2:1:(Diff.dur)
                       tmp_water = DiffSim_ToolBox.diffusion_step(Diff.N,Diff.D,Diff.dim,Diff.dT);
                       xt(:,:,cpt_t)=xt(:,:,cpt_t-1)+tmp_water;
                       waitbar(cpt_t/ Diff.dur ,h);
            end
            close(h)   
        end
            
         %% Init the water in the voxel in the given compartment
        function [xt]= Water_distribution_init(Voxel_dim,Cells_mask,Diff)
             
             
             xt=[];
             
             
             for cpt_dim=1:1:Diff.dim
                   xt(:,cpt_dim,1)= rand(Diff.N,1)*Voxel_dim(cpt_dim);  % Uniform distribution along the voxel %repmat(250*Diff.Res(cpt_dim),Diff.N,1);%
             end
             
             % If we specified a compartement we are interested in
             tic
             if ~isempty(Diff.Compartement)
                 In_out_List=[];
                 [In_out_List(:,1)] = Collision_ToolBox2.Collision_Detection_Mask(Cells_mask, squeeze(xt(:,:,1)),Diff.Res); 
                 
                 List = find(In_out_List==Diff.Compartement); % List of Molecules outside the compartment of interest
                 while List  % We reroll the molecules until everything is in the compartement of interest
                     for cpt_dim=1:1:Diff.dim
                           xt(List,cpt_dim,1)= rand(length(List),1)*Voxel_dim(cpt_dim);  % Uniform distribution along the voxel %repmat(250*Diff.Res(cpt_dim),Diff.N,1);%
                     end
                     [In_out_List(:,1)] = Collision_ToolBox2.Collision_Detection_Mask(Cells_mask, squeeze(xt(:,:,1)),Diff.Res); 
                	 List = find(In_out_List~=Diff.Compartement);
                 end
             end
             toc
        end
        
        function MRResults=local_MR_distribution(MR_struct,Water) 
         for cpt_MR=1:1:length(MR_struct.struct)
              
                MR_params=MR_struct.struct{cpt_MR};
                
                   
                    %disp([' Simulating MR distribution ' num2str(cpt_MR) ' with Diff distribution ' ]);
                    
                    %%% 0) Init
                    MRResults=[];
                                                          
                    for cpt_g=1:1:MR_params.Waveform            % Loop over waveforms
                        for cpt_b=1:1:length(MR_params.Bval)    % Loop over b-values
                           
                            %%% 1) Generate the MR waveform
                            [G, t]= Waveform_ToolBox.Factory(MR_params.Dn(cpt_g,:),MR_params.Gn(cpt_g,:),MR_params.RampTime,MR_params.dT*1e6,MR_params.Bval);
                            G=G*10^-3; % Convert from mT/m to T/m

                            %%% 2) Simulate the MR experiement 
                            [MRResults.Att_approx(cpt_g,cpt_b,:,:), MRResults.Phi_dist(cpt_g,cpt_b,:,:,:)]= DiffSim_ToolBox.Simulate_MR(Water,MR_params,G) ; % [Waveform Direction] % Grad in [T/m], dimension : NPoints x Waveform with dT corresponding to World.dT  

                            %%% 3) Calculate the ADC and tensor paramaters from the MR simulation and collect the results      
                             MRResults.ADC(cpt_g,cpt_b)    = DiffSim_ToolBox.ADC_mean(abs(squeeze(MRResults.Att_approx(cpt_g,cpt_b,:,:))),MR_params.Bval(cpt_b));
                            [~, MRResults.EigValue(cpt_g,cpt_b,:), MRResults.EigVector(cpt_g,cpt_b,:,:),MRResults.MD(cpt_g,cpt_b),MRResults.FA(cpt_g,cpt_b),~]= DiffSim_ToolBox.Tensor(abs(squeeze(MRResults.Att_approx(cpt_g,cpt_b,:,:))),MR_params.Dir,MR_params.Bval(cpt_b));     
                            
                            MRResults.Magn(cpt_g,cpt_b,:)= MRResults.Att_approx(cpt_g,cpt_b,1,:);
                            MRResults.Ph(cpt_g,cpt_b,:)=   mean( MRResults.Phi_dist(cpt_g,cpt_b,1,:,:),5);
                        end
                    end
                  
                   
                     
                    %%% 3) Save the results
                 
             end
        end
        
        %%%% Main functions %%%%% 
		function Simulate_cell_structure(Cell_struct,varargin)
		    
            
            for cpt_cell=1:1:length(Cell_struct.struct)
                tic
                disp([' Simulating cell structure ' num2str(cpt_cell)]);
                
                %%% 0) Init
                Cell_params=Cell_struct.struct{cpt_cell};
                Cells=[];
                Cells_mask=[];
                Mask_id=[];
               
                %%% 1) Simulate a basic 2D cell structure
                [Cells] = Cellule_ToolBox2.Cell_structure(Cell_params.Cell,Cell_params.Voxel,Cell_params.ECV,1000,5,Cell_params.Bridge);
                
                %%% 2) We generate a mask from structure
                [Cells_mask, Mask_id]=Cellule_ToolBox2.Cells_2_Mask(Cells,Cell_params.Voxel,Cell_params.MaskRes);
                Cells_mask(Cells_mask>=1)=1;
                
                %%% 3) Check the real ECV from the Mask
                Cell_params.ECV_calculate=1-sum(Cells_mask(:))/(size(Cells_mask,1)*size(Cells_mask,2)*size(Cells_mask,3));
                
                disp([' Cell structure generated in ' num2str(toc) ' with ECV from mask = ' num2str(Cell_params.ECV_calculate)]);
                
                %%% 4)SAve the results
                save([Cell_params.folder_out '_' num2str(Cell_params.ID) '.mat'],'Cell_params','Cells','Cells_mask','Mask_id')
            end
            
        end
        
        function Simulate_diffusion_distribution(Diff_struct,varargin)
		    
            
            for cpt_diff=1:1:length(Diff_struct.struct)
              
                Diff_params=Diff_struct.struct{cpt_diff};
                CellList = dir([Diff_params.folder_in '*.mat']);  %get list of files and folders in any subfolder
                CellList = CellList(~[CellList.isdir]);  %remove folders from list
                
                for cpt_cell=1:1:length(CellList)
                    tic
                    disp([' Simulating diff distribution ' num2str(cpt_diff) ' with cell structure ' CellList(cpt_cell).name]);
                    
                    %%% 0) Init
                    Water=[];
                    In_out_List=[];
                    DiffResults=[];
                    Diff_params.file_in=[CellList(cpt_cell).folder '/' CellList(cpt_cell).name];
                    
                    load(Diff_params.file_in); 
                    Diff_params.Res=Cell_params.Voxel./size(Cells_mask);
                    
                    Water_in =  DiffSim_ToolBox.Water_distribution_init(Cells_mask,Cell_params.Voxel,Diff_params);
                    %%% 1) Simulate the water inside the cells
                   
                    [Water,In_out_List] =  DiffSim_ToolBox.Water_distribution_mask(Water_in,Cells_mask,Cell_params.Voxel,Diff_params);
                   
                    
                    %%% 2) Calculate the diff parameters of the distribution
                                  
                    DiffResults.ADC = DiffSim_ToolBox.ADC_Water_mean(Water,Diff_params.dT);
                    [~, DiffResults.EigValue,DiffResults.EigVector,DiffResults.MD,DiffResults.FA,~]   = DiffSim_ToolBox.Tensor_Water(Water,Diff_params.dT);

                    DiffResults.ADC_water_in = DiffSim_ToolBox.ADC_Water_mean(Water(find(In_out_List(:,1)),:,:),Diff_params.dT);
                    [~, DiffResults.EigValue_in,DiffResults.EigVector_in,DiffResults.MD_in,DiffResults.FA_in,~]   = DiffSim_ToolBox.Tensor_Water(Water(find(In_out_List(:,1)),:,:),Diff_params.dT);
                    
                    DiffResults.ADC_water_out = DiffSim_ToolBox.ADC_Water_mean(Water(find(~In_out_List(:,1)),:,:),Diff_params.dT);
                    [~, DiffResults.EigValue_out,DiffResults.EigVector_out,DiffResults.MD_out,DiffResults.FA_out,~]   = DiffSim_ToolBox.Tensor_Water(Water(find(~In_out_List(:,1)),:,:),Diff_params.dT);
                    
                    
                    disp([' Diff distribution generated in ' num2str(toc)]);
                     
                    %%% 3) Save the results
                    save([Diff_params.folder_out '_' num2str(Diff_params.ID) '_Cell' num2str(Cell_params.ID) '.mat'],'Diff_params','Cell_params','Water','In_out_List','DiffResults')
                end
            end
            
        end
        
        function Simulate_displ_distribution(Disp_struct,varargin)
         
            for cpt_disp=1:1:length(Disp_struct.struct)
              
                Disp_params=Disp_struct.struct{cpt_disp};
                DiffList = dir([Disp_params.folder_in '*.mat']);  %get list of files and folders in any subfolder
                DiffList = DiffList(~[DiffList.isdir]);  %remove folders from list 
                 
                load(Disp_params.file_in); 
                
                for cpt_diff=1:1:length(DiffList)
                    tic
                    disp([' Simulating Displ distribution ' num2str(cpt_disp) ' with Diff distribution ' DiffList(cpt_diff).name]);
                    
                    Disp_params.file_in_2=[DiffList(cpt_diff).folder '/' DiffList(cpt_diff).name];
                    load(Disp_params.file_in_2); 
                    
                    
                    % We start by 
                    pointR=(nodes_DTI.Rotation*nodes_DTI.points')';
                    pointR2=(inv(nodes_DTI.Rotation)*pointR')';

                    Xt = [-1 -1 -1; ....
                           1 -1 -1; ....
                           1  1 -1; ....
                          -1  1 -1; ....
                          -1 -1 1; ....
                           1 -1 1; ....
                           1  1 1; ....
                          -1  1 1];
                      
                  
                      
                      % For all the Pixels
                      tmp_points=[];
                    for cpt_points=1:1:size(nodes_DTI.points,1)
                       tmp_points=[tmp_points;(inv(nodes_DTI.Rotation)*[pointR(cpt_points,:)+[Xt.*DTI_param.xres/2]]')'];
                        
                      %waitbar(cpt_points/ size(nodes_DTI.points,1) ,h);
                    end
                    
                    T_us=size(Water,3)*10;
                    NDiffphase=ceil(T_us./(DENSE_param(1).temporal_res*1e3))+2;
                    NCardiacPhase=floor(size(nodes_DENSE(1).points,3)/NDiffphase);
                    for cpt_t=1:1:NCardiacPhase
                    
                        [Simpoints] = ComputeF_disp_KM(nodes_DENSE, tmp_points, DTI_param.CardiacPhase,(cpt_t-1)*NDiffphase+1:1:cpt_t*NDiffphase); 
                        disp([' Displ F  generated in ' num2str(toc)]);




                        for cpt_points=1:1:size(nodes_DTI.points,1)


                           % displ=diff(Simpoints((cpt_points-1)*8+1:(cpt_points)*8,:,:),[],3);
                            displ=Simpoints((cpt_points-1)*8+1:(cpt_points)*8,:,:)-Simpoints((cpt_points-1)*8+1:(cpt_points)*8,:,1);
                            x=[0:size(displ,3)-1]*15*1e3; % ms ->us
                            xq=[0:10:max(x)]; % ms ->us
                            vq=[];

                       %  tic
                            for cpt=1:1:8
                                vq(cpt,1,:)=interp1(x,squeeze(displ(cpt,1,:)),xq);
                                vq(cpt,2,:)=interp1(x,squeeze(displ(cpt,2,:)),xq);
                                vq(cpt,3,:)=interp1(x,squeeze(displ(cpt,3,:)),xq);
                            end




                                Dwater=CubeFastInterpolation((Water(:,:,1)-0.25e-3)*4e3,vq(:,:,1:size(Water,3))*1e-3); % onvert the displacement from mm to m
                           % toc
        %                     Xt2=Xt;
        %                     Xt2=Xt2/2+0.5;
        %                     Xt2=Xt2*0.5*1e-3;
        %                     Disp2=zeros(size(Water));
        %                                     
        %                     for cpt_t=1:1:size(Water,3)
        %                         Fx = scatteredInterpolant(Xt2(:,1),Xt2(:,2),Xt2(:,3), vq(:,1,cpt_t));
        %                         Fy = scatteredInterpolant(Xt2(:,1),Xt2(:,2),Xt2(:,3),vq(:,2,cpt_t));
        %                         Fz = scatteredInterpolant(Xt2(:,1),Xt2(:,2),Xt2(:,3),vq(:,3,cpt_t));
        % 
        %                         Disp2(:,1,cpt_t)=Fx(Water(:,1,cpt_t),Water(:,2,cpt_t),Water(:,3,cpt_t));
        %                         Disp2(:,2,cpt_t)=Fy(Water(:,1,cpt_t),Water(:,2,cpt_t),Water(:,3,cpt_t));
        %                         Disp2(:,3,cpt_t)=Fz(Water(:,1,cpt_t),Water(:,2,cpt_t),Water(:,3,cpt_t));
        %                     end


                             % Water2=Water(:,:,1:end)+cumsum(Dwater,3);
                              Water2=Water(:,:,1:end)+Dwater;
                              Water3=Water2-mean(Water2);
                               DisplResults.ADC(cpt_t,cpt_points) = DiffSim_ToolBox.ADC_Water_mean(Water3,Diff_params.dT);
    %                             [~, DiffResults.EigValue,DiffResults.EigVector,DiffResults.MD,DiffResults.FA,~]   = DiffSim_ToolBox.Tensor_Water(Water2,Diff_params.dT);
    % 
    %                             DiffResults.ADC_water_in = DiffSim_ToolBox.ADC_Water_mean(Water2(find(In_out_List(:,1)),:,:),Diff_params.dT);
    %                             [~, DiffResults.EigValue_in,DiffResults.EigVector_in,DiffResults.MD_in,DiffResults.FA_in,~]   = DiffSim_ToolBox.Tensor_Water(Water2(find(In_out_List(:,1)),:,:),Diff_params.dT);
    % 
    %                             DiffResults.ADC_water_out = DiffSim_ToolBox.ADC_Water_mean(Water2(find(~In_out_List(:,1)),:,:),Diff_params.dT);
    %                                [~, DiffResults.EigValue_out,DiffResults.EigVector_out,DiffResults.MD_out,DiffResults.FA_out,~]   = DiffSim_ToolBox.Tensor_Water(Water2(find(~In_out_List(:,1)),:,:),Diff_params.dT);
                        end
                    end
                end
                
                disp([' Displ distribution generated in ' num2str(toc)]);
                %%% 3) Save the results
                save([Disp_params.folder_out '_' num2str(Disp_params.ID) '_Cell' num2str(Cell_params.ID) '_Diff' num2str(Diff_params.ID) '.mat'],'Disp_params','Diff_params','Cell_params','DisplResults')

            end
        end
        
        function Simulate_displ_MR_distribution(Disp_struct,MR_struct,varargin)
         
            % For each DENSE structure
            for cpt_disp=1:1:1 %length(Disp_struct.struct)
                % Load DENSE/DTI struct
                Disp_params=Disp_struct.struct{cpt_disp};
                load(Disp_params.file_in); 
                
                % Prepare the list of Diffusion simulation
                DiffList = dir([Disp_params.folder_in '*.mat']);  %get list of files and folders in any subfolder
                DiffList = DiffList(~[DiffList.isdir]);  %remove folders from li

                % Prepare the DENSE Data

                % Define the Corners of the voxel 
                 Xt = [-1 -1 -1; ....
                       1 -1 -1; ....
                       1  1 -1; ....
                      -1  1 -1; ....
                      -1 -1 1; ....
                       1 -1 1; ....
                       1  1 1; ....
                      -1  1 1].*DTI_param.xres/2;
                      

                % For all the Pixels, shift the corner in their location
                SimPoints_static=[];
                for cpt_points=1:1:size(nodes_DTI.points,1)
                   SimPoints_static(cpt_points,:,:)=nodes_DTI.points(cpt_points,:)+(nodes_DTI.Rotation*Xt')';      
                end


                % Reformat as [Npoint*Ncorner Ndim]
                tmp_points=reshape(SimPoints_static,size(nodes_DTI.points,1)*8,3);

                % For memory efficiency
                clear SimPoints_static;

                % For each diffusion simulation
                for cpt_diff=1:1:1 %length(DiffList)
                    tic
                    disp([' Simulating Displ distribution ' num2str(cpt_disp) ' with Diff distribution ' DiffList(cpt_diff).name]);  

                    % Load Diff simulation
                    Disp_params.file_in_2=[DiffList(cpt_diff).folder '/' DiffList(cpt_diff).name];
                    load(Disp_params.file_in_2); 
                    
                   
                    % The duration of the diffusion simulation in ms
                    T_ms=size(Water,3)*Diff_params.dT*1e3; % us ->ms
                    % The number of DENSE phase the diffusion is covering
                    NDiffphase=ceil(T_ms./(DENSE_param(1).temporal_res))+2;
                    % The number of cardiac phase we can simulate
                    NCardiacPhase=floor(size(nodes_DENSE(1).points,3)/NDiffphase);

                    
                    DisplResults=[];

                    % Prepare the Diffusion simulated Data

                    % Center the water simulation in the voxel
                    Water=Water-repmat(Cell_params.Voxel./2,size(Water,1),1,size(Water,3));

                    % Rotate the Water based on the voxel
                    % orientation then shift it to the correct
                    % location
                    for cpt_tw=1:1:size(Water,3)
                        Water_rot(:,:,cpt_tw)=(nodes_DTI.Rotation*Water(:,:,cpt_tw)')';%+repmat(nodes_DTI.points(cpt_points,:),size(Water,1),1)*1e-3;
                    end
                    DWater_rot=diff(Water_rot(:,:,:),[],3);
                    Water_init=Water_rot(:,:,1);

                    % For memory efficiency
                    clear Water Water_rot; 

                    for cpt_t=1:1:NCardiacPhase
                        disp(['Simulating Cardiac phase ' num2str(cpt_t) ]);

                        % Apply DENSE displacements to all corners of the
                        % voxels
                        [SimPoints_tmp] = ComputeF_disp_KM(nodes_DENSE, tmp_points, DTI_param.CardiacPhase,(cpt_t-1)*NDiffphase+2:1:cpt_t*NDiffphase+1); 
                       
                        % Reformat as [Npoint Ncorner Ndim  Ntime]
                        SimPoints_dynamic=reshape(SimPoints_tmp,size(nodes_DTI.points,1),8,3,size(SimPoints_tmp,3));
                        
                        % For memory efficiency
                        clear SimPoints_tmp; 

                        disp(['Time to prep one cardiac phase ' num2str(toc)]);

                        %For Each Pixel
                        for cpt_points=1:1:1:size(nodes_DTI.points,1)
                            tic

                            % Take the 8 corners of a given voxel
                            displ=squeeze(SimPoints_dynamic(cpt_points,1:8,1:3,:));

                            % reinterpolate the displacement of the corner to the simulation time
                            t=[0:size(displ,3)-1]*DENSE_param(1).temporal_res*1e3; % ms ->us
                            tq=[0:10:max(t)]; % ms ->us
                            vq=[];
                            for cpt=1:1:8
                                vq(cpt,1,:)=interp1(t,squeeze(displ(cpt,1,:)),tq,'spline');
                                vq(cpt,2,:)=interp1(t,squeeze(displ(cpt,2,:)),tq,'spline');
                                vq(cpt,3,:)=interp1(t,squeeze(displ(cpt,3,:)),tq,'spline');
                            end
                        
                            % Calculate the position of the first time
                            % point of the water molecules inside the cube
                            Water_loc=Water_init(:,:,1)+repmat(nodes_DTI.points(cpt_points,:),size(Water_init,1),1)*1e-3;
                            % Interpolate the corner displacement to the
                            % water molecules
                            Disp_Water=CubeFastInterpolation(Water_loc,vq(:,:,1:(size(DWater_rot,3)+1))*1e-3); % convert the displacement from mm to m
                            % Add the delta water displacement to it
                            Disp_Water(:,:,2:end)=Disp_Water(:,:,2:end)+DWater_rot;
                          
                                                     
                            % Simulate MR and save the data
                            MRResults=DiffSim_ToolBox.local_MR_distribution(MR_struct,Disp_Water);
                            DisplResults.ADC_disp_only(cpt_t,cpt_points) = DiffSim_ToolBox.ADC_Water_mean(Disp_Water,Diff_params.dT);
                            DisplResults.MRPhi(:,:,:,:,:,cpt_t)=MRResults.Phi_dist;
                            DisplResults.MD(cpt_t,cpt_points,:)=squeeze(MRResults.MD);
                            DisplResults.FA(cpt_t,cpt_points,:)=squeeze(MRResults.FA) ;
                            DisplResults.Magn(cpt_t,cpt_points,:,:)=  squeeze(MRResults.Magn);
                            DisplResults.Ph(cpt_t,cpt_points,:,:)=    squeeze(MRResults.Ph);
                            disp(['Time to one full pixel cardiac phase ' num2str(toc)]);
                         end 
                        disp(['Time to one full cardiac phase ' num2str(toc)]);
                    end
                    save([Disp_params.folder_out 'TEST_' num2str(Disp_params.ID) '_Cell' num2str(Cell_params.ID) '_Diff' num2str(Diff_params.ID) '_' num2str(cpt_t) '.mat'],'Disp_params','Diff_params','Cell_params','DisplResults','-v7.3')

                end

            end
        end
        
        function Simulate_MR_distribution(MR_struct,varargin)
		    
            
            for cpt_MR=1:1:length(MR_struct.struct)
              
                MR_params=MR_struct.struct{cpt_MR};
                DiffList = dir([MR_params.folder_in '*.mat']);  %get list of files and folders in any subfolder
                DiffList = DiffList(~[DiffList.isdir]);  %remove folders from list
                
                for cpt_diff=1:1:length(DiffList)
                    tic
                    disp([' Simulating MR distribution ' num2str(cpt_MR) ' with Diff distribution ' DiffList(cpt_diff).name]);
                    
                    %%% 0) Init
                    MRResults=[];
                    MR_params.file_in=[DiffList(cpt_diff).folder '/' DiffList(cpt_diff).name];
                    load(MR_params.file_in);  
                                                          
                    for cpt_g=1:1:MR_params.Waveform            % Loop over waveforms
                        for cpt_b=1:1:length(MR_params.Bval)    % Loop over b-values
                           
                            %%% 1) Generate the MR waveform
                            [G, t]= Waveform_ToolBox.Factory(MR_params.Dn(cpt_g,:),MR_params.Gn(cpt_g,:),MR_params.RampTime,MR_params.dT*1e6,MR_params.Bval);
                            if isfield(MR_params,'Rfs')
                                if ~isempty(MR_params.Rfs)
                                    [G, t]= Waveform_ToolBox.Factory2(MR_params.Dn(cpt_g,:),MR_params.Gn(cpt_g,:),MR_params.Rfs(cpt_g,:),MR_params.RampTime,MR_params.dT*1e6,MR_params.Bval);
                                end              
                            end
                            G=G*10^-3; % Convert from mT/m to T/m

                            %%% 2) Simulate the MR experiement 
                            [MRResults.Att_approx(cpt_g,cpt_b,:,:), MRResults.Phi_dist(cpt_g,cpt_b,:,:,:)]= DiffSim_ToolBox.Simulate_MR(Water,MR_params,G) ; % [Waveform Direction] % Grad in [T/m], dimension : NPoints x Waveform with dT corresponding to World.dT  

                            %%% 3) Calculate the ADC and tensor paramaters from the MR simulation and collect the results      
                             MRResults.ADC(cpt_g,cpt_b)    = DiffSim_ToolBox.ADC_mean(abs(squeeze(MRResults.Att_approx(cpt_g,cpt_b,:,:))),MR_params.Bval(cpt_b));
                            [~, MRResults.EigValue(cpt_g,cpt_b,:), MRResults.EigVector(cpt_g,cpt_b,:,:),MRResults.MD(cpt_g,cpt_b),MRResults.FA(cpt_g,cpt_b),~]= DiffSim_ToolBox.Tensor(abs(squeeze(MRResults.Att_approx(cpt_g,cpt_b,:,:))),MR_params.Dir,MR_params.Bval(cpt_b));     
                            
                        end
                    end
                  
                    disp([' MR distribution generated in ' num2str(toc)]);
                     
                    %%% 3) Save the results
                    save([MR_params.folder_out '_' num2str(MR_params.ID) '_Cell' num2str(Cell_params.ID) '_Diff' num2str(Diff_params.ID) '.mat'],'MR_params','Diff_params','Cell_params','In_out_List','DiffResults','MRResults')
                end
            end
            
        end

        
    end % Method Static
    
end % Class
