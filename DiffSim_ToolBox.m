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
                for cpt_dir=1:1:size(dwi,2)
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
            dx=abs(Water(:,:,end)-Water(:,:,1));   
            for cpt_dir=1:1:size(dir,1)
               
                % Make the projection on the corresponding axis 
                dx_dir=dx.*(repmat(abs(dir(cpt_dir,:)),size(Water,1),1));
                
                % Take the norm of the vectors
                dx_norm=sqrt(dx_dir(:,1).*dx_dir(:,1)+dx_dir(:,2).*dx_dir(:,2)+dx_dir(:,3).*dx_dir(:,3));
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
   
        KD = sqrt(D * 2 * Dim * dT);
        Dvar=KD * randn(N,1);
         
        VectorDir=randn(N,3);
        VectorDir=VectorDir./repmat(Vector_ToolBox.Norm_vect_n(VectorDir),1,3);
        
        Dstep=VectorDir.*Dvar;
      
        end
        %% Simulate MR from a vector gradient file and water 
        function [Att_approx, Phi_Dist]= Simulate_MR(Water,Seq,Grad,World)
        %%% Simulation gradients
        Phi_Dist=[];
        Nb_Water=size(Water,1);

            %%% Loop over Direction %%%
            for cpt_dir=1:1:size(Seq.Dir,1)
                %%% Loop over Waveform and or B-values %%%
                for cpt_g=1:1:size(Grad,2)          
                     %%% Intregral Diffusion + gradients  %%%           
                     phi=zeros(Nb_Water,1);   
                     for cpt=2:1:size(Grad,1)
                        phi(:)=phi(:)+ nansum( World.Gamma *World.dT*( repmat(Grad(cpt,cpt_g),Nb_Water,3).*repmat(Seq.Dir(cpt_dir,:),Nb_Water,1).*squeeze(Water(:,:,cpt))),2);
                     end

                    %%% Phase Distribution
                    Phi_Dist(cpt_g,cpt_dir,:)=phi(:);
                    %%% MR Signal attenuation
                    Att_approx(cpt_g,cpt_dir)=nanmean(cos(phi(:))+i*sin(phi(:)));  
                    %%%
                end
            end
        end
    
    end 
end
