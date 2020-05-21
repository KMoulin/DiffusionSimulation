classdef Cellule_ToolBox

    methods(Static)
                
        %% Function to generate a new cell
        function C=New_Cell(cell_dim,voxel_dim)
                    C=[];

                    a = cell_dim(1);
                    b = cell_dim(2);
                    C(1) = abs((b-a).*rand(1,1) + a); % Radius

                    a = cell_dim(3);
                    b = cell_dim(4);
                    C(2) = abs((b-a).*rand(1,1) + a); % Length
                    
                    a = 0;
                    b = voxel_dim(1);
                    C(3) = (b-a).*rand(1,1) + a;      % Pos X
                    
                    a = 0;
                    b = voxel_dim(2);
                    C(4) = (b-a).*rand(1,1) + a;      % Pos Y
                    
                    C(5) = 0;                         % Pos Z
                    
                    C(6)=C(1)*C(1)*pi;  %Surface
                    
                    C(7)=C(5)*C(2);     %Volume
        end
        
        %% Function to generate a new polynome cell with random boundary
        function C=New_Cell_Poly(cell_dim,voxel_dim,nb_poly,rng_poly)
             C=[];
            if (rng_poly(1)==rng_poly(2)) || (nb_poly<3)
                C=Cellule_ToolBox.New_Cell(cell_dim,voxel_dim);
            else
                a = cell_dim(1);
                b = cell_dim(2);
                C(1) = abs((b-a).*rand(1,1) + a); % Radius

                a = cell_dim(3);
                b = cell_dim(4);
                C(2) = abs((b-a).*rand(1,1) + a); % Length

                a = 0;
                b = voxel_dim(1);
                C(3) = (b-a).*rand(1,1) + a;      % Pos X

                a = 0;
                b = voxel_dim(2);
                C(4) = (b-a).*rand(1,1) + a;      % Pos Y


                C(5) = 0;                         % Pos Z


                C(6)=C(1)*C(1)*pi;  %Surface

                C(7)=C(5)*C(2);     %Volume

                C(8)=nb_poly;       % Number of polygon vertices

                C(9:(9+C(8)-1))=[(rng_poly(1) + (rng_poly(2)-rng_poly(1)) .* rand(1,C(8)))*C(1)];   % range of radius of each vertices
            end
        end

        %% Function to add cell in the Z dimension as a function of the voxel dimension
        function stack_cell2=Add_Cell_Z(stack_cell,cell_dim,voxel_dim)
            %  
            stack_cell2=stack_cell;
            for cpt=1:1:size(stack_cell,1)
                tmpz=stack_cell(cpt,5)+stack_cell(cpt,2);
               while tmpz<=voxel_dim(3)
                    % Create a new cell witht he same dimension and on top of the previous one
                    cell=stack_cell(cpt,:);

                    a = cell_dim(3);
                    b = cell_dim(4);
                    cell(2) = abs((b-a).*rand(1,1) + a); % Length

                    cell(5) = tmpz;  % Pos Z
                    cell(7) = cell(6)*cell(2);
                    tmpz=tmpz+cell(2);

                    stack_cell2=[stack_cell2; cell];         
               end
            end
        end
        
        %% Function to generate a new polynome cell with random boundar
        function stack_cell=Add_Cell_Poly(stack_cell,nb_poly)
            if size(stack_cell,2)<8 
                for cpt=1:1:size(stack_cell,1) 
                    stack_cell(cpt,8)=nb_poly;
                    stack_cell(cpt,9:9+stack_cell(cpt,8)-1)=stack_cell(cpt,1); 
                end
            end
        end
                            
        %% Function Collision between one cells and the stack of cells
        function col=Collision_Cell_Cells(cell,stack_cell,idx_expect)
            col=0;
            k=1;
            n_cell=size(stack_cell,1);
            while k<=n_cell 
                   if Collision_ToolBox.Z_Plane(cell,stack_cell(k,:))&&(k~=idx_expect) % Pre check if both cells are on the same Z plane and the index is not the same 
                                   if size(cell,2)>7
                                        if Collision_ToolBox.Poly_Poly(cell,stack_cell(k,:)) % True Collision detection
                                                col=1;
                                                break;
                                        
                                        end
                                   else
                                       if Collision_ToolBox.Circle_Circle(cell,stack_cell(k,:)) % True Collision detection
                                             col=1;
                                             break;
                                       end
                                   end
                                   
                   end
                   k=k+1;
            end
                    
        end
        
        
        %% Function that inflammate cells until a given ECV                 
         function stack_cell=Expend_cell_ECV(stack_cell,voxel_dim,ECV_cutoff)
                   ECV_tmp=Cellule_ToolBox.Calculate_ECV(stack_cell,voxel_dim);
                   k=1;
                   
                   h = waitbar(0,['ECV loading ']);
                   while ECV_tmp>=ECV_cutoff && k<10
                        stack_cell=Cellule_ToolBox.Inflammation(stack_cell,0.1);
                        ECV_tmp=Cellule_ToolBox.Calculate_ECV(stack_cell,voxel_dim);
                        k=k+1;   
                        
                        waitbar(ECV_tmp/ECV_cutoff,h,['ECV [' num2str(ECV_tmp) '/' num2str(ECV_cutoff) ']']);
                   end
                    
                    close(h);
         end             
         
         
         %% Function that depress cells until a given ECV                  
         function stack_cell=Depress_cell_ECV(stack_cell,voxel_dim,ECV_cutoff)
                   ECV_tmp=Cellule_ToolBox.Calculate_ECV(stack_cell,voxel_dim);
                   k=1;
                   
                   h = waitbar(0,['ECV loading ']);
                   while ECV_tmp<=ECV_cutoff && k<10
                        stack_cell=Cellule_ToolBox.Shrink(stack_cell,0.1);
                        ECV_tmp=Cellule_ToolBox.Calculate_ECV(stack_cell,voxel_dim);
                        k=k+1;   
                        
                        waitbar(ECV_tmp/ECV_cutoff,h,['ECV [' num2str(ECV_tmp) '/' num2str(ECV_cutoff) ']']);
                   end
                   
                   close(h);
         end
         
        
         %% Function that inflammate the cell                  
         function stack_cell=Inflammation(stack_cell,infl_val)
                    for cpt=1:1:size(stack_cell,1) 
                        for cpt_poly=1:1:stack_cell(cpt,8)
                            C=stack_cell(cpt,:);
                            C(8+cpt_poly)=C(8+cpt_poly)+C(1)*infl_val;
                            if(~Cellule_ToolBox.Collision_Cell_Cells(C,stack_cell,cpt))
                                stack_cell(cpt,:)=C;
                            end
                        end
                    end
         end

        %% Function that randomly shrink cells by a percentage 
        function stack_cell=Shrink_rng(stack_cell,rng_poly)
                    for cpt=1:1:size(stack_cell,1)     
                            stack_cell(cpt,9:9+stack_cell(cpt,8)-1)=stack_cell(cpt,9:9+stack_cell(cpt,8)-1,:) - [(rng_poly(1) + (rng_poly(2)-rng_poly(1)) .* rand(1,stack_cell(cpt,8))).*stack_cell(cpt,9:9+stack_cell(cpt,8)-1,:)];
                            stack_cell(cpt,:)=Cellule_ToolBox.Calculate_Volume(stack_cell(cpt,:));
                    end
        end

        %% Function that remove cell bellow a given vol
        function stack_cell2=Necrosis(stack_cell,cut_off)
                    stack_cell2=[];
                    for cpt=1:1:size(stack_cell,1)     
                            if stack_cell(cpt,6)> stack_cell(cpt,1)*stack_cell(cpt,1)*pi*cut_off
                                stack_cell2=[stack_cell2; stack_cell(cpt,:)];
                            end         
                    end
        end

        %% Function that shrink cells by a percentage
        function stack_cell=Shrink(stack_cell,prc_poly)
                    for cpt=1:1:size(stack_cell,1)     
                            stack_cell(cpt,9:9+stack_cell(cpt,8)-1)=stack_cell(cpt,9:9+stack_cell(cpt,8)-1,:) - prc_poly.*stack_cell(cpt,9:9+stack_cell(cpt,8)-1,:);
                            stack_cell(cpt,:)=Cellule_ToolBox.Calculate_Volume(stack_cell(cpt,:));
                    end
        end
        
        
        
        %% Calculate the volume for a polynome
        function Poly=Calculate_Volume(Poly)
            area=0;
            % Poly [Radius Length Pos_X Pos_Y Pos_Z %Surface %Volume Nb_Poly p1 p2 p3 p4 ..]
             for cpt_cell=1:1:Poly(8)

                   p1(:,1)=Poly(3)+Poly(8+cpt_cell)*(cos((2*pi*(cpt_cell-1))/(Poly(8))));
                   p1(:,2)=Poly(4)+Poly(8+cpt_cell)*(sin((2*pi*(cpt_cell-1))/(Poly(8))));

                   if cpt_cell==Poly(8)
                        p2(1)=Poly(3)+Poly(8+1)*cos(0);
                        p2(2)=Poly(4)+Poly(8+1)*sin(0);
                   else
                        p2(1)=Poly(3)+Poly(8+cpt_cell+1)*(cos((2*pi*(cpt_cell))/(Poly(8))));
                        p2(2)=Poly(4)+Poly(8+cpt_cell+1)*(sin((2*pi*(cpt_cell))/(Poly(8))));
                   end

                   area = area + (p2(1) + p1(1)) * (p2(2) - p1(2)); 
             end
            Poly (6)= abs(area / 2.0);
            Poly (7)= Poly(6) * Poly(2);          
        end
        
        
        %% Function create a mask from a Stack of Cells                 
         function mask_cell=Cells_2_Mask(stack_cell,voxel_dim,resolution)
             tic       
             NbPixel=round(voxel_dim./resolution);
             mask_cell=zeros(NbPixel);

             [X,Y] = ndgrid(1:1:NbPixel(1),1:1:NbPixel(2));
              vectP=[X(:)*resolution(1) Y(:)*resolution(2) zeros(500*500,1)];
              collision=Collision_ToolBox.Collision_Detection(stack_cell, vectP); 
             
              % allocation the collision to each pixel
              % it's not efficient to use loop but couldn't get the
              % indexation to work:
              % [row,col] = ind2sub([NbPixel(1),NbPixel(2)],1:1:NbPixel(1)*NbPixel(2));
              cpt=1;
              for cptX=1:1:NbPixel(1)
                  for cptY=1:1:NbPixel(2)
                      mask_cell(cptX,cptY,:)=collision(cpt);
                      cpt=cpt+1;
                  end
                  
              end
              
             toc
         end
         
        %% Calculate the ECV for a the current voxel
        function ECV=Calculate_ECV(stack_cell,voxel_dim)
            %% Approximation by randomly generated point inside the volume
            nb_point=10000; % precise enought
            count=[];
            xt=[];
            for cpt_dim=1:1:3
               xt(:,cpt_dim)= rand(nb_point,1)*voxel_dim(cpt_dim);  % Uniform distribution along the voxel
            end
            count=Collision_ToolBox.Collision_Detection(stack_cell, xt);
            count(count>0)=1;
            ECV=1-mean(count);
        end
        
        %% Calculate the ECV for a the current voxel
        function ECV=Calculate_ECV_Mask(mask_cell,voxel_dim)
            %% Approximation by randomly generated point inside the volume
            
            Resolution(1)=voxel_dim(1)./size(mask_cell,1);
            Resolution(2)=voxel_dim(2)./size(mask_cell,2);
            Resolution(3)=voxel_dim(3)./size(mask_cell,3);
            
            nb_point=10000; % precise enought
            count=[];
            xt=[];
            for cpt_dim=1:1:3
               xt(:,cpt_dim)= rand(nb_point,1)*voxel_dim(cpt_dim);  % Uniform distribution along the voxel
            end
            count=Collision_ToolBox.Collision_Detection_Mask(mask_cell, xt,Resolution);
            count(count>0)=1;
            ECV=1-mean(count);
        end
        
                          
    end
end