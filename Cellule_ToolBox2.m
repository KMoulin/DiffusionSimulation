classdef Cellule_ToolBox2

    methods(Static)
                
        %% Function to generate a new cell
        function C=New_Cell(cell_dim,voxel_dim,angle_orien,z_pos,ID)
                C=[];
                
                C.tree=ID.tree;
                C.ID=ID.current;
                C.IDprevious=[];
                
                a = cell_dim(1);
                b = cell_dim(2);
                C.Rad = abs((b-a).*rand(1,1) + a); % Radius

                a = cell_dim(3);
                b = cell_dim(4);
                C.Length = abs((b-a).*rand(1,1) + a); % Length

                a = 0;
                b = voxel_dim(1);
                C.Center(1) = (b-a).*rand(1,1) + a;      % Pos X

                a = 0;
                b = voxel_dim(2);

                C.Center(2) = (b-a).*rand(1,1) + a;      % Pos Y
                
                C.Center(3) = z_pos;                         % Pos Z

                
                C.Angle=angle_orien;
                
                C.Vect=[0 0 1]+ [angle_orien*(rand(2,1)'-0.5) 0];
                

                 C.CV = C.Center + C.Vect;
                 C.BB = C.Length/2; % The radius of the cirular bounding box 
                    
                    
                C.Sur=C.Rad*C.Rad*pi;  %Surface

                C.Vol=C.Sur*C.Length;     %Volume
                
                C.Poly=0;
                
                C.Vertex=[];
        end
        
        
         function C=New_Cell_Linked(cell_dim,voxel_dim,cell,ID)
                C=[];

                
                C.ID=ID.current;
                C.tree=cell.tree;
                C.IDprevious=cell.ID;
                
                a = cell_dim(1);
                b = cell_dim(2);
                C.Rad = abs((b-a).*rand(1,1) + a); % Radius

                a = cell_dim(3);
                b = cell_dim(4);
                C.Length = abs((b-a).*rand(1,1) + a); % Length

               
                C.Center = cell.Center+cell.Vect.*cell.Length;      % Pos X
              
%                 C.Angle=cell.Angle;
%                 C.Vect=[0 0 1]+ [cell.Angle*(rand(2,1)'-0.5) 0];
%                 C.Vect=C.Vect./norm(C.Vect);


                C.Angle=cell.Angle;
                C.Vect=cell.Vect;
                C.Vect=C.Vect./norm(C.Vect);


                C.CV = C.Center + C.Vect;
                C.BB = C.Length/2; % The radius of the cirular bounding box 
                    
                    
                C.Sur=C.Rad*C.Rad*pi;  %Surface

                C.Vol=C.Sur*C.Length;     %Volume
                
                C.Poly=0;
                
                C.Vertex=[];
         end
            function C=New_Cell_Bridged(cell_dim,voxel_dim,cell,ID)
                C=[];

                
                C.ID=ID.current;
                C.tree=cell.tree;
                C.IDprevious=cell.ID;
                
                a = cell_dim(1);
                b = cell_dim(2);
                C.Rad = abs((b-a).*rand(1,1) + a); % Radius

                a = cell_dim(3);
                b = cell_dim(4);
                C.Length = abs((b-a).*rand(1,1) + a); % Length

               
                C.Center = cell.Center;      % Pos X
              
%                 C.Angle=cell.Angle;
%                 C.Vect=[0 0 1]+ [cell.Angle*(rand(2,1)'-0.5) 0];
%                 C.Vect=C.Vect./norm(C.Vect);


                C.Angle=cell.Angle;
                C.Vect=cell.Vect;
                C.Vect=C.Vect./norm(C.Vect);


                C.CV = C.Center + C.Vect;
                C.BB = C.Length/2; % The radius of the cirular bounding box 
                    
                    
                C.Sur=C.Rad*C.Rad*pi;  %Surface

                C.Vol=C.Sur*C.Length;     %Volume
                
                C.Poly=0;
                
                C.Vertex=[];
        end
        %% Function to generate a new polynome cell with random boundary
        function C=New_Cell_Poly(cell_dim,voxel_dim,nb_poly,rng_poly)
             C=[];
            if (rng_poly(1)==rng_poly(2)) || (nb_poly<3)
                C=Cellule_ToolBox2.New_Cell(cell_dim,voxel_dim);
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
        function [stack_cell2, ID]=Add_Cell_Z(stack_cell,cell_dim,voxel_dim,ID)
            %  
            stack_cell2=stack_cell;
            for cpt=1:1:size(stack_cell,1)
                tmpz=stack_cell(cpt).Center(3)+stack_cell(cpt).Vect(3).*stack_cell(cpt).Length;
                tmp_cell=stack_cell(cpt);
               while tmpz<=voxel_dim(3)
                    % Create a new cell witht he same dimension and on top of the previous one
                    ID.current=ID.current+1;
                    cell =Cellule_ToolBox2.New_Cell_Linked(cell_dim,voxel_dim,tmp_cell, ID);
                    tmpz=tmpz+cell.Vect(3).*cell.Length;
                    
                    stack_cell2=[stack_cell2; cell];
                    tmp_cell=cell;
               end
            end
        end
        
        %% Function to add cell in the Z dimension as a function of the voxel dimension
         function [stack_cell2, ID]=Add_Cell_Bridge(stack_cell,cell_dim,voxel_dim,ID)
            %  
            stack_cell2=stack_cell;
                  
            Center=cell2mat({stack_cell(:).Center}');
            Tree=cell2mat({stack_cell(:).tree}');
            
            for cpt=1:1:size(stack_cell,1)
                
                tmp_cell=stack_cell(cpt);
                tmp_dist=sqrt( (Center(cpt,1)-Center(:,1)).^2+ (Center(cpt,2)-Center(:,2)).^2 + (Center(cpt,3)-Center(:,3)).^2);
                tmp_distZ=(Center(:,3)-Center(cpt,3));
                IdxT=find(Tree~=tmp_cell.tree & tmp_distZ>0);
                if ~isempty(IdxT)
                    
                    Idx=find(tmp_dist==min(tmp_dist(IdxT))& tmp_distZ>0);
                    Idx=Idx(1); % In case we have more than one match

                    % Create a new cell with a angled orientation doing the
                    % bridge between two trees

                    ID.current=ID.current+1;
                    tmp_cell.Vect=stack_cell(Idx).Center-tmp_cell.Center;
                    cell_dim(3)=min(tmp_dist(IdxT));
                    cell_dim(4)=cell_dim(3);

                    cell=Cellule_ToolBox2.New_Cell_Bridged(cell_dim,voxel_dim,tmp_cell, ID);

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
                            
       
        
        
        %% Function that inflammate cells until a given ECV                 
         function stack_cell=Expend_cell_ECV(stack_cell,voxel_dim,ECV_cutoff)
                   ECV_tmp=Cellule_ToolBox2.Calculate_ECV(stack_cell,voxel_dim);
                   k=1;
                   
                   h = waitbar(0,['ECV loading ']);
                   while ECV_tmp>=ECV_cutoff && k<10
                        stack_cell=Cellule_ToolBox2.Inflammation(stack_cell,0.1);
                        ECV_tmp=Cellule_ToolBox2.Calculate_ECV(stack_cell,voxel_dim);
                        k=k+1;   
                        
                        waitbar(ECV_tmp/ECV_cutoff,h,['ECV [' num2str(ECV_tmp) '/' num2str(ECV_cutoff) ']']);
                   end
                    
                    close(h);
         end             
         
         
         %% Function that depress cells until a given ECV                  
         function stack_cell=Depress_cell_ECV(stack_cell,voxel_dim,ECV_cutoff)
                   ECV_tmp=Cellule_ToolBox2.Calculate_ECV(stack_cell,voxel_dim);
                   k=1;
                   
                   h = waitbar(0,['ECV loading ']);
                   while ECV_tmp<=ECV_cutoff && k<100
                        stack_cell=Cellule_ToolBox2.Shrink(stack_cell,0.01);
                        ECV_tmp=Cellule_ToolBox2.Calculate_ECV(stack_cell,voxel_dim);
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
                            if(~Cellule_ToolBox2.Collision_Cell_Cells(C,stack_cell,cpt))
                                stack_cell(cpt,:)=C;
                            end
                        end
                    end
         end

        %% Function that randomly shrink cells by a percentage 
        function stack_cell=Shrink_rng(stack_cell,rng_poly)
                    for cpt=1:1:size(stack_cell,1)     
                            stack_cell(cpt,9:9+stack_cell(cpt,8)-1)=stack_cell(cpt,9:9+stack_cell(cpt,8)-1,:) - [(rng_poly(1) + (rng_poly(2)-rng_poly(1)) .* rand(1,stack_cell(cpt,8))).*stack_cell(cpt,9:9+stack_cell(cpt,8)-1,:)];
                            stack_cell(cpt,:)=Cellule_ToolBox2.Calculate_Volume(stack_cell(cpt,:));
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
                            stack_cell(cpt,:)=Cellule_ToolBox2.Calculate_Volume(stack_cell(cpt,:));
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
         function [mask_cell, mask_cell_idx]=Cells_2_Mask(stack_cell,voxel_dim,resolution)
                    
             NbPixel=round(voxel_dim./resolution);
             mask_cell=zeros(NbPixel);
             mask_cell_idx=zeros(NbPixel);
             [X,Y] = ndgrid(1:1:NbPixel(1),1:1:NbPixel(2));
             vectP=[X(:)*resolution(1) Y(:)*resolution(2) zeros(500*500,1)];
             vectZ= (1:1:NbPixel(3))*resolution(3);
              
             
             Center=cell2mat({stack_cell(:).Center}');
             Length=cell2mat({stack_cell(:).Length}');
           
             tol=resolution(3);
             
              for cptZ=1:1:NbPixel(3)
               %tic
                  vectP(:,3)=vectZ(cptZ);
                  
                  Idx= find( (vectZ(cptZ)+tol>=Center(:,3)) & ((Center(:,3)+Length)>=(vectZ(cptZ)-tol)));
                  [collision, ColIdx]=Collision_ToolBox2.Collision_Detection(stack_cell(Idx), vectP); 

                  % allocation the collision to each pixel
                  % it's not efficient to use loop but couldn't get the
                  % indexation to work:
                  % [row,col] = ind2sub([NbPixel(1),NbPixel(2)],1:1:NbPixel(1)*NbPixel(2));
                  cpt=1;
                  for cptX=1:1:NbPixel(1)
                      for cptY=1:1:NbPixel(2)
                          mask_cell(cptX,cptY,cptZ)=collision(cpt);
                          mask_cell_idx(cptX,cptY,cptZ)=ColIdx(cpt);
                          cpt=cpt+1;
                      end

                  end
                %  toc
              end
             
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
            count=Collision_ToolBox2.Collision_Detection(stack_cell, xt);
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
            count=Collision_ToolBox2.Collision_Detection_Mask(mask_cell, xt,Resolution);
            count(count>0)=1;
            ECV=1-mean(count);
        end 
        
         %% Function Collision between one cells and the stack of cells
        function col=Collision_Cell_Cells(cell,stack_cell,idx_expect)
            col=0;
            k=1;
            n_cell=size(stack_cell,1);
            while k<=n_cell 
                   if k~=idx_expect % Pre check if both cells are on the same Z plane and the index is not the same 
                                   if cell.Poly~=0
                                        if Collision_ToolBox2.Poly_Poly(cell,stack_cell(k)) % True Collision detection
                                                col=1;
                                                break;
                                        
                                        end
                                   else
                                       if Collision_ToolBox2.Circle_Circle(cell,stack_cell(k)) % True Collision detection
                                             col=1;
                                             break;
                                       end
                                   end
                                   
                   end
                   k=k+1;
            end
                    
        end
        
        function col=Collision_Cell_Cells2(cell,stack_cell,idx_expect)
            col=0;
            k=1;
            n_cell=size(stack_cell,1);
            
            Center=cell2mat({stack_cell(:).Center}');
            Radius=cell2mat({stack_cell(:).Rad}');
            if idx_expect~=0
                Center=Center(setdiff(1:end,idx_expect),:);
                Radius=Radius(setdiff(1:end,idx_expect));
            end
            
            
             List_col=Collision_ToolBox2.Circle_Circle2(cell.Center,cell.Rad,Center,Radius); % True Collision detection
             col=sum(List_col);                       
                    
        end
        
        function col=Collision_Cell_Cells3(cell,stack_cell,idx_expect)
            
            n_cell=size(stack_cell,1);
            
            Center=cell2mat({stack_cell(:).Center}');
            Radius=cell2mat({stack_cell(:).Rad}');
            Vect=cell2mat({stack_cell(:).Vect}');
            Length=cell2mat({stack_cell(:).Length}');
            if idx_expect~=0
                Center=Center(setdiff(1:end,idx_expect),:);
                Radius=Radius(setdiff(1:end,idx_expect));
                Vect=Vect(setdiff(1:end,idx_expect),:);
                Length=Length(setdiff(1:end,idx_expect));
            end
            Vect1(:,1)=cell.Center;
            Vect1(:,2)=cell.Center+cell.Vect.*cell.Length;
            Vect2(:,:,1)=Center;
            Vect2(:,:,2)=Center+Vect.*Length;
            
            In_out=Collision_ToolBox2.Cylinder_Cylinder(Vect1,Vect2,cell.Rad+Radius);
           % List_col=Collision_ToolBox2.Circle_Circle2(cell.Center,cell.Rad+Center,Radius); % True Collision detection
            col=sum(In_out);                       
                    
        end
        function [C] = Cell_structure(cell_dim,voxel_dim,ECV,max_iter,sup_iter,bridge)   
            %   cell_dim =[minR maxR minZ maxZ];
            %   voxel_dim=[X Y Z];
            %   ICV=[ % ]
            %   mat_iter: define the max number of try before giving up (default 1000)
            %   sup_iter: define the number of time the cell are going to be reorganize
            %   to minimize ECV (default 5)

            disp(['Generate cell Structure ']);
            tic
            ICV=1-ECV;  % [%]

            ICV_sur=ICV*voxel_dim(1)*voxel_dim(2)*voxel_dim(3); % [%] x [mm2] x [mm2]
            tmp_IVC=0;

            n_poly=0;
            vect_poly=[0.8 1];

            rmean=(cell_dim(1)+cell_dim(2))/2;
            vmean=rmean*rmean*pi*voxel_dim(3);


            nb_cell_approx=tmp_IVC/vmean;

            nb_cell_mean_x=ICV*voxel_dim(1)/(rmean*2);
            nb_cell_mean_y=ICV*voxel_dim(2)/(rmean*2);

            mean_x=voxel_dim(1)/nb_cell_mean_x;
            mean_y=voxel_dim(2)/nb_cell_mean_y;
            k=1;

            %% Solve in 2D the cell in the box
            ID.current=1;
            ID.tree=1;

            C=[];
            C=Cellule_ToolBox2.New_Cell(cell_dim,voxel_dim,0.0,0,ID);
            %C=new_cell_local_poly(cell_dim,voxel_dim,n_poly,vect_poly);

            tmp_IVC=tmp_IVC+C.Sur*voxel_dim(3);

            %% Add cells until a given ECV is reach
            h = waitbar(0,['ECV loading ']);
            for cpt_iter=1:1:sup_iter
                k=1;
                while k<max_iter && tmp_IVC < ICV_sur  
                    ID.current=ID.current+1;
                    ID.tree=ID.tree+1;
                    tmp_cell=Cellule_ToolBox2.New_Cell(cell_dim,voxel_dim,0,0,ID);
                    p=1;
                    while p<max_iter && Cellule_ToolBox2.Collision_Cell_Cells2(tmp_cell,C,0)
                        tmp_cell=Cellule_ToolBox2.New_Cell(cell_dim,voxel_dim,0,0,ID);
                        p=p+1;
                    end
                    if p<max_iter
                       C=[C ;tmp_cell]; 
                       tmp_IVC=tmp_IVC+tmp_cell.Sur*voxel_dim(3);
                    end

                    k=k+1;
                     waitbar(tmp_IVC/ICV_sur,h,['ECV [' num2str( 1-tmp_IVC/(voxel_dim(1)*voxel_dim(2)*voxel_dim(3))) '/' num2str(ECV) ']']);
                end

                Center=cell2mat({C(:).Center}');
                % If we don't have reach our ICV/ECV target yet, we need to shake
                % the box to try to gain some space
                if  tmp_IVC < ICV_sur && cpt_iter~=sup_iter

                    % Shake the box toward -X
                    for cpt=1:1:size(C,1)
                        [tt Idx]=sort(Center(:,1));

                        for cpt_x=cpt_iter:1:(10+cpt_iter)
                            tmp_x=mean_x/(cpt_x);
                            C(Idx(cpt)).Center(1)=C(Idx(cpt)).Center(1)-tmp_x;
                            
                            if  C(Idx(cpt)).Center(1)>0
                                if Cellule_ToolBox2.Collision_Cell_Cells2(C(Idx(cpt)),C,Idx(cpt))
                                    C(Idx(cpt)).Center(1)=C(Idx(cpt)).Center(1)+tmp_x;
                                end
                            else
                                 C(Idx(cpt)).Center(1)=C(Idx(cpt)).Center(1)+tmp_x;
                            end
                        end
                    end

                    % Shake the box toward -Y
                    for cpt=1:1:size(C,1)   
                       [tt Idx]=sort(Center(:,2));

                       for cpt_y=cpt_iter:1:(10+cpt_iter)
                            tmp_y=mean_y/(cpt_y);
                            C(Idx(cpt)).Center(2)=C(Idx(cpt)).Center(2)-tmp_y;
                            if  C(Idx(cpt)).Center(2)>0
                                if Cellule_ToolBox2.Collision_Cell_Cells2(C(Idx(cpt)),C,Idx(cpt))
                                    C(Idx(cpt)).Center(2)=C(Idx(cpt)).Center(2)+tmp_y;
                                end
                            else
                                    C(Idx(cpt)).Center(2)=C(Idx(cpt)).Center(2)+tmp_y;
                            end

                       end
                    end

                    % Shake the box toward -X and -Y
                    for cpt=1:1:size(C,1)    
                       [tt Idx]=sort(Center(:,2));

                       for cpt_y=cpt_iter:1:(10+cpt_iter)
                            tmp_x=mean_x/(cpt_x);
                            tmp_y=mean_y/(cpt_y);
                            C(Idx(cpt)).Center(1)=C(Idx(cpt)).Center(1)-tmp_x/2;
                            C(Idx(cpt)).Center(2)=C(Idx(cpt)).Center(2)-tmp_y/2;
                            if  C(Idx(cpt)).Center(2)>0
                                if Cellule_ToolBox2.Collision_Cell_Cells2(C(Idx(cpt)),C,Idx(cpt))
                                   C(Idx(cpt)).Center(1)=C(Idx(cpt)).Center(1)+tmp_x/2;
                                   C(Idx(cpt)).Center(2)=C(Idx(cpt)).Center(2)+tmp_y/2;
                                end
                            else
                                   C(Idx(cpt)).Center(1)=C(Idx(cpt)).Center(1)+tmp_x/2;
                                   C(Idx(cpt)).Center(2)=C(Idx(cpt)).Center(2)+tmp_y/2;
                            end

                       end
                    end 

               end     
            end
            close(h);

            %     figure
            %     DisplayStruct_ToolBox.Cells(C,voxel_dim);
            %     axis([0 (voxel_dim(1)) 0 (voxel_dim(2))])
            %     
            %     if(size(C,2)>7)
            %         figure
            %         DisplayStruct_ToolBox.Cells_Poly(C,voxel_dim);
            %         axis([0 (voxel_dim(1)) 0 (voxel_dim(2))])
            %     end
            %     
            %     %% Add Inflamation
            %     C2=C;
            %     for cpt_step=1:1:5
            %          C2=add_cell_inflammation_local(C2,20,0.1);
            %         if(size(C2,2)>7)
            %             display_cell_poly_local(C2,voxel_dim);
            %             stack_img=add_gif('test_cell_struct',false,stack_img);
            %         end
            %     end
            %     
            %     %% Add apoptosis
            %     %C3=C;
            %     for cpt_step=1:1:3
            %          C2=add_cell_apsotosis_local(C2,20,[0 0.2]);
            %         if(size(C2,2)>7)
            %             display_cell_poly_local(C2,voxel_dim);
            %            stack_img=add_gif('test_cell_struct',false,stack_img);
            %         end
            %     end
            %     
            %     %% Necrosis
            %     
            %     for cpt_step=1:1:3
            %         C2=add_cell_necrosis_local(C2,0.6+0.07*(cpt_step-1));
            %         if(size(C2,2)>7)
            %                 display_cell_poly_local(C2,voxel_dim);
            %                 stack_img=add_gif('test_cell_struct',false,stack_img);
            %         end
            %     end
            %     stack_img=add_gif('test_cell_struct',true,stack_img);
            %% Solve in Z

            
            [C, ID]=Cellule_ToolBox2.Add_Cell_Z(C,cell_dim,voxel_dim,ID);
            if bridge
                [C, ID]=Cellule_ToolBox2.Add_Cell_Bridge(C,cell_dim,voxel_dim,ID);
            end
            display(['Simulated ECV ' num2str( 1-tmp_IVC/(voxel_dim(1)*voxel_dim(2)*voxel_dim(3)))]);

            toc

        end
        
        
        
    end % Method Static
    
end % Class
