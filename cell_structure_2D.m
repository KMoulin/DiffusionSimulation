
function [C] = cell_structure_2D(cell_dim,voxel_dim,ECV,max_iter,sup_iter)   
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
    C=Cellule_ToolBox.New_Cell(cell_dim,voxel_dim);
    %C=new_cell_local_poly(cell_dim,voxel_dim,n_poly,vect_poly);
    
    tmp_IVC=tmp_IVC+C(6)*voxel_dim(3);
    
    %% Add cells until a given ECV is reach
    h = waitbar(0,['ECV loading ']);
    for cpt_iter=1:1:sup_iter
        k=1;
        while k<max_iter && tmp_IVC < ICV_sur        
            tmp_cell=Cellule_ToolBox.New_Cell_Poly(cell_dim,voxel_dim,n_poly,vect_poly);
            p=1;
            while p<max_iter && Cellule_ToolBox.Collision_Cell_Cells(tmp_cell,C,0)
                tmp_cell=Cellule_ToolBox.New_Cell_Poly(cell_dim,voxel_dim,n_poly,vect_poly);
                p=p+1;
            end
            if p<max_iter
               C=[C ;tmp_cell]; 
               tmp_IVC=tmp_IVC+tmp_cell(6)*voxel_dim(3);
            end

            k=k+1;
             waitbar(tmp_IVC/ICV_sur,h,['ECV [' num2str( 1-tmp_IVC/(voxel_dim(1)*voxel_dim(2)*voxel_dim(3))) '/' num2str(ECV) ']']);
        end
        
        
        % If we don't have reach our ICV/ECV target yet, we need to shake
        % the box to try to gain some space
        if  tmp_IVC < ICV_sur && cpt_iter~=sup_iter
            
            % Shake the box toward -X
            for cpt=1:1:size(C,1)
                [tt Idx]=sort(C(:,3));
                
                for cpt_x=cpt_iter:1:(10+cpt_iter)
                    tmp_x=mean_x/(cpt_x);
                    C(Idx(cpt),3)=C(Idx(cpt),3)-tmp_x;
                    if C(Idx(cpt),3)>0
                        if Cellule_ToolBox.Collision_Cell_Cells(C(Idx(cpt),:),C,Idx(cpt))
                            C(Idx(cpt),3)=C(Idx(cpt),3)+tmp_x;
                        end
                    else
                         C(Idx(cpt),3)=C(Idx(cpt),3)+tmp_x;
                    end
                end
            end
            
            % Shake the box toward -Y
            for cpt=1:1:size(C,1)   
               [tt Idx]=sort(C(:,4));
                
               for cpt_y=cpt_iter:1:(10+cpt_iter)
                    tmp_y=mean_y/(cpt_y);
                    C(Idx(cpt),4)=C(Idx(cpt),4)-tmp_y;
                    if C(Idx(cpt),4)>0
                        if Cellule_ToolBox.Collision_Cell_Cells(C(Idx(cpt),:),C,Idx(cpt))
                            C(Idx(cpt),4)=C(Idx(cpt),4)+tmp_y;
                        end
                    else
                             C(Idx(cpt),4)=C(Idx(cpt),4)+tmp_y;
                    end
                    
               end
            end
             
            % Shake the box toward -X and -Y
            for cpt=1:1:size(C,1)    
               [tt Idx]=sort(C(:,4));
               
               for cpt_y=cpt_iter:1:(10+cpt_iter)
                    tmp_x=mean_x/(cpt_x);
                    tmp_y=mean_y/(cpt_y);
                    C(Idx(cpt),3)=C(Idx(cpt),3)-tmp_x/2;
                    C(Idx(cpt),4)=C(Idx(cpt),4)-tmp_y/2;
                    if C(Idx(cpt),4)>0
                        if Cellule_ToolBox.Collision_Cell_Cells(C(Idx(cpt),:),C,Idx(cpt))
                            C(Idx(cpt),3)=C(Idx(cpt),3)+tmp_x/2;
                            C(Idx(cpt),4)=C(Idx(cpt),4)+tmp_y/2;
                        end
                    else
                             C(Idx(cpt),3)=C(Idx(cpt),3)+tmp_x/2;
                             C(Idx(cpt),4)=C(Idx(cpt),4)+tmp_y/2;
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
    
    
    
    display(['Simuated ECV ' num2str( 1-tmp_IVC/(voxel_dim(1)*voxel_dim(2)*voxel_dim(3)))]);
       
    toc
        
end




