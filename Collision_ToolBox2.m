classdef Collision_ToolBox2

    methods(Static)
        
                %% General function to manage Collision between Objects and Points
                 function [In_out_List,In_out_Idx] = Collision_Detection(Objects, Points)
                      In_out_List=false(size(Points,1),1);
                      In_out=false(size(Points,1),1);
                      In_out_Idx=zeros(size(Points,1),1);
                      for cpt_cell=1:1:size(Objects,1)
                          if Objects(1).Poly>0 % Cell is a polygon 
                            In_out= Collision_ToolBox2.Poly_Points_Bounding_Box(Objects(cpt_cell),Points); 
                            In_out_Idx(In_out)=Objects(cpt_cell).ID;
                          else               % Cell is a circle
                            In_out= Collision_ToolBox2.Cylinder_Points(Objects(cpt_cell),Points);  
                            In_out_Idx(In_out)=Objects(cpt_cell).ID;
                          end 
                          In_out_List=In_out_List+In_out; % List of the molecule associate to the corresponding cells 
                          
                      end
                      
                 end
                %% General function to manage Collision between Objects and Points
                 function In_out_List = Collision_Detection_Mask(Mask, Points,Resolution)
                      In_out_List=false(size(Points,1),1);
                      
                      % From absolute coordinate to pixel coordinate  from
                      % 1 to 500
                      for cpt_dim=1:1:3
                          Points(:,cpt_dim)=round(Points(:,cpt_dim)./Resolution(cpt_dim))+1;
                          Idx=find(Points(:,cpt_dim)>size(Mask,cpt_dim));
                          Points(Idx,cpt_dim)=size(Mask,cpt_dim);
                          
                          Idx=find(Points(:,cpt_dim)<1);
                          Points(Idx,cpt_dim)=1;
                      end 
                      Idx=sub2ind(size(Mask),Points(:,1),Points(:,2),Points(:,3));
                      In_out_List=Mask(Idx);
                 end
                %% Function which manage the permeability 
                function [In_out] = Permeability(In_out_before,In_out_after,Perma)
                           List_diff=In_out_before-In_out_after; % 0 nothing ; -1 From cells to extra ; 1 From extra to cells
                           List_perma=find(List_diff~=0);
                           Roll_perma=(rand(length(List_perma),1));
                           In_out=List_perma(find(Roll_perma(Roll_perma>Perma)));
                end
                
                %% Collision between Circle and a list of points
                function   In_out=Circle_Points(Circle,Points)
                      tmp_dist =  sqrt( (Circle.Center(1)-Points(:,1)).^2 + (Circle.Center(2)-Points(:,2)).^2);
                      in_xy = ( tmp_dist<=Circle.Rad );       
                      In_out= in_xy ;
                end
                 
                function   In_out=CylinderZ_Points(Cyl,Points)
                    %% Collision between a Cylinder aligned with the Z axis and a list of points
                      tmp_dist =  sqrt( (Cyl.Center(1)-Points(:,1)).^2 + (Cyl.Center(2)-Points(:,2)).^2);
                      in_xy = ( tmp_dist<=Cyl.Rad );
                      in_z =  ( Points(:,3)>= ( Cyl.Center(3)  ) ) & ( Points(:,3)< ( Cyl.Center(3) + Cyl.Length  ) );       
                      In_out= in_xy & in_z ;
                end
                
                %% Collision between Sphere and a list of points
                function   In_out=Sphere_Points(Sphere,Points)
                     In_out=(sqrt( sum(((Sphere.Center+Sphere.Vect.*Sphere.Length/2)-Points).^2,2))<=Sphere.BB); % Modif KM check here
                end
                
                 
                 function   In_out=Cylinder_Points(Cyl,Points)
                 %% Collision between a Cylinder and a list of points
                                       
                     % Bounding box test first
                       In_out_BB=Collision_ToolBox2.Sphere_Points(Cyl,Points);   
                       In_out_C=false(size(Points,1),1);
                       In_out_C(In_out_BB)= (vecnorm(cross(Points(In_out_BB,:)-Cyl.Center, Points(In_out_BB,:)-Cyl.CV)')<=Cyl.Rad);
                       
                       In_out=In_out_C & In_out_BB;
                       
                 end
                
                 
                 function   In_out=Cylinder_Cylinder(Vect1,Vect2,Rad)
                 %% Collision between a Cylinder and a list of points
                                       
                     % Bounding box test first
                     In_out=false(size(Vect2,1),1);
                     for cpt=1:1:size(Vect2,1)
                       [~, ~, rd] = Collision_ToolBox2.closestDistanceBetweenLines(Vect1(:,1), Vect1(:,2), Vect2(cpt,:,1), Vect2(cpt,:,2), 1, 1, 1, 1, 1);
                       if rd<Rad(cpt)
                          In_out(cpt)=true; 
                       end
                     end
                       %end                   
                 end
                 
                 
               
               
                
                
                  %% Collision between Circle and a list of points
                function   In_out=Z_Plane(Object1,Object2)
                      In_out =  ( Object2(5)>= ( Object1(5)  ) ) & ( Object2(5)< ( Object1(5) + Object1(2) ) );       
                end
                
                %% Collision between Polynome and a list of points
                function   In_out=Poly_Points(Poly,Points)     
                     % Poly [Radius Length Pos_X Pos_Y Pos_Z %Surface %Volume Nb_Poly p1 p2 p3 p4 ..]
                       in_xy =false(1,size(Points,1));
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
                           tmp_col=(((p1(2) >= Points(:,2) & p2(2) < Points(:,2)) | (p1(2) < Points(:,2) & p2(2) >= Points(:,2))) & ( Points(:,1) < (p2(1)-p1(1))*(Points(:,2)-p1(2)) / (p2(2)-p1(2))+p1(1))); 
                           in_xy(tmp_col)=~in_xy(tmp_col);
                       end
                      % in_z =  ( Points(:,3)>= ( Circle(5)  ) ) & ( Points(:,3)< ( Cells(5) + Cells(2) ) );   
                      in_z =  ( Points(:,3)>= ( Poly(5)  ) ) & ( Points(:,3)< ( Poly(5) + Poly(2) ) );       
                      In_out= in_xy' & in_z ;
                end
                 
                 
                %% Collision between Polynome and a list of points
                function   In_out=Poly_Points_Bounding_Box(Poly,Points)     
                     % Poly [Radius Length Pos_X Pos_Y Pos_Z %Surface %Volume Nb_Poly p1 p2 p3 p4 ..]
                        In_out =false(size(Points,1),1); 
                        In_out_z=false(size(Points,1),1);  
                        In_out_ob=false(size(Points,1),1);  
                        In_out_ib=false(size(Points,1),1); 
                        
                        tmp_Poly1=Poly;
                       
                        %%  Z solver 
                        In_out_z =  ( Points(:,3)>= ( Poly(5)  ) ) & ( Points(:,3)< ( Poly(5) + Poly(2) ) ); 
                        
                        
                        % Outerbounding box check
                        tmp_Poly1(1)=max(tmp_Poly1(9:end));
                        In_out_ob(In_out_z,:)=Collision_ToolBox2.Circle_Points(tmp_Poly1,Points(In_out_z,:)); 
                        
                        % InnerBounding box check
                        tmp_Poly1(1)=min(tmp_Poly1(9:end));       
                        In_out_ib(In_out_ob)=Collision_ToolBox2.Circle_Points(tmp_Poly1,Points(In_out_ob,:)); 
                        
                        In_out(In_out_ib)=true; % We are sure that these ones are in;
                        
                        % We test which remains
                        Idx_poly=In_out_ob&~In_out_ib;
                        
                        % Points that are outside the max box will never be in touching, check if it's a true collision
                        In_out(Idx_poly)=Collision_ToolBox2.Poly_Points(Poly,Points(Idx_poly,:));
                end
                 
                
                 %% Collision between two Polynomes
                 function   In_out=Poly_Poly(Poly1,Poly2)            
                   % Poly [Radius Length Pos_X Pos_Y Pos_Z %Surface %Volume Nb_Poly p1 p2 p3 p4 ..]
                   In_out=false;
                   pp1=[];
                   pp3=[];
                   tmp_Poly1=Poly1;
                   tmp_Poly1(1)=max(tmp_Poly1(9:end));
                   tmp_Poly2=Poly2;
                   tmp_Poly2(1)=max(tmp_Poly2(9:end));
                   
                 % if in_z =  ( Points(:,3)>= ( Poly(5)  ) ) & ( Points(:,3)< ( Poly(5) + Poly(2) ) );   
                   if Collision_ToolBox2.Circle_Circle(tmp_Poly1,tmp_Poly2) % Outerbounding box check
                        % OuterBounding box are touching, check if it's a true collision
                       tmp_Poly1(1)=min(tmp_Poly1(9:end));
                       tmp_Poly2(1)=min(tmp_Poly2(9:end));
                       if ~Collision_ToolBox2.Circle_Circle(tmp_Poly1,tmp_Poly2) % InnerBounding box check
                           % InnerBounding box are not touching, check if there is a vertice collision 
                           for cpt_pol=1:1:Poly1(8)
                               p1(:,1)=Poly1(3)+Poly1(8+cpt_pol)*(cos((2*pi*(cpt_pol-1))/(Poly1(8))));
                               p1(:,2)=Poly1(4)+Poly1(8+cpt_pol)*(sin((2*pi*(cpt_pol-1))/(Poly1(8))));
                               p3(:,1)=Poly2(3)+Poly2(8+cpt_pol)*(cos((2*pi*(cpt_pol-1))/(Poly2(8))));
                               p3(:,2)=Poly2(4)+Poly2(8+cpt_pol)*(sin((2*pi*(cpt_pol-1))/(Poly2(8))));
                               if cpt_pol==Poly1(8)
                                    p2(1)=Poly1(3)+Poly1(8+1)*cos(0);
                                    p2(2)=Poly1(4)+Poly1(8+1)*sin(0);
                                    p4(1)=Poly2(3)+Poly2(8+1)*cos(0);
                                    p4(2)=Poly2(4)+Poly2(8+1)*sin(0);
                               else
                                    p2(1)=Poly1(3)+Poly1(8+cpt_pol+1)*(cos((2*pi*(cpt_pol))/(Poly1(8))));
                                    p2(2)=Poly1(4)+Poly1(8+cpt_pol+1)*(sin((2*pi*(cpt_pol))/(Poly1(8))));
                                    p4(1)=Poly2(3)+Poly2(8+cpt_pol+1)*(cos((2*pi*(cpt_pol))/(Poly2(8))));
                                    p4(2)=Poly2(4)+Poly2(8+cpt_pol+1)*(sin((2*pi*(cpt_pol))/(Poly2(8))));
                               end
                               
                               if Collision_ToolBox2.Poly_Line(Poly2,p1,p2)
                                   In_out=true; % Vertice collision
                                   break;
                               end
                           end
                       else
                           In_out=true; % Innerbounding box are Touching, this is a collision collide
                       end
                   else
                       In_out=false; % Outerbounding box not Touching, they will never collide
                   end   
                 end

                %% Collision between Polynome and a Line
                function   In_out=Poly_Line(Poly,p3,p4)
                     % Poly [Radius Length Pos_X Pos_Y Pos_Z %Surface %Volume Nb_Poly p1 p2 p3 p4 ..]
                        In_out=false;
                           for cpt_pol=1:1:Poly(8)
                               p1(:,1)=Poly(3)+Poly(8+cpt_pol)*(cos((2*pi*(cpt_pol-1))/(Poly(8))));
                               p1(:,2)=Poly(4)+Poly(8+cpt_pol)*(sin((2*pi*(cpt_pol-1))/(Poly(8))));
                               if cpt_pol==Poly(8)
                                    p2(1)=Poly(3)+Poly(8+1)*cos(0);
                                    p2(2)=Poly(4)+Poly(8+1)*sin(0);
                               else
                                    p2(1)=Poly(3)+Poly(8+cpt_pol+1)*(cos((2*pi*(cpt_pol))/(Poly(8))));
                                    p2(2)=Poly(4)+Poly(8+cpt_pol+1)*(sin((2*pi*(cpt_pol))/(Poly(8))));
                               end
                               if Collision_ToolBox2.Line_Line(p1,p2,p3,p4)
                                  In_out=true;
                                  break
                               end
                           end
                end

                %% Collision between a Line and a Line
                function   In_out=Line_Line(p1,p2,p3,p4)
                     uA = ((p4(1)-p3(1))*(p1(2)-p3(2)) - (p4(2)-p3(2))*(p1(1)-p3(1))) / ((p4(2)-p3(2))*(p2(1)-p1(1)) - (p4(1)-p3(1))*(p2(2)-p1(2)));
                     uB = ((p2(1)-p1(1))*(p1(2)-p3(2)) - (p2(2)-p1(2))*(p1(1)-p3(1))) / ((p4(2)-p3(2))*(p2(1)-p1(1)) - (p4(1)-p3(1))*(p2(2)-p1(2))); 
                     %   // if uA and uB are between 0-1, lines are colliding
                    if (uA >= 0 && uA <= 1 && uB >= 0 && uB <= 1)    
                          In_out=true;
                    else
                          In_out=false;
                    end
                end

                            % Pos Z
                    
                    
                function   In_out=Circle_Circle(Circle1,Circle2)
                    In_out=false;
                    a = Circle1.Center(1) - Circle2.Center(1);
                    b = Circle1.Center(2) - Circle2.Center(2);
                    dist=sqrt(a*a+b*b);
                    if dist < (Circle1.Rad + Circle2.Rad) % distance too small = collision
                        In_out=true;
                    end
                end
                
                function   In_out=Circle_Circle2(Center1,Radius1,Center2,Radius2)
                    a = Center1(1) - Center2(:,1);
                    b = Center1(2) - Center2(:,2);
                    dist=sqrt(a.*a+b.*b);
                    In_out= (dist < (Radius1 + Radius2)); % distance too small = collision
                        
                end
               
                %% Collision inside a ROI for a list of points
                function   In_out=ROI_Points(ROI,Points)     
                     % Poly [Radius Length Pos_X Pos_Y Pos_Z %Surface %Volume Nb_Poly p1 p2 p3 p4 ..]
                       in_xy =false(1,size(Points,1));
                      for cpt_cell=1:1:size(ROI,1)
                           p1(:,1)=ROI(cpt_cell,1);
                           p1(:,2)=ROI(cpt_cell,2);
                           if cpt_cell<size(ROI,1)
                                p2(:,1)=ROI(cpt_cell+1,1);
                                p2(:,2)=ROI(cpt_cell+1,2);
                           else
                               p2(:,1)=ROI(1,1);
                               p2(:,2)=ROI(1,2);
                           end
                           tmp_col=(((p1(2) >= Points(:,2) & p2(2) < Points(:,2)) | (p1(2) < Points(:,2) & p2(2) >= Points(:,2))) & ( Points(:,1) < (p2(1)-p1(1))*(Points(:,2)-p1(2)) / (p2(2)-p1(2))+p1(1))); 
                           in_xy(tmp_col)=~in_xy(tmp_col);
                       end 
                      In_out= in_xy ;
                end
                 
                function [pA, pB, rd] = closestDistanceBetweenLines(a0, a1, b0, b1, clampAll,  clampA0,  clampA1, clampB0, clampB1)

                    pA=[];
                    pB=[];
                    rd=[];
                    if (clampAll)
                        clampA0 = 1;
                        clampA1 = 1;
                        clampB0 = 1;
                        clampB1 = 1;
                    end

                    A = soma3(a1, a0, -1);
                    B = soma3(b1, b0, -1);
                    A = multiplica3(A, 1 / norma3(A));
                    B = multiplica3(B, 1 / norma3(B));
                    cross = cross3(A,B);
                    denom = norma3(cross).^2;

                    if (denom == 0)
                        d0 = dot3(A, soma3(b0, a0, -1));
                        d = norma3(soma3(soma3(multiplica3(A, d0), a0, 1), b0, -1));
                        if (clampA0 || clampA1 || clampB0 || clampB1)
                             d1 = dot3(A, soma3(b1, a0, -1));
                            if (d0 <= 0 && 0 >= d1)
                                if (clampA0 && clampB1)
                                    if (abs(d0) < abs(d1))
                                        pA = b0;
                                        pB = a0;
                                        rd = norma3(soma3(b0, a0, -1));

                                    else
                                        pA = b1;
                                        pB = a0;
                                        rd = norma3(soma3(b1, a0, -1));
                                    end
                                end
                            elseif (d0 >= norma3(A) && norma3(A) <= d1)
                                if (clampA1 && clampB0)
                                    if (abs(d0) <abs(d1))
                                        pA = b0;
                                        pB = a1;
                                        rd = norma3(soma3(b0, a1, -1));

                                    else
                                        pA = b1;
                                        pB = a1;
                                        rd = norma3(soma3(b1, a1, -1));
                                    end
                                end
                            end
                        else
                            pA = NULL;
                            pB = NULL;
                            rd = d;
                        end
                    else
                        t = soma3(b0, a0, -1);
                        det0 = determinante3(t,B,cross);
                        det1 = determinante3(t,A,cross);
                        t0 = det0 / denom;
                        t1 = det1 / denom;
                        pA = soma3(a0, multiplica3(A, t0),1);
                        pB = soma3(b0, multiplica3(B, t1),1);

                        if (clampA0 || clampA1 || clampB0 || clampB1)
                            if (t0 < 0 && clampA0)
                                pA = a0;
                            elseif (t0 > norma3(A) && clampA1)
                                pA = a1;
                            end
                            if (t1 < 0 && clampB0)
                                pB = b0;
                            elseif (t1 > norma3(B) && clampB1)
                                pB = b1;
                            end
                        end

                        d = norma3(soma3(pA, pB, -1));

                        pA = pA;
                        pB = pB;
                        rd = d;
                    end


                function d=determinante3( a, v1, v2)

                    d= a(1) * (v1(2) * v2(3) - v1(3) * v2(2)) + a(2) * (v1(3) * v2(1) - v1(1) * v2(3)) + a(3) * (v1(1) * v2(2) - v1(2) * v2(1));
                end

                function v=cross3(v1, v2)
                    v=[];
                    v(1) = v1(2) * v2(3) - v1(3) * v2(2);
                    v(2) = v1(3) * v2(1) - v1(1) * v2(3);
                    v(3) = v1(1) * v2(2) - v1(2) * v2(1);

                end

                    function dv=dot3(v1, v2)
                        dv=sum(v1.*v2);
                    end

                function nv=norma3(v1)
                    nv=sqrt(sum(v1.^2));
                end

                function v2=multiplica3(v1, v)
                    for cpt_i=1:3 
                        v2(cpt_i) = v1(cpt_i) * v;

                    end
                end

                function v=soma3(v1, v2, sinal)
                v=[];
                    for cpt_i=1:3 
                        v(cpt_i) = v1(cpt_i) + sinal * v2(cpt_i);
                    end
                end

                end

    end
end