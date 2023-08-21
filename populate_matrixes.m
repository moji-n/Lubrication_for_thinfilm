% Copyright: Mojtaba Norouzisadeh

%  e-mail: norouzi.mojtaba.sade@gmail.com
%
% date: 08/2023
% 
%
%
%This file is part of Lubrication_for_thinfilm.
%
%    Lubrication_for_thinfilm is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    Lubrication_for_thinfilm is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with Diffmicro.  If not, see <https://www.gnu.org/licenses/>.
%

function [Coef_matrix, const_matrix]=populate_matrixes(num_p,grid_neigh,h_temp,x,del_x,del_mid_x,theta,mob_inter,mob_inter_temp,del_t,g_coef,dis_coef,Coef_matrix, const_matrix,period_bc)
         
 

    for i=1:num_p
        
        Coef_matrix(i,i)=Coef_matrix(i,i)+1;
        const_matrix(i,1)=const_matrix(i,1)+h_temp(i);
        for j=1:length(grid_neigh(i,:))
             

            if grid_neigh(i,j)>0
                i_n=grid_neigh(i,j);

                

                if i>i_n && i<num_p
                    mob_int_i=mob_inter(i_n,:);
                    mob_int_i_temp=mob_inter_temp(i_n,:);
                    mob_int_i_far=mob_inter(i,:);
                    mob_int_i_far_temp=mob_inter_temp(i,:);
                    d_x=del_x(i_n);
                    d_x_far=del_x(i);
                    
                elseif  i==num_p
                    mob_int_i=mob_inter(i_n,:);
                    mob_int_i_temp=mob_inter_temp(i_n,:);
                    if period_bc==0
                        mob_int_i_far=zeros(1,6);%mob_inter(i_n,:);
                        mob_int_i_far_temp=zeros(1,6);%mob_inter_temp(i_n,:);
                    elseif period_bc==1
                        mob_int_i_far=mob_inter(1,:);
                        mob_int_i_far_temp=mob_inter_temp(1,:);
                    end
                    d_x=del_x(i_n);
                    d_x_far=del_x(i_n);
                        
                elseif i<i_n && i>1
                    mob_int_i=mob_inter(i,:);
                    mob_int_i_temp=mob_inter_temp(i,:);
                    mob_int_i_far=mob_inter(i-1,:);
                    mob_int_i_far_temp=mob_inter_temp(i-1,:);
                    d_x=del_x(i);
                    d_x_far=del_x(i-1);

                elseif i==1
                    mob_int_i=mob_inter(i,:);
                    mob_int_i_temp=mob_inter_temp(i,:);
                    if period_bc==0
                        mob_int_i_far=zeros(1,6);%mob_inter(i,:);
                        mob_int_i_far_temp=zeros(1,6);%mob_inter_temp(i,:);
                    elseif period_bc==1
                        mob_int_i_far=mob_inter(num_p-1,:);
                        mob_int_i_far_temp=mob_inter_temp(num_p-1,:);
                    end

                    d_x=del_x(i);
                    d_x_far=del_x(i);

                end
    
                x1=x(i_n);
                %x0=x(i);
    
                
                a_t=-g_coef*(mob_int_i(2)/(del_mid_x(i)*d_x))*del_t ...
                    -dis_coef*(mob_int_i(3)/(del_mid_x(i)*d_x))*del_t ...
                    -dis_coef*(mob_int_i(4)/(del_mid_x(i)*d_x))*del_t ...
                    -dis_coef*(mob_int_i(5)/(del_mid_x(i)*d_x))*del_t;
                
    
                a_temp=-g_coef*(mob_int_i_temp(2)/(del_mid_x(i)*d_x))*del_t ...
                       -dis_coef*(mob_int_i_temp(3)/(del_mid_x(i)*d_x))*del_t ...
                       -dis_coef*(mob_int_i_temp(4)/(del_mid_x(i)*d_x))*del_t ...
                       -dis_coef*(mob_int_i_temp(5)/(del_mid_x(i)*d_x))*del_t;
                
    
    
                Coef_matrix(i,i)=Coef_matrix(i,i)-theta*a_t;
                Coef_matrix(i,i_n)=Coef_matrix(i,i_n)+theta*a_t;
                const_matrix(i,1)= const_matrix(i,1)-(1.0-theta)*(-a_temp*h_temp(i)+a_temp*h_temp(i_n));
    
                for jj=1:length(grid_neigh(i_n,:))
    
                    if grid_neigh(i_n,jj)>0 && grid_neigh(i_n,jj)~= i
                        i_nn=grid_neigh(grid_neigh(i,j),jj);
                        x2=x(i_nn);
                        
                        

    
                        a_2_t=del_t*(mob_int_i(1)-mob_int_i(6))/(del_mid_x(i)*d_x*del_mid_x(i_n)*abs(x1-x2));
                        a_2_temp=del_t*(mob_int_i_temp(1)-mob_int_i_temp(6))/(del_mid_x(i)*d_x*del_mid_x(i_n)*abs(x1-x2));
    
                        a_t=-(del_t/del_mid_x(i))*(     (mob_int_i(1)    -mob_int_i(6))/d_x*(1/(del_mid_x(i_n)*d_x)+(1/(del_mid_x(i)*d_x)))...
                                                    +   (mob_int_i_far(1)-mob_int_i_far(6))/(d_x_far*del_mid_x(i)*d_x));
                        a_temp=-(del_t/del_mid_x(i))*(  (mob_int_i_temp(1)-mob_int_i_temp(6))/d_x*(1/(del_mid_x(i_n)*d_x)+(1/(del_mid_x(i)*d_x)))...
                                                    +   (mob_int_i_far_temp(1)-mob_int_i_far_temp(6))/(d_x_far*del_mid_x(i)*d_x));
    
                        
                        Coef_matrix(i,i_nn)=Coef_matrix(i,i_nn)+theta*a_2_t;
                        Coef_matrix(i,i_n)=Coef_matrix(i,i_n)+theta*(-a_2_t+a_t);
                        
                        Coef_matrix(i,i)=Coef_matrix(i,i)-theta*a_t;
    
                        const_matrix(i,1)= const_matrix(i,1)-(1.0-theta)*(-a_temp*h_temp(i)+(-a_2_temp+a_temp)*h_temp(i_n)+a_2_temp*h_temp(i_nn));
                    elseif grid_neigh(grid_neigh(i,j),jj)<0 && grid_neigh(i_n,jj)~= i
                        %left bc

                       if period_bc==0
    
                            a_t=-(del_t/del_mid_x(i))*(     (mob_int_i(1)-mob_int_i(6))/d_x*(1/(del_mid_x(i_n)*d_x)+(1/(del_mid_x(i)*d_x)))...
                                                          + (mob_int_i_far(1)-mob_int_i_far(6))/(d_x_far*del_mid_x(i)*d_x)    );
                            
                            a_temp=-(del_t/del_mid_x(i))*(  (mob_int_i_temp(1)-mob_int_i_temp(6))/d_x*(1/(del_mid_x(i_n)*d_x)+(1/(del_mid_x(i)*d_x)))...
                                                        +   (mob_int_i_far_temp(1)-mob_int_i_far_temp(6))/(d_x_far*del_mid_x(i)*d_x)        );
        
    
                            Coef_matrix(i,i_n)=Coef_matrix(i,i_n)+theta*(a_t);
                            
                            Coef_matrix(i,i)=Coef_matrix(i,i)-theta*(+a_t);
        
                            const_matrix(i,1)= const_matrix(i,1)-(1.0-theta)*((-a_temp)*h_temp(i)+(a_temp)*h_temp(i_n));
                       elseif period_bc==1
                            a_2_t=del_t*(mob_int_i(1)-mob_int_i(6))/(del_mid_x(i)*d_x*del_mid_x(i_n)*abs(del_x(i)));

                            a_2_temp=del_t*(mob_int_i_temp(1)-mob_int_i_temp(6))/(del_mid_x(i)*d_x*del_mid_x(i_n)*abs(del_x(i)));

                            a_t=-(del_t/del_mid_x(i))*(     (mob_int_i(1)-mob_int_i(6))/d_x*(1/(del_mid_x(i_n)*d_x)+(1/(del_mid_x(i)*d_x)))...
                                                          + (mob_int_i_far(1)-mob_int_i_far(6))/(d_x_far*del_mid_x(i)*d_x)    );
                            
                            a_temp=-(del_t/del_mid_x(i))*(  (mob_int_i_temp(1)-mob_int_i_temp(6))/d_x*(1/(del_mid_x(i_n)*d_x)+(1/(del_mid_x(i)*d_x)))...
                                                        +   (mob_int_i_far_temp(1)-mob_int_i_far_temp(6))/(d_x_far*del_mid_x(i)*d_x)        );
                            
                            if grid_neigh(grid_neigh(i,j),jj)==-1
                                i_nn=num_p-1;
                            elseif grid_neigh(grid_neigh(i,j),jj)==-2
                                i_nn=2;
                            end
                            Coef_matrix(i,i_nn)=Coef_matrix(i,i_nn)+theta*a_2_t;

                            Coef_matrix(i,i_n)=Coef_matrix(i,i_n)++theta*(-a_2_t+a_t);
                            
                            Coef_matrix(i,i)=Coef_matrix(i,i)-theta*(+a_t);
                            
                            const_matrix(i,1)= const_matrix(i,1)-(1.0-theta)*((-a_temp)*h_temp(i)+(-a_2_temp+a_temp)*h_temp(i_n)+a_2_temp*h_temp(i_nn));


                       end

                    end
                 
    
                end
    
    
            elseif grid_neigh(i,j)<0
                if period_bc==0
                    if  i==num_p
                        mob_int_i=mob_inter(i-1,:);
                        mob_int_i_temp=mob_inter_temp(i-1,:);
                        mob_int_i_far=mob_inter(i-1,:);
                        mob_int_i_far_temp=mob_inter_temp(i-1,:);
    
                    elseif i==1
                        mob_int_i=mob_inter(i,:);
                        mob_int_i_temp=mob_inter_temp(i,:);
                        mob_int_i_far=mob_inter(i,:);
                        mob_int_i_far_temp=mob_inter_temp(i,:);

                        


                    end
    
    
                    %left bc
                    for jjj=1:length(grid_neigh(i,:))
                        if grid_neigh(i,jjj)>0 && grid_neigh(i,jjj)~=i
                            i_mirror=grid_neigh(i,jjj);
                            break
                        end
                    end


                elseif period_bc==1
                    if  i==num_p
                        mob_int_i=mob_inter(1,:);
                        mob_int_i_temp=mob_inter_temp(1,:);
                        %mob_int_i_far=mob_inter(i-1,:);
                        %mob_int_i_far_temp=mob_inter_temp(i-1,:);
                        i_mirror=2;
    
                    elseif i==1
                        mob_int_i=mob_inter(num_p-1,:);
                        mob_int_i_temp=mob_inter_temp(num_p-1,:);
                        %mob_int_i_far=mob_inter(i,:);
                        %mob_int_i_far_temp=mob_inter_temp(i,:);
                        i_mirror=num_p-1;
                    end
    
                end
    
    



    


                a_t=0;%-g_coef*(mob_int_i(2)/(del_mid_x(i)*abs(del_x(i))))*del_t ...
                    %-dis_coef*(mob_int_i(3)/(del_mid_x(i)*abs(del_x(i))))*del_t ...
                    %-hyd_coef*(mob_int_i(4)/(del_mid_x(i)*abs(del_x(i))))*del_t ...
                    %-edl_coef*(mob_int_i(5)/(del_mid_x(i)*abs(del_x(i))))*del_t;
                a_temp=0;%-g_coef*(mob_int_i_temp(2)/(del_mid_x(i)*abs(del_x(i))))*del_t...
                       %-dis_coef*(mob_int_i_temp(3)/(del_mid_x(i)*abs(del_x(i))))*del_t ...
                       %-hyd_coef*(mob_int_i_temp(4)/(del_mid_x(i)*abs(del_x(i))))*del_t ...
                       %-edl_coef*(mob_int_i_temp(5)/(del_mid_x(i)*abs(del_x(i))))*del_t;
    
                Coef_matrix(i,i)=Coef_matrix(i,i)-theta*a_t;
                Coef_matrix(i,i_mirror)=Coef_matrix(i,i_mirror)+theta*a_t;
                const_matrix(i,1)= const_matrix(i,1)-(1.0-theta)*(-a_temp*h_temp(i)+a_temp*h_temp(i_mirror));


                if period_bc==0
    
                            
               elseif period_bc==1

                   if grid_neigh(i,j)==-1
                       mob_int_i=mob_inter(num_p-1,:);
                        mob_int_i_temp=mob_inter_temp(num_p-1,:);
                        mob_int_i_far=mob_inter(i,:);
                        mob_int_i_far_temp=mob_inter_temp(i,:);
                        i_n=num_p-1;
                        i_nn=num_p-2;
                    elseif grid_neigh(i,j)==-2
                        mob_int_i=mob_inter(1,:);
                        mob_int_i_temp=mob_inter_temp(1,:);
                        mob_int_i_far=mob_inter(num_p-1,:);
                        mob_int_i_far_temp=mob_inter_temp(num_p-1,:);
                        i_n=2;
                        i_nn=3;
                   end

                    d_x=del_x(i);
                    a_2_t=del_t*(mob_int_i(1)-mob_int_i(6))/(del_mid_x(i)*d_x*del_mid_x(i_n)*abs(del_x(i)));

                    a_2_temp=del_t*(mob_int_i_temp(1)-mob_int_i_temp(6))/(del_mid_x(i)*d_x*del_mid_x(i_n)*abs(del_x(i)));

                    a_t=-(del_t/del_mid_x(i))*(     (mob_int_i(1)    -mob_int_i(6))/d_x*(1/(del_mid_x(i_n)*d_x)+(1/(del_mid_x(i)*d_x)))...
                                                    +   (mob_int_i_far(1)-mob_int_i_far(6))/(d_x_far*del_mid_x(i)*d_x));
                    a_temp=-(del_t/del_mid_x(i))*(  (mob_int_i_temp(1)-mob_int_i_temp(6))/d_x*(1/(del_mid_x(i_n)*d_x)+(1/(del_mid_x(i)*d_x)))...
                                                    +   (mob_int_i_far_temp(1)-mob_int_i_far_temp(6))/(d_x_far*del_mid_x(i)*d_x));
    
                    
                    
                    Coef_matrix(i,i_nn)=Coef_matrix(i,i_nn)+theta*a_2_t;

                    Coef_matrix(i,i_n)=Coef_matrix(i,i_n)++theta*(-a_2_t+a_t);
                    
                    Coef_matrix(i,i)=Coef_matrix(i,i)-theta*(+a_t);
                    
                    const_matrix(i,1)= const_matrix(i,1)-(1.0-theta)*((-a_temp)*h_temp(i)+(-a_2_temp+a_temp)*h_temp(i_n)+a_2_temp*h_temp(i_nn));
                end
        end
    
    end

end
