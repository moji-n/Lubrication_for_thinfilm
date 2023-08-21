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

function output=read_input(case_dir,file)

if strcmp(file,'parameters')
    fileID_parameters = fopen(strcat(case_dir,'/input/parameters.txt'),'r');

    if fileID_parameters<0 
        fprintf("input file %s could not be found \n",strcat(case_dir,'/input/parameters.txt'))
    else

        while ~feof(fileID_parameters)
          st = fgetl(fileID_parameters);           
           if  contains(st,'cluster')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(1)=str2double(cell2mat(stread{1,1}(2)));
           end
    
           if  contains(st,'n_mobility')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(2)=str2double(cell2mat(stread{1,1}(2)));
           end
    
           if  contains(st,'precursor')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(3)=str2double(cell2mat(stread{1,1}(2)));
           end
    
           if  contains(st,'plot')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(4)=str2double(cell2mat(stread{1,1}(2)));
           end
    
           if  contains(st,'laplace_loop_max')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(5)=str2double(cell2mat(stread{1,1}(2)));
           end
    
           if  contains(st,'roughness')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(6)=str2double(cell2mat(stread{1,1}(2)));
           end
    
           if  contains(st,'alpha_rough')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(7)=str2double(cell2mat(stread{1,1}(2)));
           end
    
           if  contains(st,'del_t_i')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(8)=str2double(cell2mat(stread{1,1}(2)));
           end
    
           if  contains(st,'del_t_max')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(9)=str2double(cell2mat(stread{1,1}(2)));
           end
    
           if  contains(st,'end_time')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(10)=str2double(cell2mat(stread{1,1}(2)));
           end
    
           if  contains(st,'gravity')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(11)=str2double(cell2mat(stread{1,1}(2)));
           end
    
           if  contains(st,'surf_tension')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(12)=str2double(cell2mat(stread{1,1}(2)));
           end
    
           if  contains(st,'density')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(13)=str2double(cell2mat(stread{1,1}(2)));
           end
    
           if  contains(st,'x_c')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(14)=str2double(cell2mat(stread{1,1}(2)));
           end
    
           if  contains(st,'h_c')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(15)=str2double(cell2mat(stread{1,1}(2)));
           end
           if  contains(st,'viscosity')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(16)=str2double(cell2mat(stread{1,1}(2)));
           end
    
           if  contains(st,'expanson_fac')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(17)=str2double(cell2mat(stread{1,1}(2)));
           end
    
            if  contains(st,'x_dim')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(18)=str2double(cell2mat(stread{1,1}(2)));
           end
            if  contains(st,'num_threads')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(19)=str2double(cell2mat(stread{1,1}(2)));
            end
    
            if  contains(st,'start_from_zero')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(20)=str2double(cell2mat(stread{1,1}(2)));
            end

            if  contains(st,'delta_x_max')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(21)=str2double(cell2mat(stread{1,1}(2)));
            end


            if  contains(st,'init_profile')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(22)=str2double(cell2mat(stread{1,1}(2)));
            end

            if  contains(st,'crank_nichelson')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(23)=str2double(cell2mat(stread{1,1}(2)));
            end

            if  contains(st,'potential_c')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(24)=str2double(cell2mat(stread{1,1}(2)));
            end

            if  contains(st,'potential_g')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(25)=str2double(cell2mat(stread{1,1}(2)));
            end

            if  contains(st,'potential_vdw')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(26)=str2double(cell2mat(stread{1,1}(2)));
            end

            if  contains(st,'potential_edl')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(27)=str2double(cell2mat(stread{1,1}(2)));
            end

            if  contains(st,'potential_hyd')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(28)=str2double(cell2mat(stread{1,1}(2)));
            end

            if  contains(st,'refinement')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(29)=str2double(cell2mat(stread{1,1}(2)));
            end

            if  contains(st,'period_bc')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(30)=str2double(cell2mat(stread{1,1}(2)));
            end

            if  contains(st,'mobility_slop')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(31)=str2double(cell2mat(stread{1,1}(2)));
            end

            if  contains(st,'slip_length')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(32)=str2double(cell2mat(stread{1,1}(2)));
            end
            
            if  contains(st,'delta_x_min')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(33)=str2double(cell2mat(stread{1,1}(2)));
            end


    
        
    
    
        end
    
        fprintf("successfully read the simulation parameters \n")
        fclose(fileID_parameters);
    end

end


if strcmp(file,'timesteps')
    fileID_timesteps = fopen(strcat(case_dir,'/input/timesteps.txt'),'r');

    if fileID_timesteps<0 
        fprintf("input file %s could not be found \n",strcat(case_dir,'/input/timesteps.txt'))
    else
        if fileID_timesteps>0 
            time_step=fscanf(fileID_timesteps,'%E');
        end
        output=time_step;

        fprintf("successfully read %d number of time steps to write  \n",length(time_step))
        fclose(fileID_timesteps);
    end

end



if strcmp(file,'output')
    fileID_sum = fopen(case_dir,'r');
    
    if fileID_sum<0 
        fprintf("input file %s could not be found \n",case_dir)
    else

        [~,name,~] = fileparts(case_dir);
    
        frs_line= textscan(fileID_sum, '%f %f %f', 1);
        output{1}=frs_line{1};
        output{2}=frs_line{2};
        output{3}=[];
        output{4}=[];
        output{5}=[];
        output{6}=[];
        output{7}=[];
        output{8}=frs_line{3};


        temp_read_output = textscan(fileID_sum, '%f %f %f %f %d %d %d %d');
        
            output{3}=temp_read_output{2};
            output{4}=temp_read_output{3};
            output{7}=temp_read_output{1};
            output{9}={temp_read_output{5},temp_read_output{6},temp_read_output{7},temp_read_output{8} };
                    output{10}=temp_read_output{4};

        output{5}=str2double(name);
        %output{6}=temp_read_output{}
        fprintf("successfully read the date from time step %d \n",str2double(name))

        fclose(fileID_sum)  ;
        
        
    end

end



if strcmp(file,'dis_parameters')

    fileID_dis_para = fopen(strcat(case_dir,'/input/disjoining_parametres.txt'),'r');

    if fileID_dis_para<0 
        fprintf("input file %s could not be found \n",strcat(case_dir,'/input/timesteps.txt'))
    else
        while ~feof(fileID_dis_para)
          st = fgetl(fileID_dis_para);           
           if  contains(st,'ham_cte1')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(1)=str2double(cell2mat(stread{1,1}(2)));
           end

           if  contains(st,'ham_cte2')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(2)=str2double(cell2mat(stread{1,1}(2)));
           end

           if  contains(st,'born_length')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(3)=str2double(cell2mat(stread{1,1}(2)));
           end

           if  contains(st,'k_b')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(4)=str2double(cell2mat(stread{1,1}(2)));
           end

           if  contains(st,'Temp')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(5)=str2double(cell2mat(stread{1,1}(2)));
           end

           if  contains(st,'epsilon0')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(6)=str2double(cell2mat(stread{1,1}(2)));
           end

           if  contains(st,'epsilon')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(7)=str2double(cell2mat(stread{1,1}(2)));
           end

           if  contains(st,'electron_charge')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(8)=str2double(cell2mat(stread{1,1}(2)));
           end

           if  contains(st,'avagadro_number')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(9)=str2double(cell2mat(stread{1,1}(2)));
           end

           if  contains(st,'nb')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(10)=str2double(cell2mat(stread{1,1}(2)));
           end

           if  contains(st,'zeta_1')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(11)=str2double(cell2mat(stread{1,1}(2)));
           end

           if  contains(st,'zeta_2')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(12)=str2double(cell2mat(stread{1,1}(2)));
           end

           if  contains(st,'lambda_0')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(13)=str2double(cell2mat(stread{1,1}(2)));
           end

           if  contains(st,'w_o')
                stread = textscan(st,'%s','Delimiter',{')','(',' '},'MultipleDelimsAsOne',1);
                output(14)=str2double(cell2mat(stread{1,1}(2)));
           end

        end
    
        fprintf("successfully read disjoining pressure parameters \n")
        fclose(fileID_dis_para);
    end


end


end

