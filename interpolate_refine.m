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

function output=interpolate_refine(h_t,num_p,num_p_old,x,x_old)

    h_t_new=interp1(x_old(1:num_p_old),h_t(1:num_p_old),x(1:num_p));
    b=h_t_new(end);
    i=0;
    while isnan(b)
        i=i+1;
       b=h_t_new(end-i);
    end
    for j=1:i
        h_t_new(end-j+1)=b;
    end


    for i=1:num_p
        if h_t_new(i)<0
            jb=i-1;
            while h_t_new(jb)<0
                jb=jb-1;
            end

            ja=i+1;
            while h_t_new(ja)<0
                ja=ja+1;
            end
            h_t_new(i)=(x(i)-x(jb))/(x(ja)-x(jb))*h_t_new(jb)+(x(ja)-x(i))/(x(ja)-x(jb))*h_t_new(ja);                
        end

    end
    output=h_t_new;

    for i=1:num_p
        if h_t_new(i)<0
            fprintf("error in refinement interpolation")
            output=-1;
            break
        end
    end
    

    

end
