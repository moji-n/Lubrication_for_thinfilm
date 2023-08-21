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

function [grid_neigh, x]=grid_gen_refin(del_x_init,x_dim,rc,r_end,x,b,ex_f,if_refinement)
                         

   
    if if_refinement==1

    
        
        dxm=0.05;
        xm=rc-dxm;
        a=(b*ex_f-del_x_init)/(xm^2-rc^2);
        a_=2*a*(xm-rc);
        b=-2*a*rc;
        c=del_x_init+a*xm^2;
        ym=a*xm^2+b*xm+c;
        i=1;
        while x(i)<= xm
            del_x(i)= a_*x(i)+del_x_init;
            i=i+1;
        end
        
        while x(i)<= rc
            del_x(i)=a*x(i)^2+b*x(i)+c;
            i=i+1;
        end


        del_xend=del_x(i-1);

        
        while x(i)<= r_end
            del_x(i)=del_xend;
            i=i+1;
        end

        while x(i)<= r_end+(dxm)
            del_x(i)=a*(x(i)-(r_end-rc))^2+b*(x(i)-(r_end-rc))+c;
            i=i+1;
        end

        while i<= length(x) && x(i)<= x_dim
            del_x(i)= -2*a*(xm-rc)*((x(i)-(r_end-rc))-(rc+(rc-xm)))+ym;
            i=i+1;
        end



        x_sel(1)=x(1);

        delx(1)=del_x(1);
        len=del_x(1);
        j=1;
        while len<x_dim %&& delx(j)<=del_x_init
            
            %indic_x=find(abs(x-len)<1e-6);
        
            vq = spline(x,del_x,len) ;
            %vq=a.*len.^2+b_.*len+del_x_init;
            x_sel(j+1)=len;

            if vq>del_x_init
                break
            else
                j=j+1;
                
                delx(j)=vq;
                len=len+vq;
            end
        end
        

        while len+del_x_init<x_dim
            j=j+1;
            x_sel(j)=len+del_x_init;
            delx(j)=del_x_init;
            len=len+del_x_init;
        end

        x=x_sel;
%         if x(end)<x_dim
%             x(end+1)=x_dim;
% 
%         end


        num_p=length(x);
%          delx=[delx delx(end)];
%          figure()
%          plot(x,gradient(delx,x),'.')
%         hold on 
%         plot(x,gradient(gradient(delx,x),x),'.')
%         hold on 
%         plot(x,gradient(gradient(gradient(delx,x),x),x),'.')
%         hold on 
%         plot(x,gradient(gradient(gradient(gradient(delx,x),x),x),x),'.')
        
        
    else
    

        %num=ceil(x_dim/del_x_init);
        num=del_x_init;
        del_x_init=x_dim/num;
        
        for i=1:num
           x(i)=(i-1)*del_x_init;
        end

        if(x(end)<x_dim)
            x=[x x_dim];
            num=num+1;
        end

        num_p=num;
    end

    grid_neigh=zeros(num_p,4);
        for i=2:num_p-1
            grid_neigh(i,1)=i+1;
            grid_neigh(i,2)=i-1;
        end
        i=1;
        grid_neigh(i,1)=i+1;
        grid_neigh(i,2)=-1;

        
        i=num_p;
        grid_neigh(i,1)=-2;
        grid_neigh(i,2)=i-1;


end

