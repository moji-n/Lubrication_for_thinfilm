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

function output = interpolate_Mobility(h1,h2,n_mobility,b,s,phi_der1,mobility_slop,c,h_grad1,h_grad2,slip_length,born_length,h_c,x_c)

    output=zeros(1,6);

   slop_mobiliy1=0;
    if mobility_slop==1 

          slop_mobiliy1=(slip_length^4/born_length^2)*(born_length^2/((h1+h2)*h_c*0.5+born_length)^3-born_length^8/((h1+h2)*h_c*0.5+born_length)^9);
    end


    
    if s(1)==1
        if n_mobility==1
            if h1==h2
                output(1) = h1^5/(h1*b^(10/3)+h1^4);
            else
                output(1) = (h2-h1)/(log(h2/h1)-(b^(10/3)/3.0)*(1.0/(h2^3)-1.0/(h1^3)));
                
            end

        elseif n_mobility==2
            if h1==h2
                output(1) = h1^6/(h1^2*b^(10/3)+h1^4);
            else
                output(1) = h1*h2/(1-(b^(10/3)/3)*(1/h1^2+1/(h1*h2)+1/h2^2));
            end
        elseif n_mobility==3
            if h1==h2 && (h1>1e-20&&h2>1e-20)
                output(1) = h1^7/(h1^3*b^(10/3)+h1^4);
            elseif (h1>1e-20&&h2>1e-20)
                output(1) = h1*h2/(0.5*(1/h1+1/h2)-(b^(10/3)/3)*(1/h1^2+1/(h1*h2)+1/h2^2));
            else
                output(1)=0;
            end
        end

        %output(1)=output(1)+(slop_mobiliy1*h1^2+slop_mobiliy2*h2^2)*0.5;
        output(1)=output(1)+(slop_mobiliy1*((h1+h2)*0.5)^2);

    end


     


    if s(2)==1
        %output(2) = (h1^3+h2^3)/2+(slop_mobiliy1*h1^2+slop_mobiliy2*h2^2)*0.5;
        output(2) = ((0.5*h1+0.5*h2)^3)+(slop_mobiliy1*(h1*0.5+h2*0.5)^2);
    end

    if s(3)==1
        
        %output(3)=0.5*(h1^3+slop_mobiliy1*h1^2)*phi_der1(1)+0.5*(h2^3+slop_mobiliy2*h2^2)*phi_der2(1);
        output(3)=((0.5*h1+0.5*h2)^3+slop_mobiliy1*(h1*0.5+h2*0.5)^2)*phi_der1(1);
    elseif s(3)==2

        %output(3)=0.5*(h1^3+slop_mobiliy1*h1^2)*phi_der1(1)+0.5*(h2^3+slop_mobiliy2*h1^2)*phi_der2(1);
        output(3)=((0.5*h1+0.5*h2)^3+slop_mobiliy1*(h1*0.5+h2*0.5)^2)*phi_der1(1);

        if h1==h2 && (h1>1e-20 && h2>1e-20)
            output(6) =c*h1^5/(h1*c*b^(10/3)+h1^4);

        elseif (h1>1e-20 && h2>1e-20)
            output(6) = c*(h2-h1)/(log(h2/h1)-(b^(10/3)/3.0)*c*(1.0/(h2^3)-1.0/(h1^3)));
        else
            output(6) =0;
        end

        %output(6)=output(6)+c*(slop_mobiliy1+slop_mobiliy2)*0.5;
        output(6)=output(6)+c*(slop_mobiliy1);

    end

    if s(4)==1
        
        %output(4)=0.5*(h1^3+slop_mobiliy1*h1^2)*phi_der1(2)+0.5*(h2^3+slop_mobiliy2*h2^2)*phi_der2(2);
        output(4)=((0.5*h1+0.5*h2)^3+slop_mobiliy1*(h1*0.5+h2*0.5)^2)*phi_der1(2);
        
    end


    if s(5)==1
        

        %output(5)=0.5*(h1^3+slop_mobiliy1*h1^2)*phi_der1(3)+0.5*(h2^3+slop_mobiliy2*h2^2)*phi_der1(3);
        output(5)=0.5*(h1^3+slop_mobiliy1*(h1*0.5+h2*0.5)^2)*phi_der1(3);
        
    end

  


    


end
