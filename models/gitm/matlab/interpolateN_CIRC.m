function [y_new, xodd, xodd2] = interpolateN_CIRC(y_old,x_old,x_new)
% Linear N-Dimensional interpolation for data that is defined over
% rectangular grid (might work on weirder grids as well, but it definitely 
% doesnot work with data at random locations). 
%
% CIRC variation assumes first coordinate is circular (last element is connected to
% first, so interpolation on the edges of the first coordinate makes use of
% this assumption.
%
% y_new = interpolateN_CIRC(y_old,x_old,x_new)
%
% y_new is an array of inter-(extra-)polated values 
% y_old is an data defined over a grid defined by x_old{k}
% x_old is a structure with coordinates of the grid (1st coord. is circular!)
% x_new is a structure with coordinates of locations to be
% inter-(extra-)polated to (1st coordinate is circular!)
%
% REQUIREMENTS:
% - y_old must have N dimensions, NONE of which can have length 1 (please
%   use squeeze if you need to)
% -- DO NOT create 3D y_old via MESHGRID (or if you do, do 
%    permute(y_old, [2 1 3]) before giving it to this routine).
% - both, x_old and x_new, must be a structure of Nx1 size.
% -- Each element in x_old and x_new structure must be a vector.
% --- Length of vectors comprising x_old must equal to the corresponding
%     sizes(dimensions) of y_old and they must be increasing arrays
%     (coordinate-arrays corresponding to edges (axis) of y_old = grid-locations
%     of where datapoints in y_old are taken
% --- The vectors comprising x_new must have equal length (what locations
%     do you want to interpolate y_old to?)
%
% NOTE: -it DOES extrapolate (and tells you how many points it had to
%        extrapolate - run or see the last 5 lines of this file.
%
% Algorithm details (optional to read): One interpolation issue is that (1,1,...,1) column is already in
% A (see code), so need to check that A is full rank
% So to do that just take the whole hyper-cube verticies (all 2^N points).
%
% USAGE:
% Example: N=2, y_old is 2x2 
%
% interpolateN( [1 2; 2 3], {[1 3];[1 3]}, {2;2} ) 
%
% Example: N=3, y_old is 2x3x4
% y_old(:,:,1)=[1 2 3;  
%               2 3 4];
% y_old(:,:,2)=[2 3 4; 
%               3 4 5];
% y_old(:,:,3)=[3 4 5; 
%               4 5 6];
% y_old(:,:,4)=[4 5 6; 
%               5 6 7];
% x_old{1,1}=[1 2]; %row-dimension coordinates
% x_old{2,1}=[1 2 3]; %column-dimension coordinates
% x_old{3,1}=[1 2 3 4]; %depth-dimension coordinates
% x_new{1,1}=[1.5 1.999];
% x_new{2,1}=[1   2.999];
% x_new{3,1}=[1   3.5];

% DART $Id$
% CREDIT: Alexey Morozov

N=length(size(y_old)); %number of dims
ni=length(x_new{1}); %number of locations to be interpolated to

y_new=nan(1,ni);
fl=0;
fl2=0;

for k=1:N 
    if size(y_old,k)~=length(x_old{k}) %check if dimensions are right
        error('interpolateN:dimensions','y_old must have the same dimensions in the same order as x_old. Are you using meshgrid? Yes->read ''help interpolateN''')
    end
end

for i=1:ni
   
    a=nan(1,N+1); % the location to be interpolated to
    xi=nan(1,N); %indecies of the cube where the current interpolation location resides wrt the (1,1,1,...) location in y_old
    for k=1:N
        temp=find(x_old{k}(:)>x_new{k}(i));
%         disp([isempty(temp) length(temp) length(x_old{k}) ])
        if isempty(temp) %if we are on the or to the right of the rightmost domain member
            if k==1 %if we are on the circular coordinate
                xi(1,k)=1;
                fl2=fl2+1;
                xodd2(fl2)=i;
            else
                xi(1,k)=length(x_old{k});
                fl=fl+1;
                xodd(fl)=i;
            end
        elseif length(temp)==length(x_old{k}) %if we are on the or to the left of the leftmost domain member
            if k==1 %if we are on the circular coordinate
                xi(1,k)=1;
                fl2=fl2+1;
                xodd2(fl2)=i;
            else
                xi(1,k)=2;
                fl=fl+1;
                xodd(fl)=i;
            end
        else
            xi(1,k)=temp(1);
        end
        a(1,k)=x_new{k}(i);
    end
    a(1,N+1)=1;
    

        h=dec2bin(0:(2^N-1) )-48; %cube vertix coordinates in 0s and 1s
        in=ones(2^N,1)*xi-h; %index of grid, cube vertix coordinates
        % of the cube containing the location to be interpolated to
        
        in(in==0)=length(x_old{1}); %circular modification - wrap the cube only along the first coordinate
        
        
        A=nan(2^N,N+1); % the matrix to be pseudo-inverted, as in Ax=b
        b=nan(2^N,1); % b in Ax=b
        
        for j=1:2^N;
            ind='y_old(';
            for k=1:N
                A(j,k)=x_old{k}(in(j,k));
                ind=[ind 'in(' num2str(j) ',' num2str(k) '),'];
            end
            ind=[ind(1:(end-1)) ')']; %remove the last comma and close the paranthesis
%             j
%             ['[' ind(7:end-1) ']']
%             size(y_old)
%             eval(['[' ind(7:end-1) ']'])
%             eval(ind)
            b(j,1)=eval(ind);
        end
        A(:,N+1)=1;
        
        x = A \ b; %the coefficients of the linear object (line, or plane, or volume, or hyper...) can be found via x=inv(A)*b
        y_new(i)= a*x; %evaluate the linear object at the desired location

end

if fl>0
    disp(['Extrapolation was used at ' num2str(fl) ' locations, out of ' num2str(ni) 'x' num2str(N) ' total, which is ' num2str(100*fl/(ni*N),'%.1f') ' percent.' ])
    disp(['In particular, the troublesome x_new indecies are ' num2str(xodd)]) 
    disp('(the above entries can be duplicated multiple times if the troublesome x_new location extends past domain on multiple dimensions)')
end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
