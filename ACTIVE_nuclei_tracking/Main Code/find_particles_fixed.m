function [particles pellipses weird_vec] = find_particles_fixed( parent_info, parent_vec, x_vector, y_vector, removeflagged, area_thresh, min_area, fit_height )
% function [particles pellipses weird_vec area_vec children] = find_particles_fixed( parent_info, parent_vec, x_vector, y_vector, removeflagged, area_thresh, min_area, fit_height )
%
% Find unique particles from a set of contour data.
%
% Program description:
%   Get the particles: those contours that have no other contours inside them.
%   Fit them with ellipses at half-height, and find if any share common
%   lower-level contours (siblings).
%   Delete "particles" which are actually mutiple peaks within the same particle.
%   Criteria: the underlying sibling has very low aspect ratio
%
%  1/11/2012
%  M. Brasch, R. Baker, L. Manning
%   
% Conditions and terms of use:
% The software packages provided here are M-files executable in MATLAB, a 
% proprietary numerical computing enviornment developed by MathWorks.
% You are free to use this software for research purposes, but you should 
% not redistribute it without the consent of the authors. In addition, end 
% users are expected to include adequate citations and acknowledgments 
% whenever results or derivatives that are based on the software are presented or published.
%
% Citation to ACTIVE should include the following:
% Baker RM, Brasch ME, Manning ML, and Henderson JH. Automated, 
%        contour-based tracking and analysis of cell behavior over long 
%        timescales in environments of varying complexity and cell density.
%        Journal information to be updated when available.
%
% Citations to work foundational to ACTIVE are suggested to include the following, at a minimum:
%
% Idema T. A new way of tracking motion, shape, and divisions. European 
%        Biophysics Journal. 2013:1-8.
% Crocker JC, Grier DG. Methods of digital video microscopy for colloidal 
%        studies. Journal of Colloid and Interface Science. 1996;179(1):298-310.
% Gao Y, Kilfoil ML. Accurate detection and complete tracking of large 
%        populations of features in three dimensions. Optics Express. 
%        2009;17(6):4685-704.
%
%  INPUTS: 
%  parent_info, parent_vec, x_vector, y_vector, are from 
%    makecontourparentarray
%  removeflagged: if set to 1, the function removes "junk" particles, 
%    particles that have been identified in multi-peak noise cases but are 
%    not indicative of a cell
%  area_thresh: pixel area threshold used to determine if cells are 
%    dividing; values below this threshold can be considered divisions
%  min_area: particles with an area lower than min_area are deleted
%  fit_height: contour level to fit the ellipse to
%
%  OUTPUTS:
%  particles: matrix containing [x, y, longaxis, shortaxis, angle, 
%    toplevel, time, topcontour_id, sibling, weird]
%  pellipses: cell array with the ellipse contours for each particle
%  weird_vec: vector of flagged particles


if nargin < 6
    removeflagged = 1;
end
 
%  Local Parameters 
  nParticles = length(parent_info);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%$$Start change
  particles = zeros(nParticles, 10); %xc, yc, a, b, angle, split_pos, contour # at split pos, flagged, multi-body interaction number
  %%$$End change
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  pellipses = cell(nParticles,1);
  flagged = zeros(nParticles, 1);
  weird_vec = zeros(nParticles,3);  % Added 4/20/12 - List of flagged
  % particles with corresponding area
  area_vec = zeros(nParticles,1);    % List of area values for 'non-weird' 
  % particles (num_extra = 0; see loop below)
  children = zeros(nParticles,1);
  
%identify particles and fit ellipses
for i = 1:nParticles
    pivec = parent_info{i};
    %Cells will be fit at half height for tracking. 
    hhindex = ceil(length(pivec)*fit_height);
    hhvalue = pivec(hhindex);
    
    %Divisions and collisions will be related to shared contours. Shared
    %contours can be found by looking at the number of children each
    %contour has (num_children) and comparing it to the number of children
    %each contour should have (1) if you are looking at a single cell only.
    num_children = parent_vec(pivec,2);
    num_extra = nnz(num_children > 1);
    %Number of 'extra' children (anywhere where >1 child exists)
    children(i,1) = num_extra;

    % if num_extra equals zero (nothing weird)
    if num_extra == 0 
        %fit contour with an ellipse
        [xc yc xspread yspread angle es] = fit_ellipse(x_vector{pivec(hhindex)},y_vector{pivec(hhindex)});
        particles(i,1:5) = [xc, yc, xspread, yspread, angle];
        pellipses{i} = es;
        area_vec(i) = pi*xspread*yspread;   % Added 4/20/12 - Transfers area value of 
        % non- weird particles to area vector. 

    %For cases where multiple children exist, first find lowest shared
    %contour (dictates sibling match up for all cases). Then find highest
    %shared contour and determine best fit. 
    else
        %Find lowest split position - should encompass all related children
        split_pos = find(num_children > 1, 1, 'last'); 
        particles(i,6) = split_pos;
        split_contour = pivec(split_pos);        
        particles(i,7) = split_contour;
        flagged(i) = 1;

        %Now fit at top shared split position
        split_pos = find(num_children >1, 1, 'first');
        
        %fix contour fit if cell division is above half max (otherwise
        %fit at half height)
        if split_pos <= hhindex                    
            %fit ellipse at contour directly above split postion
            [xc yc xspread yspread angle es] = fit_ellipse(x_vector{pivec(split_pos-1)},y_vector{pivec(split_pos-1)});
        else
            [xc yc xspread yspread angle es] = fit_ellipse(x_vector{pivec(hhindex)},y_vector{pivec(hhindex)});
        end
        
        %Now use area to determine if junk particles exist:
        area = pi*xspread*yspread;
       
        if (area > min_area) && (area < area_thresh)
            flagged(i) = 1;
        else
            flagged(i) = 2;
        end
    end
    
    %Transfer particle information to overall particle list
    particles(i,1:5) = [xc, yc, xspread, yspread, angle];
    particles(i,8) = flagged(i); 

    if parent_vec(hhvalue, 3) == 1 
        flagged(i) = 3;
        particles(i,8) = 3;
        pellipses{i,:} = es;
        fprintf('Ring particle found; particle %d \n',i);
    else  
        pellipses{i} = es;
    end
end
  
  %remove unwanted particles (multipeak cases)
 if removeflagged ==1 
      
    %first, figure out which contours the junk particles share  
    numbering = 1:nParticles;  
    boolarea = (particles(:,3).*particles(:,4))*2*pi()<min_area; %find particles with lower area than min area
    booljunk = particles(:,8) ==2;  
    junkids = numbering(booljunk);
    junknbrs = cell(length(junkids),1);
    cid = length(parent_vec)*ones(length(junkids),1); %list of the lowest shared contour ids for junk particles
    ictr =0;

    for i=junkids
        ictr = ictr +1;
        jctr = 0;
        for j=junkids
            jctr = jctr +1;
            if i ~=j
                %find out if these two particles share a parent
                [a b c] = intersect(parent_info{i}, parent_info{j});
                if ~isempty(a)
                    junknbrs{ictr} = [junknbrs{ictr}, j];
                    junknbrs{jctr} = [junknbrs{jctr}, i];
                    % keep the largest shared contour between these two
                    cmax = max(a);
                    % keep the smallest of all contours that this particle
                    % shares with other junk particles
                    cid(ictr) = min(cmax, cid(ictr));
                    cid(jctr) = min(cmax, cid(jctr));
                end
            end
        end
    end

    % id1 is an index into junkids that lists the contours we've
    %decided to keep
    [contours id1 id2] = unique(cid);
    ngood = nParticles - length(junkids) + length(contours) - sum(boolarea);

    [nP ncol] = size(particles);
    %make new arrays for particles with junk removed
    tempparticles = zeros(ngood,ncol);
    tempellipses = cell(ngood,1);

    %mapid is a map from the particle index to the parent_contour index
    mapid = zeros(ngood,1);
    igood = 0;

    for i=1:nParticles
        if ((booljunk(i) == 0) || (~isempty(intersect(junkids(id1),i)))) && boolarea(i) == 0
            igood = igood+1;
            tempparticles(igood,:) = particles(i,:);
            tempellipses{igood} = pellipses{i};
            mapid(igood) = i;
        else        
                disp(['Deleting particle ', num2str(i)]);

        end
    end
    
  particles = tempparticles;
  pellipses = tempellipses;
  else
      mapid = zeros(nParticles,1);
 end
 
%Identify dividing particles
cell_divisions = find(particles(:,8) == 1); 
num_cell_div = length(cell_divisions);        

%Check for cell divisions
if num_cell_div>1

    %Find any multi-particle collisions - flag separately
    div_info = particles(cell_divisions,:);
    Sibset=unique(div_info(:,7));
    count=[hist(div_info(:,7),Sibset)]';
    if length(Sibset) == 1
        bool_count = count~=0;
        count = count(bool_count);
    end
    Sibset(:,2) = count;
    multisib = find(Sibset(:,2)>2);
    numweird = length(multisib);
    for i = 1:numweird
        wc_index = multisib(i);     %Find weird contour index value
        wc = Sibset(wc_index,1);    %Find weird contour
        bool_wc = particles(:,7) == wc; %Find all weird particles with this split contour
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%$$ Start change
        particles(bool_wc,10) = i;  
        %%$$ End change
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        particles(bool_wc,8) = 4;   %New flag (4) indicates multi-body cases
    end

    %Update cell division information (remove multiparticle collision info)
    cell_divisions = find(particles(:,8) == 1); 
    num_cell_div = length(cell_divisions);

    %Match up sibling cells
    for i = 1:num_cell_div
        for j = 1:num_cell_div
            if i ~= j && particles(cell_divisions(i),7) == particles(cell_divisions(j), 7)
            sibling = cell_divisions(j);
            particles(cell_divisions(i), 9) = sibling;
            end
        end
    end

end

end
