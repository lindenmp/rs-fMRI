function [Vox_coords, MNI_coords] = GetROICoords(filename)
	% This script will spit out all the coordinates of a mask in MNI mm coordinates. 
	% 
	% The only thing to watch for is whether you want neurological versus radiological orientation. 
	% Just delete the line with x = -x to toggle between the two orientations. 
	% Left brain is negative if you keep x = -x. Left brain is positive if you delete x = -x. 

	orientation = 1; % 1 for left being negative; 0 for left being positive.

	voxSize = 2;

	[hdr,data] = read(filename); %read in mask filename
	[x,y,z] = ind2sub(size(data),find(data)); %voxel coordinates
	Vox_coords = [x,y,z];

	MNIorigin = [45,63,36]; % voxel coords of origin of MNI space i.e. x=0mm is voxel coordinate 46
		
	x = voxSize * (x - MNIorigin(1)); %multiply by voxdim and subtract origin
	y = voxSize * (y - MNIorigin(2));
	z = voxSize * (z - MNIorigin(3));

	if orientation == 1
	    x = -x; %neurological versus radiological
	end

	%x, y and z are now MNI coordinated expressed in mm 
	MNI_coords = [x,y,z];
end