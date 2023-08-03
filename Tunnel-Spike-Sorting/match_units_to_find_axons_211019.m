function [unit_match_mat] = match_units_to_find_axons_211019(nmi_mat,nmi_thresh)
%This function finds unit pairs which optimizing based maximum NMI
%   The algorithm finds the max NMI value in the matrix, stores the row and
%   col unit numbers, and then looks for the same in the remaining
%   rows,col. 
%   Note: the current version does not care about the duplicate NMI values,
%   it will just choose the first pair.
%
%   Input - 
%   nmi_mat = 2D matrix containing NMI values (double matrix)
%   nmi_thresh = NMI threshold (double)
%   
%   Output - 
%   unit_match_mat = matched unit pairs with 1st col containing row
%   info for nmi_mat, and 2nd col containing col info
%
% Update 10/19/21: Updated to ouput all cluster matches above the threshold

% no_pairs = min(size(nmi_mat));
% unit_match_mat = [];
% for pairi = 1:no_pairs
    % apply nmi threshold 
    threshi = nmi_mat>nmi_thresh;
    nmi_above_thresh = nmi_mat;
    nmi_above_thresh(~threshi) = 0;
    if ~sum(nmi_above_thresh(:))
        unit_match_mat=[];
    end
    
    % find max nmi
    [row, col] = find(nmi_above_thresh~=0);
%     if length(row)>2
%         warning 'Non-uniques NMI values found'
%     end
    
    unit_match_mat = [row, col];
    % substitute rows and cols with -1
%     nmi_mat(row,:)=-1;
%     nmi_mat(:,col)=-1;
%end

