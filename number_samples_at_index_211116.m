function [rep_idx,num_el_in_column]=number_samples_at_index_211116(binvec,histcount)
% feed in a matrix, calculates number of non zero elements in each column
% and matches them to a time, returns repmat of times with non zero
% elements in each column

rep_idx=[];
for i=1:length(binvec(1,:))
    num_el_in_column(i)=sum(histcount(:,i)~=0);
    rep_idx=[rep_idx,repmat(binvec(i),1,num_el_in_column(i))];
end

end