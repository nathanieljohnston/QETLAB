%%  UNIQUE_PERMS    Computes the unique permutations of a vector
%   This function has one required argument:
%     V: a vector
%
%   PERM_LIST = unique_perms(V) is list of unique permutations of a vector.
%   
%   This function does the same thing as unique(perms(v),'rows'), but is 
%   much faster and less memory-intensive in most cases.
%
%   URL: http://www.qetlab.com/unique_perms

%   requires: nothing
%   author: John D'Errico
%	package: QETLAB 
%	last updated: November 27, 2014

function perm_list = unique_perms( v )

uv = unique(v);
n = length(v);
nu = length(uv);

if nu <= 1
    perm_list = v;
elseif n == nu
    perm_list = perms(v);
else
    perm_list = cell(nu,1);
    for j = 1:nu
        vt = v;
        vt(find(vt==uv(j),1)) = [];
        t = unique_perms(vt);
        perm_list{j} = [repmat(uv(j),size(t,1),1),t];
    end
    perm_list = cell2mat(perm_list);
end

end