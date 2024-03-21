function [Mobs,flag,r,rate] = LTIobsv(A,C)
%{
This function checks the observability of a given LTI system and provides 
a measure of the degree of observability as exposed in: "Road roughness 
identification in vehicle vertical dynamics: behind the scenes".
Authors: A. Leanza, S. De Carolis, and L. Soria
Department of Mechanics, Mathematics and Management
Polytechnic University of Bari, Italy
%}


% ============= Check observability =============

N = length(A(1,:));
Mobs = obsv(A,C);
r = rank(Mobs);
if r < N
    flag = 'NOT observable';
elseif r == N
    flag = 'observable';
end


% ============= Check the correctness of observability =============

allminors1elemnt = nchoosek(Mobs(:,1),N);
Nminors = length(allminors1elemnt);
ranghi = zeros(Nminors,1);
for i = 1:Nminors
    [row,~] = find(Mobs(:,1)==allminors1elemnt(i,:));
    ranghi(i) = rank(Mobs(row,:));
end
if  ~isempty(find(ranghi==N,1))
    r = N;
    flag = 'observable';
end


% ============= Degree of observability =============

size_Mobs = size(Mobs);
MobsN = zeros(size_Mobs);
for i = 1:size_Mobs(1)
    MobsN(i,:) = Mobs(i,:)/norm(Mobs(i,:));
end
[VV,DD] = eig(MobsN'*MobsN);
DD = diag(DD);
[eig_min,I] = min(DD);
rate.eig_min = eig_min;
rate.VVmin = abs(VV(:,I));
[~,rate.lessObsvState] = max(rate.VVmin);  % Less observable state;
rate.UBobsv = 1./sum(abs(pinv(MobsN)),2);  % Upper-bound observability
rate.VARobsv = 1./sum((pinv(MobsN)).^2,2);  % Variance observability
rate.STDobsv = 1./sqrt(sum((pinv(MobsN)).^2,2));  % Standard-deviation observability

end
