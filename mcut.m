function D = mcut(C,i)
% To remove ith row and ith column from C (size: NxN) and
% return D (size: N-1xN-1)
% This function is useful in imposing boundary condition onto
% the global stiffness and force vectors.
[m,n] = size(C);
d1 = C(1:i-1,1:i-1);
d2 = C(1:i-1,i+1:n);
d3 = C(i+1:m,1:i-1);
d4 = C(i+1:m,i+1:n);
D = [d1 d2; d3 d4]; 