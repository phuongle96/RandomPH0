
%% INITIALIZATION STEP %%
n = 2;
% There are 2n+2 pts, 2n endpts of n intervals and \alpha and \omega
% with 2n+1 small segments
m = 2*n;
% From bottom to top, probs is the probability that start from currently
% considered ith position, it will go down to the (i-1)th positions. 
% Note that the first position start from the second lowest endpt among 
% 2n endpts, and the last point we consider is the last endpt of n intervals 


%TRANSITION PROB%
syms b1 b2 b3 d1 d2 d3
syms P1 P2 P3
pts = [b1, b2,d2,d1];
%Prob=[P1,P2,P3];
probs = sym([1 m-1]);
for i = 1:(m-2)
    probs(i) = (pts(i)/pts(i+1) - pts(i)/pts(i+2))/(1-pts(i)/pts(i+2));
    %probs(i)=Prob(i);
end
probs(m-1) = pts(m-1)/pts(m); 

%probs(m-1) = Prob(m-1); 


%LABEL OF EDGES
% From bottom to top, labels is the label that we used to mark each edge as
% xi or yi (Convention: yi lower than xi) or xiyj or n (no label)
% We also start from second lowest endpt. We label the edge from i
% to (i-1) labels(i)
labels = {'y1', 'x2y2', 'x1'};

%{
Explicit probs and labels for overlapping case:
probs = [(b1/b2-b1/d1)/(1-b1/d1), (b2/d1-b2/d2)/(1-b2/d2), d1/d2];
labels = {'y1', 'x1y2', 'x2'};
%}

% f are two arrays storing coefficient (a multimonial with vars xi,yj) 
% so that F_current is alternated between f1 and f2 with the final result
% being f2
f = sym([2 4^n]);
for t1=1:2
    f(t1, 1) = 1;
    for t2 = 2:4^n
        f(t1,t2) = 0;
    end
end


%% RECURSIVE STEP %%

% Now we store indices to store the left polynomial coefficient without
% taking care of index i j of xi or yj to keep the recursive step simple.
% However, we also keep track of those index wrt the step it is
% inserted by an array indices: xi index is i, yj index is j+n
indices = []; l = 0;
for step = 1:(m-1)
    i = 0; j = 0;       
    if length(labels{step}) == 4
        i = str2num(labels{step}(2));
        j = str2num(labels{step}(4));
        indices(l+1) = i;
        indices(l+2) = j+n;
    elseif length(labels{step}) == 2
        if labels{step}(1) == 'x'
            i = str2num(labels{step}(2));
            indices(l+1) = i;
        else
            j = str2num(labels{step}(2));
            indices(l+1) = j+n;
        end
    end
    curr = mod(step, 2);
    if i == 0 && j == 0
        f(curr+1,1:2^l) = probs(step)/(probs(step)-1)*f(curr+1,1:2^l)+ ...
            1/(1-probs(step))*f(2-curr,1:2^l);
    elseif i ~= 0 && j ~= 0
        f(curr+1,1+3*2^l:2^(l+2)) = f(curr+1,1:2^l)*probs(step)/(probs(step)-1);
        f(curr+1,1:2^l) = f(2-curr,1:2^l)*1/(1-probs(step));
        l = l+2;
    else
        f(curr+1,(1+2^l):2^(l+1)) = f(curr+1,1:2^l)*probs(step)/(probs(step)-1);
        f(curr+1,1:2^l) = f(2-curr,1:2^l)*1/(1-probs(step));
        l = l+1;
    end
end

disp(f);
 

%% BUILDING MATRIX STEP %%
z = sym('z%d', [1 n]);
M = sym([2^n 2^n]);
for i1 = 1:2^n
    for j1 = 1:2^n
        M(i1, j1) = 0;
    end
end

for YCoef = 0:(2^n-1)
    for ind = 1:4^n
        [X, Y] = trueIndex(indices, ind-1, n);
        [coef, I] = mult(z, X, Y, YCoef, n);
        M(I+1, YCoef+1) = simplify(M(I+1, YCoef+1)) + f(2, ind)*coef;
    end
end
%% DEBUG %%
disp(M)

%% Solving for sum Q(z) and take nth partial derivatives to get result r %%
r = det(M);
for i1 = 1:2^n
    M(1, i1) = 1;
end
r = det(M)/r;

%{
for i = 1:n
    r = subs(diff(r, z(i)), z(i), 1);
end
%}

%% RESULT %%
r = simplify(r)
%% HELPER FUNCTIONS %%
function [X,Y] = trueIndex(lookup, index, k)
    X = 0; Y = 0;
    for i = 1:2*k
        if bitget(index,i)
            if(lookup(i) <= k)
                X = bitset(X, lookup(i));
            else
                Y = bitset(Y, lookup(i)-k);
            end
        end
    end     
end
% If we have some prodx_i (X) * prody_j (Y) this fct return the 
% coefficient and the corresponding entry when we multiply the above product
% with a specific prod_yk
function [coef, I] = mult(z, X, Y, YCoef, k)
    coef = 1;
    U = bitor(Y, YCoef);
    V = bitand(X,U);
    I = U-V;
    for i = 1:k
        if bitget(V, i)
            coef = coef * z(i);
        end
    end
end