using LinearAlgebra
using SparseArrays
# using Laplacians
using Random
# using Arpack

function lap(G :: Graph)
    F = zeros(G.n, G.n);
    for i = 1 : G.m
        F[G.u[i], G.v[i]] -= 1
        F[G.v[i], G.u[i]] -= 1
        F[G.u[i], G.u[i]] += 1
        F[G.v[i], G.v[i]] += 1
    end
    return F
end

function lapsp(G :: Graph)
	d=zeros(G.n);
	for i=1:G.m
		x=G.u[i];y=G.v[i];
		d[x]+=1;d[y]+=1;
	end
	uu=zeros(2*G.m+G.n);
	vv=zeros(2*G.m+G.n);
	ww=zeros(2*G.m+G.n);
	a=zeros(G.n);
	for i=1:G.n
		a[i]=i;
	end
	uu[1:G.m]=G.u;uu[G.m+1:2*G.m]=G.v;
	uu[2*G.m+1:2*G.m+G.n]=a;
	vv[1:G.m]=G.v;vv[G.m+1:2*G.m]=G.u;
	vv[2*G.m+1:2*G.m+G.n]=a;
	ww[1:2*G.m].=-1;ww[2*G.m+1:2*G.m+G.n]=d;
    return sparse(uu,vv,ww)
end


function adjsp(G :: Graph)
	uu=zeros(2*G.m);
	vv=zeros(2*G.m);
	ww=zeros(2*G.m);
	uu[1:G.m]=G.u;uu[G.m+1:2*G.m]=G.v;
	vv[1:G.m]=G.v;vv[G.m+1:2*G.m]=G.u;
	ww[1:2*G.m].=1;
	return sparse(uu,vv,ww)
end
function mppinv(L)
    n = size(L,2)
	J1=ones(n,n)/n;
    # let L = copy(L)
    #     L .-= 1.0 / n
    #     L .= inv(L)
    #     L .+= 1.0 / n
    #     return L
    # end
	invL=inv(L+J1)-J1;
	return invL
end
