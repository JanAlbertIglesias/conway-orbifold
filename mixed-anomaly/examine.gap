cr:=centralizerconjugacyclassreps;;
findclass := function(group,cc,elem)
	local i;
	for i in [1..Length(cc)] do
		if IsConjugate(group,cc[i],elem) then
			Print(-i);
			return i;
		fi;
	od;
end;

findComPair:=function(a,b)
	local i,c,d,j;
	i:=findclass(Co1,conjrep,a);
	c:=RepresentativeAction(Co1,conjrep[i],a);
	d:=c*b*Inverse(c);
	j:=findclass(centralizers[i],cr[i],d);
	return [i,j];
end;

g:=cr[1][4];
h:=cr[4][21];
result:=List([0..12],x->List([0..25],y->findComPair(g^(2*x),h*(g^y))));
PrintTo("result.txt",result);