Read("co1-g1g2.txt");
Co1:=Group(g1,g2);
Read("ccrr.txt");
co1centralizers:=[];
co1centralizers[1]:=Co1;
for i in [2..101] do
        co1centralizers[i]:=Centralizer(Co1,ccrr[1][i]);
od;

bosh:=function(p)
	if p>98280 then
		return p-98280;
	else
		return 98281-p;
	fi;
end;

toCo1:=function(perm)
	return PermList(List([98281..196560],x->bosh(x^perm)));
end;
	
findclass := function(group,cc,elem)
	local i;
	for i in [1..Length(cc)] do
		if IsConjugate(group,cc[i],elem) then
			Print(-i);
			return i;
		fi;
	od;
	Print("shouldn't happen");
end;



i := $j;

e := toCo1(centralizerconjugacyclassreps[1][i]);
x := findclass(Co1,ccrr[1],e);
a := RepresentativeAction(Co1,ccrr[1][x],e);
inva := Inverse(a);
pairs :=[];

for k in [1..Length(centralizerconjugacyclassreps[i])] do
	q := toCo1(centralizerconjugacyclassreps[i][k]);
	s := a * q * inva;
	t := findclass(co1centralizers[x],ccrr[x],s);
	pairs[k]:=[x,t];
od;

PrintTo("pairs-$j.txt",pairs);
quit;
EOF
