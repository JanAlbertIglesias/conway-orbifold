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
sqaction:=[];
trySq := function(i)
	sqaction[i]:=List(cr[i], x-> findclass(centralizers[i],cr[i], x ^2) );
end;
paction:=[];
tryP := function(i)
	paction[i]:=List(cr[i], x-> findclass(centralizers[i],cr[i], x ^-1) );
end;
taction:=[];
tryT := function(i)
	taction[i]:=List(cr[i], x-> findclass(centralizers[i],cr[i], x * cr[1][i]) );
end;
zaction:=[];
tryZ := function(i)
	zaction[i]:=List(cr[i], x-> findclass(centralizers[i],cr[i], x * cr[1][2]) );
end;


findclassS := function(a,b)
	local i,c,d,j;
	i:=findclass(Co0,conjrep,b);
	c:=RepresentativeAction(Co0,conjrep[i],b);
	d:=c*Inverse(a)*Inverse(c);
	j:=findclass(centralizers[i],cr[i],d);
	return [i,j];
end;
saction:=[];
tryS:=function(i)
	Print(i," with ", Length(cr[i]), "cc's\n");
	saction[i]:=List(cr[i], x-> findclassS(conjrep[i],x) );
	Display(saction[i]);
end;

for i in [1..167] do
	tryP(i);
	trySq(i);
	tryT(i);
	tryZ(i);
	PrintTo(Concatenation("p-sq-t-z-action-",String(i),".txt"),[paction[i],sqaction[i],taction[i],zaction[i]]);
	tryS(i);
	PrintTo(Concatenation("s-action-",String(i),".txt"),saction[i]);	
od;
quit;