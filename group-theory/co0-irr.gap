for i in [2..167] do
	g:=centralizers[$j];
	i:=Irr(g);
	m:=List(i,chi->ValuesOfClassFunction(chi));
	f:=Concatenation("irr-",String(i),".txt");
	PrintTo(f,m);
od;
quit;