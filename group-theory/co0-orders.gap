indivorder := function(a,b)
	local g;
	g:=Group(a,b);
	return [Order(a),Order(b),Size(g)];
end;

orders:=[];
for i in [1..167] do
	orders[i]:=List(cr[i], x-> indivorder(conjrep[i],x));
	filename:=Concatenation("orders-",String(i),".txt");
	PrintTo(filename,orders[i]);
od;


quit;
EOF
