Read("co0-g1g2.txt");
Co0:=Group(g1,g2);
cc:=ConjugacyClasses(Co0);
conjrep:=List(cc,Representative);
centralizers:=[];
centralizerconjugacyclasses:=[];
for j in [2..167] do
	i:=169-j;
	Display(i);
	centralizers[i]:=Centralizer(Co0,conjrep[i]);
od;
centralizers[1]:=Co0;

for i in [2..167] do
	Display(i);
	centralizerconjugacyclasses[i]:=ConjugacyClasses(centralizers[i]);
od;

centralizerconjugacyclasses[1]:=cc;
sizes:=List(cc,Size);
sizesizes:=List(centralizerconjugacyclasses,x->List(x,Size));
centralizerconjugacyclassreps:=	List(centralizerconjugacyclasses,x->List(x,Representative));

ba:=[110989, 100395, 99081, 98574, 98511, 98454, 98409, 98363, 98336, 98325, 98320, 98306, 98303, 98298, 98297, 98296, 98295, 98292, 98290, 98289, 98286, 98284, 98282, 98281];

images:=List(centralizerconjugacyclassreps,x->List(x,y->List(ba,z->z^y));

PrintTo("co0sizes.txt",sizesizes);
PrintTo("co0images.txt",images);

SaveWorkspace("conway0");
quit;
