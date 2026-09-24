Read("co1-g1g2.txt");
Co1:=Group(g1,g2);
cc:=ConjugacyClasses(Co1);
conjrep:=List(cc,Representative);
centralizers:=[];
centralizerconjugacyclasses:=[];
for i in [2..101] do
	Display(i);
	centralizers[i]:=Centralizer(Co1,conjrep[i]);
od;
centralizers[1]:=Co1;
for i in [2..101] do
	Display(i);
	centralizerconjugacyclasses[i]:=ConjugacyClasses(centralizers[i]);
od;
centralizerconjugacyclasses[1]:=cc;
sizes:=List(cc,Size);
sizesizes:=List(centralizerconjugacyclasses,x->List(x,Size));
centralizerconjugacyclassreps:=	List(centralizerconjugacyclasses,x->List(x,Representative));

ba:=[12709, 2115, 801, 294, 231, 174, 129, 83, 56, 45, 40, 26, 23, 18, 17, 16, 15, 12, 10, 9, 6, 4, 2, 1];

images:=List(centralizerconjugacyclassreps,x->List(x,y->List(ba,z->z^y));

PrintTo("co1sizes.txt",sizesizes);
PrintTo("co1images.txt",images);

SaveWorkspace("conway1");
quit;
