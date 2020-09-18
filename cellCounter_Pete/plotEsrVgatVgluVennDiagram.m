%plot PMC Esr Venn diagram
% For three-circle venn diagrams, A is a three element vector [c1 c2 c3], 
% and I is a four element vector [i12 i13 i23 i123], specifiying the 
% two-circle intersection areas i12, i13, i23, and the three-circle
% intersection i123.

Esr = 1;
Vgat = 0.05;
Vglut = 0.88;

A = [300 200]; I = 150;
venn(A,I,'FaceColor',{'r','y'},'FaceAlpha',{1,0.6},'EdgeColor','black')

