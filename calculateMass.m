function  [Mass]=calculateMass(density,Elements,Sections,L)

Mass = 0;
for i = 1:length(Elements)
    Aelem = Sections( Elements(i,3), 2);
    Mass = Mass + (Aelem*L)*density;
end



