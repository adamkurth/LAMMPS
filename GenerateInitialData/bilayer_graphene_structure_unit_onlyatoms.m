function bilayer_graphene_structure_unit_onlyatoms(name,angle)

clearvars -except name n angle
close all

arot=angle/180*pi;
z0=3.34;
l0=1.42;

ID=fopen(append('data.mobile.',name),'r');
fgets(ID);
atoms=sscanf(fgets(ID),'%i');
bonds=sscanf(fgets(ID),'%i');
angles=sscanf(fgets(ID),'%i');
dihedrals=sscanf(fgets(ID),'%i');
for i=1:36
    fgets(ID);
end
x=zeros(1,atoms);
y=zeros(1,atoms);
z=zeros(1,atoms);
for i=1:atoms
    A=sscanf(fgets(ID),'%i%i%i%f%f%f%f');
    x(i)=A(5);
    y(i)=A(6);
    z(i)=A(7);
end
for i=1:3
    fgets(ID);
end
b1=zeros(1,bonds);
b2=zeros(1,bonds);
for i=1:bonds
    A=sscanf(fgets(ID),'%i%i%i%i');
    b1(i)=A(3);
    b2(i)=A(4);
end
for i=1:3
    fgets(ID);
end
a1=zeros(1,angles);
a2=zeros(1,angles);
a3=zeros(1,angles);
for i=1:angles
    A=sscanf(fgets(ID),'%i%i%i%i%i');
    a1(i)=A(3);
    a2(i)=A(4);
    a3(i)=A(5);
end
for i=1:3
    fgets(ID);
end
d1=zeros(1,dihedrals);
d2=zeros(1,dihedrals);
d3=zeros(1,dihedrals);
d4=zeros(1,dihedrals);
for i=1:dihedrals
    A=sscanf(fgets(ID),'%i%i%i%i%i%i');
    d1(i)=A(3);
    d2(i)=A(4);
    d3(i)=A(5);
    d4(i)=A(6);
end
fclose(ID);

ID=fopen(append('data.substrate.',name),'r');
fgets(ID);
satoms=sscanf(fgets(ID),'%i');
for i=1:39
    fgets(ID);
end

x=x-sum(x)/atoms;
y=y-sum(y)/atoms;

xs=zeros(1,satoms);
ys=zeros(1,satoms);
zs=zeros(1,satoms);
for i=1:satoms
    A=sscanf(fgets(ID),'%i%i%i%f%f%f%f');
    xs(i)=A(5);
    ys(i)=A(6);
    zs(i)=A(7);
end
fclose(ID);

xs=xs-sum(xs)/satoms;
ys=ys-sum(ys)/satoms;

tatoms=atoms+satoms;
X=zeros(1,tatoms);
Y=zeros(1,tatoms);
Z=zeros(1,tatoms);

X(1:atoms)=x;
Y(1:atoms)=y;
Z(1:atoms)=z0;

for i=1:satoms
    X(atoms+i)=cos(arot)*xs(i)-sin(arot)*ys(i);
    Y(atoms+i)=sin(arot)*xs(i)+cos(arot)*ys(i);
    Z(atoms+i)=zs(i);
end

l=sqrt((X(b2(1))-X(b1(1)))^2+(Y(b2(1))-Y(b1(1)))^2);

X=X/l*l0;
Y=Y/l*l0;

ID=fopen(append('../initial_data_bi_',name,'.txt'),'w');

fprintf(ID,'%s\n\n%i%s\n\n%s\n\n%8.4f%s%8.4f%s\n%8.4f%s%8.4f%s\n%8.4f%s%8.4f%s\n\n%s\n\n','Lammps Square Lattice',...
    tatoms,' atoms','1 atom types',-2000,' ',2000,' xlo xhi',-2000,' ',2000,' ylo yhi',-10,' ',10,' zlo zhi','Atoms');

for i=1:atoms
    %if abs(X(i))<n && abs(Y(i))<n
        fprintf(ID,'%12d%12d%14.7f%14.7f%14.7f%12d%12d%12d\n',i,1,X(i),Y(i),Z(i),0,0,0);
    %else
    %    fprintf(ID,'%12d%12d%12d%14.7f%14.7f%14.7f%14.7f\n',i,2,1,0.0,X(i),Y(i),Z(i));
    %end
end

for i=atoms+1:atoms+satoms
    %fprintf(ID,'%12d%12d%12d%14.7f%14.7f%14.7f%14.7f\n',i,2,X(i),Y(i),Z(i),0,0,0);
    fprintf(ID,'%12d%12d%14.7f%14.7f%14.7f%12d%12d%12d\n',i,2,X(i),Y(i),Z(i),0,0,0);
end

fprintf(ID,'\n\n%s\n\n','Velocities');

for i=1:atoms+satoms
    %if abs(X(i))<n && abs(Y(i))<n
        fprintf(ID,'%12d%14.7f%14.7f%14.7f\n',i,0,0,0);
    %else
    %    fprintf(ID,'%12d%12d%12d%14.7f%14.7f%14.7f%14.7f\n',i,2,1,0.0,X(i),Y(i),Z(i));
    %end
end

fclose(ID);
% 
% ID=fopen(append('../../RESULTS/lammps_data_init_',name,'.csv'),'w');
% 
% for i=1:atoms
%     if abs(X(i))<n && abs(Y(i))<n
%     	fprintf(ID,'%1d,%10f,%10f,%10f\n',1,X(i),Y(i),Z(i));
%     end
% end

%fclose(ID);
