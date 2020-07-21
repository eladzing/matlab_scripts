clust=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];

% no. of streams 
streamno=[4 8 7 9 5 8 6 8 7 4 5 6 1 8 5 10];

%maximal penetrating stream in % of rvir
maxpenH=[7 40 25 40 23 9 25 29 53 13 39 24 40 16 16 45];
maxpen=[7 20  25 45 23 9  9 35 40  4 39 24 40 6  49  5];

maxpen2=[17 50 52 46 62 31 40 36 66 29 42 51 49 16 53 55];

avgpen=[37 67 68 55 57 41 47 54 74 37 50 46 46 52 50 36];

maxpen06=[11 25 29 33 20 7 5 5 10 18 60 38 36 25 12 9];

%        101 102 103 104 105 106 107 003 005 006 007 009 010 011 014 024];
% 0 designates unrelaxed, 1 relaxed for each projection
penet=maxpen<=median(maxpen);


relaxedness=[0,0,0;0,0,0;0,0,0;1,1,1;0,0,1;0,0,0;1,0,0;1,1,1;...
    1,1,1;0,0,0;1,1,1;0,0,0;1,1,1;0,0,0;1,1,1;0,1,0];
unrelaxedness=~relaxedness;

urlx=any(unrelaxedness,2)';
rlx=~urlx;


rlx06=[1 0 0 1 0 0 0 0 0 0 1 0 1 0 0 0];
urlx06=~rlx06;
rlx06=~urlx06;

%rlxID=find(rlx);
%urlxID=find(urlx);

clustRLX=clust(rlx);
clustURLX=clust(urlx);