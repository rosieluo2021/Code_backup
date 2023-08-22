% compare [m] genes between models
mGeneCompareTable2=cell(5,5);
mGeneCompareTable2{1,1}='cross match [m] genes';
type=fieldnames(Allmodels);
type=type(~contains(type,'constrain'));% ignore unconstrained models
% ASYN and ASYNPD
index1=contains(type,'ASYN');
group1=type(index1);
models1=group1{contains(group1,'PD')};
sets1=fieldnames(Allmodels.(models1));
models2=group1{~contains(group1,'PD')};
sets2=fieldnames(Allmodels.(models2));

% SYN and SYNPD
index2=~contains(type,'ASYN');
group2=type(index2);
models3=group2{contains(group2,'PD')};
sets3=fieldnames(Allmodels.(models3));
models4=group2{~contains(group2,'PD')};
sets4=fieldnames(Allmodels.(models4));
% model name
mGeneCompareTable2{1,2} = [sets1{1} '('  num2str(size(Allmodels.(models1).(sets1{1}).mgenes,1)) ')'];
mGeneCompareTable2{1,3} = [sets1{2} '('  num2str(size(Allmodels.(models1).(sets1{2}).mgenes,1)) ')'];
mGeneCompareTable2{1,4} = [sets3{1} '('  num2str(size(Allmodels.(models3).(sets3{1}).mgenes,1)) ')'];
mGeneCompareTable2{1,5} = [sets3{2} '('  num2str(size(Allmodels.(models3).(sets3{2}).mgenes,1)) ')'];

mGeneCompareTable2{2,1} = [sets2{1} '('  num2str(size(Allmodels.(models2).(sets2{1}).mgenes,1)) ')'];
mGeneCompareTable2{3,1} = [sets2{2} '('  num2str(size(Allmodels.(models2).(sets2{2}).mgenes,1)) ')'];
mGeneCompareTable2{4,1} = [sets4{1} '('  num2str(size(Allmodels.(models4).(sets4{1}).mgenes,1)) ')'];
mGeneCompareTable2{5,1} = [sets4{2} '('  num2str(size(Allmodels.(models4).(sets4{2}).mgenes,1)) ')'];
% add overlapped gene num
mGeneCompareTable2{2,2} = sum(ismember(Allmodels.(models2).(sets2{1}).mgenes,Allmodels.(models1).(sets1{1}).mgenes));
mGeneCompareTable2{2,3} = sum(ismember(Allmodels.(models2).(sets2{1}).mgenes,Allmodels.(models1).(sets1{2}).mgenes));
mGeneCompareTable2{2,4} = sum(ismember(Allmodels.(models2).(sets2{1}).mgenes,Allmodels.(models3).(sets3{1}).mgenes));
mGeneCompareTable2{2,5} = sum(ismember(Allmodels.(models2).(sets2{1}).mgenes,Allmodels.(models3).(sets3{2}).mgenes));

mGeneCompareTable2{3,2} = sum(ismember(Allmodels.(models2).(sets2{2}).mgenes,Allmodels.(models1).(sets1{1}).mgenes));
mGeneCompareTable2{3,3} = sum(ismember(Allmodels.(models2).(sets2{2}).mgenes,Allmodels.(models1).(sets1{2}).mgenes));
mGeneCompareTable2{3,4} = sum(ismember(Allmodels.(models2).(sets2{2}).mgenes,Allmodels.(models3).(sets3{1}).mgenes));
mGeneCompareTable2{3,5} = sum(ismember(Allmodels.(models2).(sets2{2}).mgenes,Allmodels.(models3).(sets3{2}).mgenes));

%
mGeneCompareTable2{4,2} = sum(ismember(Allmodels.(models4).(sets4{1}).mgenes,Allmodels.(models1).(sets1{1}).mgenes));
mGeneCompareTable2{4,3} = sum(ismember(Allmodels.(models4).(sets4{1}).mgenes,Allmodels.(models1).(sets1{2}).mgenes));
mGeneCompareTable2{4,4} = sum(ismember(Allmodels.(models4).(sets4{1}).mgenes,Allmodels.(models3).(sets3{1}).mgenes));
mGeneCompareTable2{4,5} = sum(ismember(Allmodels.(models4).(sets4{1}).mgenes,Allmodels.(models3).(sets3{2}).mgenes));

mGeneCompareTable2{5,2} = sum(ismember(Allmodels.(models4).(sets4{2}).mgenes,Allmodels.(models1).(sets1{1}).mgenes));
mGeneCompareTable2{5,3} = sum(ismember(Allmodels.(models4).(sets4{2}).mgenes,Allmodels.(models1).(sets1{2}).mgenes));
mGeneCompareTable2{5,4} = sum(ismember(Allmodels.(models4).(sets4{2}).mgenes,Allmodels.(models3).(sets3{1}).mgenes));
mGeneCompareTable2{5,5} = sum(ismember(Allmodels.(models4).(sets4{2}).mgenes,Allmodels.(models3).(sets3{2}).mgenes));
