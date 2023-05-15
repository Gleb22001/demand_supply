clear
cd 'C:\Users\123\Desktop\Проекты\Диплом 2\Расчёты'
addpath('functions');
 %% Грузим даннные

filename="Продовольственные товары_декомпозиции.xlsx";
for i=1:17
[hd_q, hd_p]=fun_uni(2*i, 2*i-1, "Продовольственные товары.xlsx");
writematrix(hd_q,filename,'Sheet',2*i-1,'Range','B1');
writematrix(hd_p,filename,'Sheet',2*i,'Range','B1');
i
end


filename="Неродовольственные товары_декомпозиции.xlsx";
for i=1:12
[hd_q, hd_p]=fun_uni(2*i, 2*i-1, "Непродовольственные товары.xlsx");
writematrix(hd_q,filename,'Sheet',2*i-1,'Range','B1');
writematrix(hd_p,filename,'Sheet',2*i,'Range','B1');
i
end



filename="Услуги_декомпозиции.xlsx";
for i=1:13
[hd_q, hd_p]=fun_uni(2*i, 2*i-1, "Услуги2.xlsx");
writematrix(hd_q,filename,'Sheet',2*i-1,'Range','B1');
writematrix(hd_p,filename,'Sheet',2*i,'Range','B1');
i
end

filename="Регионы_декомпозиции1.xlsx";
for i=1:8
[hd_q, hd_p]=fun_uni(2*i, 2*i-1, "Регионы.xlsx");
writematrix(hd_q,filename,'Sheet',2*i-1,'Range','B1');
writematrix(hd_p,filename,'Sheet',2*i,'Range','B1');
i
end

