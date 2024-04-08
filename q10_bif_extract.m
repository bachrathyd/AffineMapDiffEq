clear all
uiopen('C:\Users\Bachrathy\Downloads\q10.fig',1)
fig = gcf;
axObjs = fig.Children
dataObjs = axObjs.Children
xH = dataObjs(1).XData
yH = dataObjs(1).YData
plot(xH,yH, 'r*')
xFup = dataObjs(2).XData
yFup = dataObjs(2).YData
hold on
plot(xFup,yFup, "r--",LineWidth=3)
xFdown = dataObjs(3).XData
yFdown = dataObjs(3).YData
hold on
plot(xFdown,yFdown, "r--",LineWidth=3)
xFix = dataObjs(4).XData
yFix = dataObjs(4).YData
hold on
plot(xFix,yFix, "r--",LineWidth=3)
xFup2 = dataObjs(5).XData
yFup2 = dataObjs(5).YData
hold on
plot(xFup2,yFup2, "r--",LineWidth=3)
xFdown2 = dataObjs(6).XData
yFdown2 = dataObjs(6).YData
hold on
plot(xFdown2,yFdown2, "r--",LineWidth=3)
save("q10_bif.mat")