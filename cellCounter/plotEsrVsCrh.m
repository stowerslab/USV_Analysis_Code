%plot PMC Esr vs. Crh sections
esrMeanCounts = fliplr([16 34.16666667 53.83333333 64.33333333 35 15.33333333]);
esrStdCounts = fliplr([11.9498954 16.58211888 9.042492282 6.976149845 14.5327217 11.91077943]);
crhMeanCounts = fliplr([78.66666667 100.6666667 109.5 121 88.16666667 45.66666667]);
crhStdCounts = fliplr([14.94880151 12.4365054 11.53689733 10.67707825 24.83076049 22.16002407]);
overlapMeanCounts = fliplr([5.666666667 9.333333333 10.33333333 10.83333333 5.5 3.166666667]);
overlapStdCounts = fliplr([5.202563471 6.889605697 4.844240567 2.483277404 3.271085447 2.316606714]);

fontSz = 20;

hPmcSections = figure;
hold on;
errorbar([-125 -75 -25 25 75 125], esrMeanCounts, esrStdCounts, 'mx', 'LineWidth', 2);
plot([-125 -75 -25 25 75 125], esrMeanCounts, 'o-', 'Color', 'm', 'MarkerSize', 5, 'LineWidth', 3);
errorbar([-125 -75 -25 25 75 125], crhMeanCounts, crhStdCounts, 'cx', 'LineWidth', 2);
plot([-125 -75 -25 25 75 125], crhMeanCounts, 'o-', 'Color', 'c', 'MarkerSize', 5, 'LineWidth', 3);
errorbar([-125 -75 -25 25 75 125], overlapMeanCounts, overlapStdCounts, 'kx', 'LineWidth', 2);
plot([-125 -75 -25 25 75 125], overlapMeanCounts, 'o-', 'Color', 'k', 'MarkerSize', 5, 'LineWidth', 3);
hold off;
axis([-140 140 0 140]);
set(gca, 'FontSize', fontSz);
set(gca,'XTick', [-125 -75 -25 25 75 125], 'XTickLabel',{'-125' '-75' '-25' '25' '75' '125'})
ylabel('# cells', 'FontSize', fontSz)
xlabel('rostral-caudal distance from PMC center (um)', 'FontSize', fontSz)




