function [p, k] = diamond_bot_curve0
p0 = ...
[1.8985    0.0003
1.8690    0.0001
1.8259         0
1.7726    0.0001
1.7125    0.0005
1.6490    0.0013
1.5855    0.0026
1.5255    0.0045
1.4724    0.0071
1.4296    0.0106
1.4005    0.0149
1.4005    0.0149
1.2540    0.0783
1.1079    0.1536
0.9647    0.2379
0.8265    0.3282
0.6956    0.4217
0.5743    0.5153
0.4649    0.6061
0.3697    0.6912
0.2909    0.7676
0.2308    0.8324
0.2308    0.8324
0.1839    0.8905
0.1437    0.9489
0.1114    1.0076
0.0877    1.0666
0.0736    1.1259
0.0702    1.1856
0.0783    1.2458
0.0989    1.3064
0.1329    1.3674
0.1813    1.4289
0.1813    1.4289
0.2939    1.5184
0.4119    1.6087
0.5333    1.6989
0.6561    1.7881
0.7783    1.8752
0.8978    1.9593
1.0126    2.0392
1.1207    2.1142
1.2201    2.1831
1.3086    2.2449
1.3086    2.2449
1.3862    2.3011
1.4547    2.3533
1.5155    2.4011
1.5699    2.4440
1.6191    2.4816
1.6645    2.5133
1.7073    2.5388
1.7488    2.5576
1.7904    2.5692
1.8333    2.5732
1.8333    2.5732
1.8561    2.5733
1.8809    2.5736
1.9078    2.5740
1.9371    2.5745
1.9691    2.5750
2.0040    2.5756
2.0421    2.5761
2.0837    2.5765
2.1291    2.5769
2.1784    2.5770
2.1784    2.5770
2.2291    2.5769
2.2780    2.5765
2.3249    2.5761
2.3695    2.5756
2.4114    2.5750
2.4504    2.5745
2.4862    2.5740
2.5183    2.5736
2.5467    2.5733
2.5708    2.5732
2.5708    2.5732
2.6137    2.5692
2.6553    2.5576
2.6969    2.5388
2.7397    2.5133
2.7851    2.4816
2.8343    2.4440
2.8886    2.4011
2.9494    2.3533
3.0180    2.3011
3.0955    2.2449
3.0955    2.2449
3.1841    2.1831
3.2834    2.1142
3.3915    2.0392
3.5063    1.9593
3.6258    1.8752
3.7480    1.7881
3.8708    1.6989
3.9923    1.6087
4.1103    1.5184
4.2228    1.4289
4.2228    1.4289
4.2712    1.3674
4.3052    1.3064
4.3258    1.2458
4.3339    1.1856
4.3305    1.1259
4.3165    1.0666
4.2928    1.0076
4.2604    0.9489
4.2203    0.8905
4.1734    0.8324
4.1734    0.8324
4.1132    0.7676
4.0344    0.6912
3.9392    0.6061
3.8298    0.5153
3.7085    0.4217
3.5777    0.3282
3.4395    0.2379
3.2962    0.1536
3.1502    0.0783
3.0036    0.0149
3.0036    0.0149
2.9746    0.0106
2.9318    0.0071
2.8786    0.0045
2.8186    0.0026
2.7552    0.0013
2.6916    0.0005
2.6315    0.0001
2.5782    0.0000
2.5351    0.0001
2.5057    0.0003] .* [3.85e1,3.5e1] + [0,7];

x  = linspace(0,1,100);
p0 = flipud(unique(p0,'rows','stable'));
p = interparc(x, p0(:,1), p0(:,2));
k = curvature(p); 
end