function colmap = bounce(n)
colmap = ...
    [0.7882    0.1569    0.5529
    0.7893    0.1643    0.5564
    0.7903    0.1717    0.5598
    0.7914    0.1792    0.5632
    0.7925    0.1866    0.5667
    0.7935    0.1940    0.5701
    0.7946    0.2015    0.5735
    0.7956    0.2089    0.5770
    0.7967    0.2163    0.5804
    0.7978    0.2238    0.5838
    0.7988    0.2312    0.5873
    0.7999    0.2386    0.5907
    0.8009    0.2461    0.5941
    0.8020    0.2535    0.5976
    0.8030    0.2609    0.6010
    0.8041    0.2684    0.6045
    0.8052    0.2758    0.6079
    0.8062    0.2833    0.6113
    0.8073    0.2907    0.6148
    0.8083    0.2981    0.6182
    0.8094    0.3056    0.6216
    0.8104    0.3130    0.6251
    0.8115    0.3204    0.6285
    0.8126    0.3279    0.6319
    0.8136    0.3353    0.6354
    0.8147    0.3427    0.6388
    0.8157    0.3502    0.6422
    0.8168    0.3576    0.6457
    0.8178    0.3650    0.6491
    0.8189    0.3725    0.6525
    0.8200    0.3799    0.6560
    0.8210    0.3873    0.6594
    0.8221    0.3948    0.6628
    0.8231    0.4022    0.6663
    0.8242    0.4096    0.6697
    0.8252    0.4171    0.6731
    0.8263    0.4245    0.6766
    0.8274    0.4319    0.6800
    0.8284    0.4394    0.6834
    0.8295    0.4468    0.6869
    0.8305    0.4543    0.6903
    0.8316    0.4617    0.6937
    0.8326    0.4691    0.6972
    0.8337    0.4766    0.7006
    0.8348    0.4840    0.7040
    0.8358    0.4914    0.7075
    0.8369    0.4989    0.7109
    0.8379    0.5063    0.7143
    0.8390    0.5137    0.7178
    0.8400    0.5212    0.7212
    0.8411    0.5286    0.7246
    0.8422    0.5360    0.7281
    0.8432    0.5435    0.7315
    0.8443    0.5509    0.7349
    0.8453    0.5583    0.7384
    0.8464    0.5658    0.7418
    0.8474    0.5732    0.7452
    0.8485    0.5806    0.7487
    0.8496    0.5881    0.7521
    0.8506    0.5955    0.7555
    0.8517    0.6029    0.7590
    0.8527    0.6104    0.7624
    0.8538    0.6178    0.7658
    0.8548    0.6253    0.7693
    0.8559    0.6327    0.7727
    0.8570    0.6401    0.7761
    0.8580    0.6476    0.7796
    0.8591    0.6550    0.7830
    0.8601    0.6624    0.7864
    0.8612    0.6699    0.7899
    0.8622    0.6773    0.7933
    0.8633    0.6847    0.7968
    0.8644    0.6922    0.8002
    0.8654    0.6996    0.8036
    0.8665    0.7070    0.8071
    0.8675    0.7145    0.8105
    0.8686    0.7219    0.8139
    0.8696    0.7293    0.8174
    0.8707    0.7368    0.8208
    0.8718    0.7442    0.8242
    0.8728    0.7516    0.8277
    0.8739    0.7591    0.8311
    0.8749    0.7665    0.8345
    0.8760    0.7739    0.8380
    0.8770    0.7814    0.8414
    0.8781    0.7888    0.8448
    0.8792    0.7962    0.8483
    0.8802    0.8037    0.8517
    0.8813    0.8111    0.8551
    0.8823    0.8186    0.8586
    0.8834    0.8260    0.8620
    0.8844    0.8334    0.8654
    0.8855    0.8409    0.8689
    0.8866    0.8483    0.8723
    0.8876    0.8557    0.8757
    0.8887    0.8632    0.8792
    0.8897    0.8706    0.8826
    0.8908    0.8780    0.8860
    0.8918    0.8855    0.8895
    0.8929    0.8929    0.8929];

if nargin>0
    colmap = clrmapping(colmap,n);
end

end

function colmap = clrmapping(colmap,arg)

if arg > 1
    x = linspace(0,1,length(colmap));
    xq = linspace(0,1,arg);
    colmap = interp1(x,colmap,xq);
elseif arg < 0
    if abs(arg) == 1, arg = 100; end
    x = linspace(0,1,length(colmap));
    xq = linspace(0,1,abs(arg));
    colmap = interp1(x,colmap,xq);
    colmap = flipud(colmap);
elseif arg == 0
    colmap = [colmap;flipud(colmap)];
end

end