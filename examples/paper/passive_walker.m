clr;

% load preset fem
fem = preset.fem.suzumori_walker('n',1250);

% solve dynamic simulation
fem = fem.simulate();

