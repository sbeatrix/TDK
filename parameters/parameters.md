# Model Parameters

| Description                                   | Parameter  | Value                    | Reference                                                 |
|-----------------------------------------------|------------|--------------------------|-----------------------------------------------------------|
| Cdc20 degradation rate                            | kdeg       | 1 min<sup>-1</sup>       | based on Cdc20 degradation kinetics [[1]](#currbio) |
| Cdc20 synthesis rate                              | ksyn_c20   | 37 min<sup>-1</sup>      | fitted to have total Cdc20 60 nM in wild-type cells       |
| MCC formation rate                                | kassMC     | 0.04 (mol. \cdot min)<sup>-1</sup> | in the range of Mad2 dimerization [[2]](#simonetta_influence_2009) |
| MCC dissociation rate                             | kDMC       | 28.5 mol.                | kDMC = kdissmc/kassmc [[3]](#faesen_basis_2017)   |
| APCCdc20 formation rate                           | kassAC     | 0.02 (mol. \cdot min)<sup>-1</sup> | in the same range of APC/C-Ubch10 assoc rate [[4]](#chang_molecular_2014) |
| APCCdc20 dissociation rate                        | kDAC       | 28.5 mol.                | kDAC = kdissmc/kassmc [[3]](#faesen_basis_2017)   |
| ACMC formation rate                               | kassACMC   | 0.04 (mol. \cdot min)<sup>-1</sup> | in the same range of APC/C-Ubch10 assoc rate [[4]](#chang_molecular_2014) |
| ACMC dissociation rate                            | kDACMC     | 10 mol.                  | based on [[5]](#foster_apc/c_2012)                |
| Background degradation rate for Cdc20, MC, AC, ACMC| kdegBG     | 0.02 min<sup>-1</sup>    | fitted to be ~10x smaller than main degradation           |
| Activation rate of M                               | kact       | 0.0015 min<sup>-1</sup>  | -                                                         |
| Inactivation rate of M                             | kinact     | 0.45 (mol. \cdot min)<sup>-1</sup> | -                                                         |
| Cyclin B synthesis rate                           | ksyn_cy    | 317 (mol. \cdot min)<sup>-1</sup> | fitted to have total Cyclin B 400 mol./cell in wild-type cells |
| Cyclin B degradation rate                         | kdeg_cy    | 0.024 (mol. \cdot min)<sup>-1</sup> | fitted to have total Cyclin B 400 mol./cell in wild-type cells |
| Cyclin B background degradation rate               | kdegbg_cy  | 0.068 min<sup>-1</sup>   | fitted to be ~10x smaller than main degradation           |
| Kinetochore attachment rate                   | knUK_att   | 0.01 min<sup>-1</sup>    | fitted to simulate drug washout time interval             |
| Saturation function for nUK signal                             | k_nuk      | 78                       | from [[6]](#joglekar)                           |
| Saturation function for nUK signal                             | J1         | 0.5                      | from [[6]](#joglekar)                           |
| M activation and inactivation: MM constant    | J          | 0.06                     | -                                                         |

## References

- <a name="currbio"></a>[1] Bonaiuti P, Chiroli E, Gross F, et al. Cells Escape an Operational Mitotic Checkpoint through a Stochastic Process. Curr Biol. 2018;28(1):28-37.e7. doi:10.1016/j.cub.2017.11.031
- <a name="simonetta_influence_2009"></a>[2] Simonetta M, Manzoni R, Mosca R, et al. The influence of catalysis on mad2 activation dynamics. PLoS Biol. 2009;7(1):e10. doi:10.1371/journal.pbio.1000010
- <a name="faesen_basis_2017"></a>[3] Faesen AC, Thanasoula M, Maffini S, et al. Basis of catalytic assembly of the mitotic checkpoint complex. Nature. 2017;542(7642):498-502. doi:10.1038/nature21384
- <a name="chang_molecular_2014"></a> [4] Chang LF, Zhang Z, Yang J, McLaughlin SH, Barford D. Molecular architecture and mechanism of the anaphase-promoting complex. Nature. 2014;513(7518):388-393. doi:10.1038/nature13543
- <a name="foster_apc/c_2012"></a>[5] Foster SA, Morgan DO. The APC/C subunit Mnd2/Apc15 promotes Cdc20 autoubiquitination and spindle assembly checkpoint inactivation. Mol Cell. 2012;47(6):921-932. doi:10.1016/j.molcel.2012.07.031
- <a name="joglekar"></a>[6] Aravamudhan P, Chen R, Roy B, Sim J, Joglekar AP. Dual mechanisms regulate the recruitment of spindle assembly checkpoint proteins to the budding yeast kinetochore. Mol Biol Cell. 2016;27(22):3405-3417. doi:10.1091/mbc.E16-01-0007


