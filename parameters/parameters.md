# Model Parameters

| Description                                   | Parameter  | Value                    | Reference                                                 |
|-----------------------------------------------|------------|--------------------------|-----------------------------------------------------------|
| Cdc20 degradation rate                            | kdeg       | 1 min<sup>-1</sup>       | based on Cdc20 degradation kinetics [[20]](#currbio) |
| Cdc20 synthesis rate                              | ksyn_c20   | 37 min<sup>-1</sup>      | fitted to have total Cdc20 60 nM in wild-type cells       |
| MCC formation rate                                | kassMC     | 0.04 (mol. \cdot min)<sup>-1</sup> | in the range of Mad2 dimerization [[33]](#simonetta_influence_2009) |
| MCC dissociation rate                             | kDMC       | 28.5 mol.                | kDMC = kdissmc/kassmc [[13]](#faesen_basis_2017)   |
| APCCdc20 formation rate                           | kassAC     | 0.02 (mol. \cdot min)<sup>-1</sup> | in the same range of APC/C-Ubch10 assoc rate [[34]](#chang_molecular_2014) |
| APCCdc20 dissociation rate                        | kDAC       | 28.5 mol.                | kDAC = kdissmc/kassmc [[13]](#faesen_basis_2017)   |
| ACMC formation rate                               | kassACMC   | 0.04 (mol. \cdot min)<sup>-1</sup> | in the same range of APC/C-Ubch10 assoc rate [[34]](#chang_molecular_2014) |
| ACMC dissociation rate                            | kDACMC     | 10 mol.                  | based on [[35]](#foster_apc/c_2012)                |
| Background degradation rate for Cdc20, MC, AC, ACMC| kdegBG     | 0.02 min<sup>-1</sup>    | fitted to be ~10x smaller than main degradation           |
| Activation rate of M                               | kact       | 0.0015 min<sup>-1</sup>  | -                                                         |
| Inactivation rate of M                             | kinact     | 0.45 (mol. \cdot min)<sup>-1</sup> | -                                                         |
| Cyclin B synthesis rate                           | ksyn_cy    | 317 (mol. \cdot min)<sup>-1</sup> | fitted to have total Cyclin B 400 mol./cell in wild-type cells |
| Cyclin B degradation rate                         | kdeg_cy    | 0.024 (mol. \cdot min)<sup>-1</sup> | fitted to have total Cyclin B 400 mol./cell in wild-type cells |
| Cyclin B background degradation rate               | kdegbg_cy  | 0.068 min<sup>-1</sup>   | fitted to be ~10x smaller than main degradation           |
| Kinetochore attachment rate                   | knUK_att   | 0.01 min<sup>-1</sup>    | fitted to simulate drug washout time interval             |
| Saturation function for nUK signal                             | k_nuk      | 78                       | from [[25]](#joglekar)                           |
| Saturation function for nUK signal                             | J1         | 0.5                      | from [[25]](#joglekar)                           |
| M activation and inactivation: MM constant    | J          | 0.06                     | -                                                         |

## References

- <a name="currbio"></a>[20] Bonaiuti P, Chiroli E, Gross F, et al. Cells Escape an Operational Mitotic Checkpoint through a Stochastic Process. Curr Biol. 2018;28(1):28-37.e7. doi:10.1016/j.cub.2017.11.031
- <a name="simonetta_influence_2009"></a>[33] Simonetta M, Manzoni R, Mosca R, et al. The influence of catalysis on mad2 activation dynamics. PLoS Biol. 2009;7(1):e10. doi:10.1371/journal.pbio.1000010
- <a name="faesen_basis_2017"></a>[13] Faesen AC, Thanasoula M, Maffini S, et al. Basis of catalytic assembly of the mitotic checkpoint complex. Nature. 2017;542(7642):498-502. doi:10.1038/nature21384
- <a name="chang_molecular_2014"></a> [34] Chang LF, Zhang Z, Yang J, McLaughlin SH, Barford D. Molecular architecture and mechanism of the anaphase-promoting complex. Nature. 2014;513(7518):388-393. doi:10.1038/nature13543
- <a name="foster_apc/c_2012"></a>[35] Foster SA, Morgan DO. The APC/C subunit Mnd2/Apc15 promotes Cdc20 autoubiquitination and spindle assembly checkpoint inactivation. Mol Cell. 2012;47(6):921-932. doi:10.1016/j.molcel.2012.07.031
- <a name="joglekar"></a>[25] Aravamudhan P, Chen R, Roy B, Sim J, Joglekar AP. Dual mechanisms regulate the recruitment of spindle assembly checkpoint proteins to the budding yeast kinetochore. Mol Biol Cell. 2016;27(22):3405-3417. doi:10.1091/mbc.E16-01-0007

# Species Concentrations

| Species                      | Name in model | Concentration (nM) | Concentration (mol./cell) | Reference                                                                             |
|------------------------------|---------------|---------------------|---------------------------|---------------------------------------------------------------------------------------|
| free Cdc20                   | C             | 10                  | 25                        | C=Ctot-AC-MC-ACMCx2                                                                    |
| free APC                     | A             | 32                  | 80                        | FCCS experiments (unpublished)                                                         |
| APC/CCdc20                   | AC            | 11.2                | 28                        | AC=Atot-A-ACMC                                                                         |
| Mad comp. (active)           | Mstar         | 42.4                | 106                       | Mstar=Mtot-MC-ACMC-M                                                                   |
| Mad comp. (inactive)         | M             | 0                   | 0                         | 0 during an arrest                                                                     |
| Cyclin B                     | CycB          | 172                 | 430                       | mass spectrophtometry data [[28]](#matafora_2017)                                        |
| Checkpoint core complex      | MC            | 15.2                | 38                        | FCCS experiments (unpublished)                                                         |
| APC:MCC                      | ACMC          | 12                  | 30                        | FCCS experiments (unpublished)                                                         |
| APC total                    | Atot          | 60                  | 150                       | FCCS experiments (unpublished)                                                         |
| Mad total                    | Mtot          | 70                  | 175                       | from Cdc23 concentration in [[20]](#currbio) which is a subunit of APC/C          |
| Cdc20 total                  | Ctot          | 60                  | 150                       | limiting species Mad3 [[20]](#currbio)                                            |

## References

- <a name="matafora_2017"></a>[28] Matafora V, Corno A, Ciliberto A, Bachi A. Missing Value Monitoring Enhances the Robustness in Proteomics Quantitation. J Proteome Res. 2017;16(4):1719-1727. doi:10.1021/acs.jproteome.6b01056
- <a name="currbio"></a>[20] Bonaiuti P, Chiroli E, Gross F, et al. Cells Escape an Operational Mitotic Checkpoint through a Stochastic Process. Curr Biol. 2018;28(1):28-37.e7. doi:10.1016/j.cub.2017.11.031
