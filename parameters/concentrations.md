# Species Concentrations

| Species                      | Name in model | Concentration (nM) | Concentration (mol./cell) | Reference                                                                             |
|------------------------------|---------------|---------------------|---------------------------|---------------------------------------------------------------------------------------|
| free Cdc20                   | C             | 10                  | 25                        | C=Ctot-AC-MC-ACMCx2                                                                    |
| free APC                     | A             | 32                  | 80                        | FCCS experiments (unpublished)                                                         |
| APC/CCdc20                   | AC            | 11.2                | 28                        | AC=Atot-A-ACMC                                                                         |
| Mad comp. (active)           | Mstar         | 42.4                | 106                       | Mstar=Mtot-MC-ACMC-M                                                                   |
| Mad comp. (inactive)         | M             | 0                   | 0                         | 0 during an arrest                                                                     |
| Cyclin B                     | CycB          | 172                 | 430                       | mass spectrophtometry data [[1]](#matafora_2017)                                        |
| Checkpoint core complex      | MC            | 15.2                | 38                        | FCCS experiments (unpublished)                                                         |
| APC:MCC                      | ACMC          | 12                  | 30                        | FCCS experiments (unpublished)                                                         |
| APC total                    | Atot          | 60                  | 150                       | FCCS experiments (unpublished)                                                         |
| Mad total                    | Mtot          | 70                  | 175                       | from Cdc23 concentration in [[2]](#currbio) which is a subunit of APC/C          |
| Cdc20 total                  | Ctot          | 60                  | 150                       | limiting species Mad3 [[2]](#currbio)                                            |

## References

- <a name="matafora_2017"></a>[1] Matafora V, Corno A, Ciliberto A, Bachi A. Missing Value Monitoring Enhances the Robustness in Proteomics Quantitation. J Proteome Res. 2017;16(4):1719-1727. doi:10.1021/acs.jproteome.6b01056
- <a name="currbio"></a>[2] Bonaiuti P, Chiroli E, Gross F, et al. Cells Escape an Operational Mitotic Checkpoint through a Stochastic Process. Curr Biol. 2018;28(1):28-37.e7. doi:10.1016/j.cub.2017.11.031
