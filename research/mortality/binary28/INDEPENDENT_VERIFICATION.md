# Independent binary28 aggregate and percentile audit

109/109 checks passed; 73 scientific Python files checked.

- No participant records or individual OOF predictions were accessed; point Brier/AUC cannot be re-estimated from raw records in this audit.
- All500 paired contrast rows were independently parsed; percentile endpoints were recomputed arithmetically without fitting any model.
- Recorded source fingerprints were checked against manifest and qualified27; source file bytes were not downloaded locally.
- Per-replicate fit parameters/gradients are not exported, so convergence guards are reviewed in code and point-fold diagnostics, with replica failures verified from exported status counts.
- Core calibration intercept/slope intervals use498 defined replicates with2 undefined; they are not nominal full500-draw coverage. Primary Brier/AUC contrasts have500 defined draws.
- Nested CV remains internal and conditional on fixed partition, reported assays, prespecified forms and complete profiles; it is not external or clinical validation.
