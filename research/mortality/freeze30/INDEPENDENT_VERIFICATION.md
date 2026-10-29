# Independent freeze30 aggregate and synthetic-scoring audit

106/106 checks passed;84 scientific files checked and80 historical hashes preserved.

Artifact SHA256 `70d29bf6e60d0967a11534634ebd1769a8b513b97ecb1fe8c66c6a5392534d43`.
Manifest SHA256 `f33df537280ead416406139412c4155f05ba8215e37408873f8442df7c4ab488`.

Synthetic reference maximum absolute differences: {"routine": 1.3877787807814457e-17, "routine_all3": 1.1102230246251565e-16}.

- No participant records, raw archive requests or model refits were performed. Whole-PSU membership and training fit equations were reviewed in code; aggregate counts alone cannot independently reconstruct participant-level splits, gradients or training Brier.
- The two artifacts and manifests were compared byte-for-byte against both outputs of an independently fetched successful CI job. This is same-environment deterministic replay, not a cross-platform guarantee.
- The reference scoring equations and raw-coefficient reconstruction use stdlib math on six explicitly invented fixtures. The adapter was called only on those synthetic fixtures; no training participants were scored by this reviewer.
- All13 scientific/documentation source hashes match the preserved29 hashes, and the recorded current source digest was independently recomputed. Source bytes were not downloaded locally. The new full retrieval-manifest bytes are not exported, so its SHA was cross-checked between artifact/source audit/manifest rather than independently rehashed from that unexported file.
- Zero maximum training replay differences are successful-run aggregate assertions, not local participant-level replay.
- Qualification results replay metadata assertions only. Empty external mapping registry, zero qualified candidates, historical assay limits and lawful-access requirements remain; these artifacts are not independent external or clinical validation.
