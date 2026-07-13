# Provenance

Authoritative base artifact:

- filename: `apple2_radia_toolkit(1).zip`
- SHA-256: `0fdf5421da768c3e3a240ce6ed643792612f8f659ce6b69755c22fe409048662`
- package version: 0.1.0

This 0.2.0 source tree is a corrected derivative implementing the eight-item
review roadmap. No published-device validation data were added. The illustrative
CSV in `examples/` is a schema example only.


## Version 0.2.1 correction

Authoritative 0.2.0 base source artifact:

- filename: `apple2_radia_toolkit_v0.2.0.zip`
- SHA-256: `93d2b9027c2f1745d1cd16d68d441220fc3c6c781dbaac14922e1272aae53212`

Version 0.2.1 changes only validation-process robustness and release metadata.
The APPLE-II geometry, analysis, optimization, and reference-comparison algorithms
are unchanged from 0.2.0.


## Version 0.2.2 correction

Authoritative 0.2.1 base source artifact:

- filename: `apple2_radia_toolkit_v0.2.1.zip`
- SHA-256: `fbdc3a1111e8bdd02b79dc0bafce16d0d2a480af197b86778d28e8c787bfe562`

Version 0.2.2 changes only validation dependency preflight and release metadata.
The APPLE-II geometry, analysis, optimization, and reference-comparison algorithms
are unchanged from 0.2.1.



## Version 0.2.3 correction

Authoritative 0.2.2 base source artifact:

- filename: `apple2_radia_toolkit_v0.2.2.zip`
- SHA-256: `ee985ed4d3dc464521593ab0d1df46a8a2d6ee07f79b55f17f2b455947b0f7e1`

The target-Mac run exposed a common-cause audit defect: `validation.log` was
hashed while the tee stream was still active, and subsequent status output
changed the file. Version 0.2.3 changes only evidence finalization and release
metadata. It closes the log before hashing, records the log digest externally,
and anchors the output manifest with a separate digest. It also corrects the
public `apple2.__version__`, which had remained at 0.2.0 despite later package
metadata. The APPLE-II geometry, analysis, optimization, and reference-comparison
algorithms are unchanged from 0.2.2.


## Version 0.2.4 correction

Authoritative 0.2.3 base source artifact:

- filename: `apple2_radia_toolkit_v0.2.3.zip`
- SHA-256: `3884f799d936e4f3f6df07d5c936533095feb7efbf94a4d233e7f2f3612cd2f6`

The target Intel Mac exposed a portability defect in the regression test rather
than in the validation script: the test replaced `PATH` with `/usr/bin` (the
location of `shasum`), while macOS provides `bash` in `/bin`. Consequently
`#!/usr/bin/env bash` could not launch and returned 127. Version 0.2.4 preserves
the host `PATH` and overrides only the dummy `PYTHON_BIN` used by the test. The
validation script and all APPLE-II geometry, analysis, optimization, field, and
reference-comparison algorithms are byte-identical to version 0.2.3.


## Version 0.2.5 correction

Authoritative 0.2.4 base source artifact:

- filename: `apple2_radia_toolkit_v0.2.4.zip`
- SHA-256: `3cb4cb68d2092fec22a69241c080ccc5f881ee9ba2ffb561d64e2409d444f5f7`

The target Intel Mac exposed an evidence-format defect: `SHA256SUMS.txt` used a
nonstandard three-column layout (`digest`, byte length, filename), so
`shasum -c` treated the byte-length text as part of each filename. Version
0.2.5 moves byte lengths to `BYTES.txt`, emits a standards-compatible checksum
manifest, and adds an executable regression check. Scientific algorithms are
unchanged from 0.2.4.
