# Evaluated PARSEC common-radiation raw cleanup

Operator's standing instruction and renewed instruction on2026-09-10:
remove test raw outputs after successful test AND stage evaluation.
Both runs ended; all82 physical datasets exactly matched. Inputs, logs,
binary/source identities, physical input packages and compact evaluation
remain. Only these three exact HDF5 targets are authorized here.

Root `/gpfs/kjhan/LRD_JWST/.parsec-sed.ohf9MX/`.

| Relative target | bytes | SHA256 before deletion |
| --- | --- | --- |
| live/output_00001/data_00001.h5 | 12649888 | 067fc563cb87e9fa23c02734c671dc3fa049538afb886707e98c70ed10c8f645 |
| live/output_00002/data_00002.h5 | 12649888 | dea558a659c4fac306fda6cab61b97be74f2e9769b9594e8f8922921712cb7fa |
| restart/output_00002/data_00002.h5 | 12649888 | c774e47cdbc56dfa2e42b6cf28eaf7c06a5d42c0240a2c9c247f65368b263c86 |

Total37949664bytes. `restart/output_00001` is only a symlink to the first
fresh checkpoint, not a fourth dump. These raw files cannot be recovered
from trash after direct deletion; they can be regenerated using retained
inputs/build. Do not use these directories as active restart checkpoints.

Status: DONE. All three exact targets removed; absence confirmed.
