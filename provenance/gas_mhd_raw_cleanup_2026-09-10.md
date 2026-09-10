# Evaluated MHD raw-output cleanup

Completed removal: 23 raw HDF5 files, 1,504,265,912 bytes (1.401 GiB).
The local test folder retains 46 MiB of inputs/logs/executables/metadata.
This is not a cleanup of unrelated project runs or scientific production data.

Final cleanup: bundle 4 native parity and driver evaluation completed; the
previously retained comparison snapshots and unused duplicate restart fixture
are now removed too. No active test depends on them. Results and exact final
binary identities are in the completion plan. Input GRAFIC files are retained.

Final removal inventory (SHA256, then byte sizes):

```
921c8866df27ae45eb52a7e385c89dcb159ac2b82717dc94fd506b78d0712196  .mhd-completion.MDLk0E/alfven32/output_00001/data_00001.h5
2860ab1974d77640e9f960e88dc50bbb8e9de7dd5f0d88276341a38538d6cd5e  .mhd-completion.MDLk0E/integrated-amr-static-ranks/output_00002/data_00002.h5
20c77bdfc21a25f529c635769cebed643a6fdec1c5d6b3419ab21dc3bb954981  .mhd-completion.MDLk0E/alfven-hybrid/output_00001/data_00001.h5
f5b104ba3dfe69612b1ac15abfe6defa12a635b8c4d98e3e96b6fa157871acd0  .mhd-completion.MDLk0E/mpi1-restart/output_00001/data_00001.h5
0d4901bb3b093c9f7acce59af69b6c49a4eb3f067995ec3f42fd08652d07fe94  .mhd-completion.MDLk0E/gas-remap/output_00002/data_00002.h5
22f5cfe4c2ca4fe0cb7a991075bbbfc3fb2e8c441c8f4568945f9e0fa35b10a9  .mhd-completion.MDLk0E/integrated-hybrid/output_00001/data_00001.h5
9644585239299c4e95da06357407c3becfbf86cf1485e611b35fdd3653ac2627  .mhd-completion.MDLk0E/gas-omp4/output_00001/data_00001.h5
d1185f0301f575a70e361550b3c6fd2bf5f2cadd38ff7a05b3d77aedb6b24038  .mhd-completion.MDLk0E/integrated-omp4/output_00001/data_00001.h5
fe570c6d6be65f3269eb87f0c2d3ad408155cb8566219d02157c77026d6082fe  .mhd-completion.MDLk0E/gas-hybrid/output_00001/data_00001.h5
7673936 .mhd-completion.MDLk0E/alfven32/output_00001/data_00001.h5
204740200 .mhd-completion.MDLk0E/integrated-amr-static-ranks/output_00002/data_00002.h5
11292928 .mhd-completion.MDLk0E/alfven-hybrid/output_00001/data_00001.h5
1801552 .mhd-completion.MDLk0E/mpi1-restart/output_00001/data_00001.h5
1801552 .mhd-completion.MDLk0E/gas-remap/output_00002/data_00002.h5
204740200 .mhd-completion.MDLk0E/integrated-hybrid/output_00001/data_00001.h5
2643456 .mhd-completion.MDLk0E/gas-omp4/output_00001/data_00001.h5
204740200 .mhd-completion.MDLk0E/integrated-omp4/output_00001/data_00001.h5
2643456 .mhd-completion.MDLk0E/gas-hybrid/output_00001/data_00001.h5
```

Operator requested cleanup after test AND stage evaluation. Raw HDF5 files
below are removed permanently; preserve effective namelists, logs, executable
identity and results in gas_mhd_completion_plan_2026-09-10.md. Reproduction
requires rerunning the retained inputs. Keep Alfven32, gas-remap final and
integrated-AMR final until the active OpenMP/CUDA comparison ends.

Bundle 3 evaluation before cleanup: both mixed-AMR runs completed with
positive thermal energy, all 36 stored fields finite. Coupled restart has
7194 particles and coordinate-sorted final hydro/face fields exactly equal
to uninterrupted run on every level 1..4. Minimum heat 2.6091578904e-8,
max divB 2.3852447795e-18. MPI 2->1 gas restart differs by at most
2.2204460493e-16 from uninterrupted run; max divB 5.3290705182e-15.
Gas remapping executes at steps 2 and 4; persistent-IR remapping remains
explicitly unsupported. No new physical-input approval is inferred.

SHA256 and exact sizes/targets, relative to /gpfs/kjhan/LRD_JWST:

```
cf0e657eeba911552be56157733fc05e4bea58ccb4ffaeb0d4a3539b88207e98  .mhd-completion.MDLk0E/stellar-dust-cpu/output_00001/data_00001.h5
5c4424a2a466ca466176365e26913fc1f5b615cf323a1c52a0a1300e64a90ea2  .mhd-completion.MDLk0E/stellar-dust-cpu/output_00002/data_00002.h5
35242534764f26429133e4c4309602caad54039584f9476075a9c533da003355  .mhd-completion.MDLk0E/expansion-grafic/output_00001/data_00001.h5
869732fb9f59730a4619485e8f72373e5c389d7ddaf0d28b307dab0b70e03355  .mhd-completion.MDLk0E/zero/output_00001/data_00001.h5
d059c65bb1953bb3e12dc2506163775b49ae8f2c299f85425c88fffa3ce5e5b0  .mhd-completion.MDLk0E/stellar-agn-cpu/output_00001/data_00001.h5
fb40cda9d12524f5e40bb73dbc971491b7d3484c5fa578f37c08dda2e0e8343e  .mhd-completion.MDLk0E/integrated-amr-restart/output_00001/data_00001.h5
fb40cda9d12524f5e40bb73dbc971491b7d3484c5fa578f37c08dda2e0e8343e  .mhd-completion.MDLk0E/integrated-amr-static-ranks/output_00001/data_00001.h5
149d02f9bbfe0d4f98157c8461c08efb0a40a3c294da9316447f95747d2a4bb4  .mhd-completion.MDLk0E/stellar-agn-cpu/output_00002/data_00002.h5
f21abe24f46405ffbc5762c9189db8f941ae47d043e59438c74124a341d36180  .mhd-completion.MDLk0E/integrated-amr-restart/output_00002/data_00002.h5
fb06d158bcc8f71fa083377532794b4d34da52b7dd05c9d4fb1d3a6ff790a095  .mhd-completion.MDLk0E/brio32/output_00001/data_00001.h5
08f6e3f05c3fc6597a96ae4cd5718cd9b9810d668fa43c2d6b35ce14e08a5c62  .mhd-completion.MDLk0E/alfven16/output_00001/data_00001.h5
e876fee70d65469748462cac6ea43cf105e7a5776ded6654fdd146c4832c6083  .mhd-completion.MDLk0E/gas-remap-restart1/output_00002/data_00002.h5
c6b8a8491173f2baad0fa9d0d6bf08b294ff412f74acea237c9c96fe87de0822  .mhd-completion.MDLk0E/gas-remap-restart1/output_00001/data_00001.h5
c6b8a8491173f2baad0fa9d0d6bf08b294ff412f74acea237c9c96fe87de0822  .mhd-completion.MDLk0E/gas-remap/output_00001/data_00001.h5
58415656 .mhd-completion.MDLk0E/stellar-dust-cpu/output_00001/data_00001.h5
58530344 .mhd-completion.MDLk0E/stellar-dust-cpu/output_00002/data_00002.h5
166480 .mhd-completion.MDLk0E/expansion-grafic/output_00001/data_00001.h5
166480 .mhd-completion.MDLk0E/zero/output_00001/data_00001.h5
58665112 .mhd-completion.MDLk0E/stellar-agn-cpu/output_00001/data_00001.h5
204337672 .mhd-completion.MDLk0E/integrated-amr-restart/output_00001/data_00001.h5
204337672 .mhd-completion.MDLk0E/integrated-amr-static-ranks/output_00001/data_00001.h5
58739256 .mhd-completion.MDLk0E/stellar-agn-cpu/output_00002/data_00002.h5
204740200 .mhd-completion.MDLk0E/integrated-amr-restart/output_00002/data_00002.h5
7673936 .mhd-completion.MDLk0E/brio32/output_00001/data_00001.h5
1010968 .mhd-completion.MDLk0E/alfven16/output_00001/data_00001.h5
1801552 .mhd-completion.MDLk0E/gas-remap-restart1/output_00002/data_00002.h5
1801552 .mhd-completion.MDLk0E/gas-remap-restart1/output_00001/data_00001.h5
1801552 .mhd-completion.MDLk0E/gas-remap/output_00001/data_00001.h5
```
