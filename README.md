# wave-haste
- Relying on the example TSE and HASTE code by Juergen Hennig and Maxim Zaitsev, this toolbox creates wave-TSE and wave-HASTE pulse sequences and performs image reconstruction.
- Specifically:
- script_wave_tse_sequence.m: example script to generate fully-sampled TSE data w/ and w/o wave encoding which can be used as calibration scans for the wave-HASTE recon.
- script_wave_tse_recon.m: loads the raw data acquired on a head phantom to compute the wave point spread function (PSF).
- script_wave_haste_sequence.m : creates a HASTE acquisition at R=3-fold acceleration including a GRE-based ACS scan at the end of the scan, w/ and w/o wave encoding.
- script_wave_haste_recon.m: loads the acquired HASTE data as well as the pre-computed PSF from the TSE calibration, performs Espirit coil sensitivity calibration and iterative Sense recon for wave-HASTE at R=3.
- raw data can be downloaded from:
- https://www.dropbox.com/scl/fo/ebfo8zbh2w3jezcsuyehm/AIosSCpvu0b-jsGzbIpOV_E?rlkey=v0aec0ftbuklqpvx48oprbimb&dl=0
- this includes 4 twix files, R=1 TSE w/ and w/o wave and R=3 HASTE w/ and w/o wave including ACS data.
