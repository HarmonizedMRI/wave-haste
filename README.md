# wave-haste
- Relying on the example TSE and HASTE code by Juergen Hennig and Maxim Zaitsev, this toolbox creates wave-TSE and wave-HASTE pulse sequences and performs image reconstruction.
- Specifically:
- script_wave_tse_sequence.m: example script to generate fully-sampled TSE data w/ and w/o wave encoding which can be used as calibration scans for the wave-HASTE recon.
- script_wave_tse_recon.m: loads the raw data acquired on a head phantom to compute the wave point spread function (PSF)
