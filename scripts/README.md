Usage
------

The scripts in this directory show the minimal working examples typically
used by astronomers.

run_nupipeline.sh is the base call that run nupipleine from Stage 1 to
Stage 2 (calibrated event files). This is the standard usage.

Use of the other scripts should be done with additional testing to ensure 
that the filtering is now introducing significant electronic nosie into
the results

run_nupipeline_highrate.sh applies a more lenient electronics filtering
that can be useful for source rates in excess of 500 cps. Please see
[this calibration coffee PDF]{https://github.com/NuSTAR/calibration_coffee/blob/main/pdf/NuSTAR_CalCoffee_20241002.pdf}
for additional info.


