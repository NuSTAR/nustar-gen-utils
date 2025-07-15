Basic Usage
-----------

The scripts in this directory show the minimal working examples typically
used by astronomers.

[run_nupipeline.sh](run_nupipeline.sh) is the base call that run nupipleine from Stage 1 to
Stage 2 (calibrated event files). This is the standard usage.

The sytnax here assumes that you are in a directory like:

```
$HOME/science/nustar/CasA
```

...and that you have unpacked an observation from the HEASARC. In this case, copy
`run_nupipeline.sh` to the parent directory:

```
cd $HOME/science/nustar/CasA
cp $HOME/sciece/local/git/nustar-gen-utils/scripts/run_nupipeline.sh .
```



NB: You may need to change the #! shell invocation at the top to make it work in your preferred
shell.

Syntax
-------

You may invoke the script by doing this:

```
cd $HOME/science/nustar/CasA
./run_nupipeline.sh 40021011002
```

...which will put output and all log files into this directory:

```
$HOME/science/nustar/CasA/40021011002/event_cl
```

Note that this may take a few minutes to run based on the length of your observation
and the speed of your computer.

Advanced Usage
---------------

Use of the other scripts should be done with additional testing to ensure 
that the filtering is now introducing significant electronic nosie into
the results


[run_nupipeline_highrate.sh](run_nupipeline_highrate.sh) applies a more lenient
electronics filtering that can be useful for source rates in excess of 500 cps.
Please see
[this calibration coffee PDF](https://github.com/NuSTAR/calibration_coffee/blob/main/pdf/NuSTAR_CalCoffee_20241002.pdf)
for additional info. This version retains veto'ing of events in the anti-coincidence
shields. It puts the output into an `event_cl_highrate` subdirectory.


[run_nupipeline_highrate_noshield.sh](run_nupipeline_highrate_noshield.sh) applies an even more lenient
electronics filtering that ignores the on-board anti-coincidence radiation
shield. This can be useful when the source rate is high enough than
chance coincidences with the events in the anit-coincidence that can be
useful for source rates in excess of a few thousand counts per second. This
can be useful, for example, to prevent artificially veto'ing events during
bright Type I X-ray bursts. It puts the output into an `event_cl_noshield` subdirectory.

