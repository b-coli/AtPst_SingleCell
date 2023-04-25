#!/bin/bash
# specify path to container


#!/usr/bin/sh

TMPDIR=~/rstudio-tmp
rm -R $TMPDIR/tmp/rstudio-server

mkdir -p $TMPDIR/tmp/rstudio-server
uuidgen > $TMPDIR/tmp/rstudio-server/secure-cookie-key
chmod 600 $TMPDIR/tmp/rstudio-server/secure-cookie-key

mkdir -p $TMPDIR/var/lib
mkdir -p $TMPDIR/var/run

singularity exec -B $TMPDIR/var/lib:/var/lib/rstudio-server -B $TMPDIR/var/run:/var/run/rstudio-server -B $TMPDIR/tmp:/tmp docker://bcoli/atpst_singlecell_r:1.0.1 rserver --server-user $USER
