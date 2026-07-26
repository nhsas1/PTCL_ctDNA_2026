#!/bin/bash
# Submits fragment extraction in chunks of 10 array tasks, waiting for each chunk before
# submitting the next, to stay within queue limits.
#
# This is superseded by submit_chunked.sh, which is the same driver with the target script
# taken as an argument rather than hardcoded, and which now also checks that each chunk
# actually completed rather than merely leaving the queue. Prefer submit_chunked.sh.
SCRIPT=/scratch/alice/n/nhsas1/ptcl_ctdna_thesis/scripts/fragmentomics/run_fragment_extraction.slurm
CHUNK=10
START=$1
END=$2

current=${START}
while [ ${current} -le ${END} ]; do
    stop=$((current + CHUNK - 1))
    if [ ${stop} -gt ${END} ]; then
        stop=${END}
    fi

    echo "Submitting array ${current}-${stop}"
    JOBID=$(sbatch --array=${current}-${stop}%4 --parsable ${SCRIPT})

    if [ -z "${JOBID}" ]; then
        echo "ERROR: submission failed for range ${current}-${stop}, stopping"
        exit 1
    fi
    echo "Job ID: ${JOBID}"

    # poll until every task of this job has left the queue
    while squeue -j ${JOBID} -h 2>/dev/null | grep -q .; do
        sleep 15
    done
    echo "Range ${current}-${stop} complete"

    current=$((stop + 1))
done

echo "All requested ranges submitted and completed"
